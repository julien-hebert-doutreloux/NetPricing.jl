# Split a node into 2 given an arc
function split(sess::BBSession, node::BBNode, arc::Int; minwidth=1f)
    # If min_width is reached, don't split
    node.tmax[arc] - node.tmin[arc] >= minwidth || return nothing, nothing

    new_min_tolls = copy(node.tmin)
    new_max_tolls = copy(node.tmax)

    midpoint = (new_min_tolls[arc] + new_max_tolls[arc]) / 2
    new_min_tolls[arc] = new_max_tolls[arc] = midpoint

    new_min_node = BBNode(sess, node.tmin, new_max_tolls)
    new_max_node = BBNode(sess, new_min_tolls, node.tmax)

    return new_min_node, new_max_node
end

# Determine if a node is prunable
function is_prunable(sess::BBSession, node::BBNode; objtol=1e-3)
    return node.bound <= sess.lowerbound + objtol
end
is_prunable(::BBSession, ::Nothing) = true

# Step a node
function step!(sess::BBSession; minwidth=1.0, objtol=1e-3, sloperatio=0.01, numcandidates=20)
    prob = sess.prob

    node = dequeue!(sess.queue)
    node.bound > sess.lowerbound + objtol || return      # Prune if UB < LB

    # Check if the node reached minwidth, then use the exact solver
    if max_width(node) < minwidth
        set_Tminmax!(sess.exact_model, node.tmin, node.tmax)
        point = lowerbound_from_exact_model(sess.exact_model)
        update_lowerbound!(sess, point)
        return node, nothing, nothing, nothing
    end

    # Choose variable
    # Filter out 10 arcs with max widths (and >= minwidth)
    ws = widths(node)
    arcs = first.(sort!(filter!(w -> w[2] >= minwidth, collect(ws)), by=last, rev=true)[1:min(numcandidates,end)])

    # Estimate the root slope and the mid point in each direction
    midrevdiff = Float64[]
    expected_midrevdiff = Float64[]

    min_point = lowerbound_from_tolls(prob, node.tmin)
    update_lowerbound!(sess, min_point)

    for a in arcs
        close_tolls = copy(node.tmin)
        mid_tolls = copy(node.tmin)
        
        close_tolls[a] += ws[a] * sloperatio
        mid_tolls[a] += ws[a] / 2

        close_point = lowerbound_from_tolls(prob, close_tolls)
        mid_point = lowerbound_from_tolls(prob, mid_tolls)

        update_lowerbound!(sess, close_point)
        update_lowerbound!(sess, mid_point)

        push!(midrevdiff, mid_point.rev - min_point.rev)
        push!(expected_midrevdiff, (close_point.rev - min_point.rev) / sloperatio / 2)
    end

    # Choose the arc with the most absolute diff in midrevdiff
    midrevdiff_error = abs.(expected_midrevdiff - midrevdiff)
    _, index = findmax(midrevdiff_error)
    branch_arc = arcs[index]

    # Branch
    left_node, right_node = split(sess, node, branch_arc, minwidth=minwidth)

    # Prune
    for child in (left_node, right_node)
        is_prunable(sess, child, objtol=objtol) || enqueue!(sess.queue, child => child.bound)
    end

    return node, branch_arc, left_node, right_node
end

## Automatic step
function solve!(sess::BBSession; objtol=1e-3, minwidth=1, log=nothing)
    isnothing(log) || empty!(log)
    while !isempty(sess.queue)
        # If LB >= UB, stop
        is_prunable(sess, peek(sess.queue)[1], objtol=objtol) && break

        parent, branch_arc, left, right = step!(sess, objtol=objtol, minwidth=minwidth)
        println(parent)
        isnothing(log) || push!(log, BBLogEntry(parent, branch_arc, left, right, sess.bestsolution))
    end
    return sess
end

solve(prob::Problem, tmax; kwargs...) = solve!(BBSession(prob, tmax); kwargs...)