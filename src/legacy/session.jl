struct BBNode
    id::Int
    tmin::TollsDict
    tmax::TollsDict
    bound::Float64
end

mutable struct BBSession
    prob::Problem
    graph::AbstractGraph
    upperbound_model::Model
    exact_model::Model

    queue::PriorityQueue{BBNode,Float64}
    lowerbound::Float64
    bestsolution::FeasiblePoint
    nodecounter::Int
end

## Create a node from a min and a max point
function BBNode(sess::BBSession, tmin::AbstractTollsDict, tmax::AbstractTollsDict)
    prob = sess.prob

    # Calculate the bound
    model = sess.upperbound_model
    set_Tminmax!(model, tmin, tmax)
    optimize!(model)
    bound = objective_value(model)

    sess.nodecounter += 1
    return BBNode(sess.nodecounter, tmin, tmax, bound)
end

## Branch-and-bound session
function BBSession(prob::Problem, tmax::AbstractTollsDict)
    graph = build_graph(prob)
    ub_model = make_upperbound_model(prob)
    exact_model = make_exact_model(prob)
    queue = PriorityQueue{BBNode,Float64}(Base.Order.Reverse)
    trivial_sol = lowerbound_from_tolls(prob, zero_tolls(prob), graph=graph)

    sess = BBSession(prob, graph, ub_model, exact_model, queue, trivial_sol.rev, trivial_sol, 0)

    root = BBNode(sess, zero_tolls(prob), tmax)
    enqueue!(queue, root => root.bound)

    return sess
end

## Pretty print
function Base.show(io::IO, node::BBNode)
    ws = widths(node)
    minwidth = round(minimum(values(ws)), digits=2)
    maxwidth = round(maximum(values(ws)), digits=2)
    print(io, "Node $(node.id) with bound $(round(node.bound, digits=3)), width $minwidth ~ $maxwidth")
end

function Base.show(io::IO, sess::BBSession)
    lb = round(sess.lowerbound, digits=3)
    ub = round(peek(sess.queue)[2], digits=3)
    print(io, "Session with {lb = $lb, ub = $ub, queue = $(length(sess.queue))}")
end

## Widths of a node
width(node::BBNode, i) = node.tmax[i] - node.tmin[i]
widths(node::BBNode) = Dict(i => width(node, i) for i in keys(node.tmin))
max_width(node::BBNode) = maximum(values(widths(node)))

## Update session lower bound
function update_lowerbound!(sess::BBSession, point::FeasiblePoint)
    if point.rev > sess.lowerbound
        sess.lowerbound = point.rev
        sess.bestsolution = point
    end
end