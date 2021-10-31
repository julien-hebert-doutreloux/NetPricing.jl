## Plot 2D
function plot2d(prob::Problem, range_x, range_y=range_x; log=nothing, func=revenue)
    length(tolled_arcs(prob)) == 2 || throw(ArgumentError("Problem must have exactly 2 tolled arcs"))
    arc_x =  tolled_arcs(prob)[1]
    arc_y =  tolled_arcs(prob)[2]
    graph = build_graph(prob)
    nk = length(prob.K)

    all_revs = [func(prob, Dict(arc_x => toll_x, arc_y => toll_y), graph=graph)
        for toll_x in range_x, toll_y in range_y]

    fig = Figure()

    # Main axis
    main_axis = Axis(fig[1,1], aspect=DataAspect())
    pl = heatmap!(main_axis, range_x, range_y, all_revs)
    Colorbar(fig[1, 2], pl)

    # Log
    if !isnothing(log)
        plog = process_log(log)

        rects = Node([Rect(minimum(range_x), minimum(range_y), maximum(range_x), maximum(range_y))])

        transparent = colorant"transparent"
        prunedcolor = RGBA{Float32}(1,0,0,0.4)
        colors = Node(Colorant[transparent])

        optimal_point = Node(Point2f0[])

        poly!(main_axis, rects, color=colors, strokecolor=:black, strokewidth=1)
        scatter!(main_axis, optimal_point, color = :red)

        slider = Slider(fig[2, 1], range=1:length(plog), startvalue=1)
        on(slider.value) do step
            entry = plog[step]
            rects.val = [Rect(node.tmin[1], node.tmin[2], node.tmax[1] - node.tmin[1], node.tmax[2] - node.tmin[2]) for node in entry.all_nodes]
            colors.val = [node.id in entry.pruned ? prunedcolor : transparent for node in entry.all_nodes]
            bestsol = entry.original.bestsolution.tolls
            optimal_point[] = [Point2f0(bestsol[1], bestsol[2])]
            try
                notify.((colors,rects))
            catch
                notify.((rects,colors))
            end
        end
    end

    current_figure()
end

## Process log
struct ProcessedLogEntry
    original::BBLogEntry
    all_nodes::Vector{BBNode}
    pruned::Vector{Int}
end

function process_log(log)
    processed = ProcessedLogEntry[]
    all_nodes = BBNode[log[1].parent]
    skipped = Int[]

    for entry in log
        # Check if there're children, replace the parent with its children
        if !isnothing(entry.left)
            filter!(n -> n.id != entry.parent.id, all_nodes)
            push!(all_nodes, entry.left)
            push!(all_nodes, entry.right)
        # Otherwise, add it to skipped nodes
        else
            push!(skipped, entry.parent.id)
        end
        # Recheck the lower bound, eliminate the nodes with smaller bound
        lowerbound = entry.bestsolution.rev
        pruned = [node.id for node in all_nodes if node.bound < lowerbound + 1e-3]
        push!(processed, ProcessedLogEntry(entry, copy(all_nodes), union(pruned, skipped)))
    end

    return processed
end