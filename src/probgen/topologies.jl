grid_graph(width, height) = grid((width, height))

function delaunay_graph(num_nodes)
    return _delaunay_graph_with_coors(num_nodes)[1]
end

function voronoi_graph(num_nodes)
    num_seeds = num_nodes ÷ 2 - 1

    tess = DelaunayTessellation(num_seeds)
    width = max_coord - min_coord
    seeds = Point2D[Point(min_coord + rand() * width, min_coord + rand() * width) for _ in 1:num_seeds]
    push!(tess, seeds)

    graph = SimpleGraph(num_nodes)
    vertices = Dict()
    for edge in voronoiedges(tess)
        a, b = geta(edge), getb(edge)
        ai = get!(vertices, a, length(vertices) + 1)
        bi = get!(vertices, b, length(vertices) + 1)
        add_edge!(graph, ai, bi)
    end
    return graph
end

function _delaunay_graph_with_coors(num_nodes)
    width = max_coord - min_coord
    seeds = Point2D[Point(min_coord + rand() * width, min_coord + rand() * width) for _ in 1:num_nodes]
    sort!(seeds, by=getx)

    return _delaunay_graph_from_seeds(seeds)
end

function _delaunay_graph_from_seeds(seeds)
    tess = DelaunayTessellation(length(seeds))
    push!(tess, seeds)

    seeds_dict = Dict(s => i for (i, s) in enumerate(seeds))
    graph = SimpleGraph(length(seeds))
    for edge in delaunayedges(tess)
        a, b = geta(edge), getb(edge)
        add_edge!(graph, seeds_dict[a], seeds_dict[b])
    end
    return graph, seeds
end