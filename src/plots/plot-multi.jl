## Plot 2D from chosen arcs
function plot2d_multi(prob::Problem, base_tolls, arc_x, range_x, arc_y, range_y; func=revenue)
    graph = [build_graph(prob) for t in 1:Threads.nthreads()]

    all_revs = zeros(length(range_x), length(range_y))
    p = Progress(length(all_revs))

    # display(collect(Iterators.product(1:length(range_x), 1:length(range_y))))
    Threads.@threads for (i, j) in collect(Iterators.product(1:length(range_x), 1:length(range_y)))
        all_revs[i, j] = func(prob, Dict(base_tolls..., arc_x => range_x[i], arc_y => range_y[j]), graph=graph[Threads.threadid()])
        next!(p)
    end

    fig = Figure()

    # Main axis
    main_axis = Axis(fig[1,1], xlabel="Arc $(arc_x)", ylabel="Arc $(arc_y)")
    pl = heatmap!(main_axis, range_x, range_y, all_revs)
    # Colorbar(fig[1, 2], pl)

    current_figure()
end

function plot2d_multi_random(prob::Problem, base_tolls = zero_tolls(prob); x=nothing, y=nothing, func=revenue, resolution=100)
    a1 = tolled_arcs(prob)
    _, bigN = calculate_bigM(prob);

    inds = rand(1:length(a1), 2)

    !isnothing(x) && (inds[1] = findfirst(a -> a == x, a1))
    !isnothing(y) && (inds[2] = findfirst(a -> a == y, a1))

    arcs = a1[inds]
    maxes = bigN[inds]

    ranges = range.(0, maxes, length = resolution + 1)

    if !(base_tolls isa Dict)
        base_tolls = Dict(zip(a1, base_tolls))
    end

    fig = plot2d_multi(prob, base_tolls, arcs[1], ranges[1], arcs[2], ranges[2]; func=func)    
    fig
end

function write_tolls(filename, tolls)
    open(filename, "w") do io
        for t in tolls
            println(io, t)
        end
    end
end

function read_tolls(filename)
    tolls = [parse(Float64, line) for line in readlines(filename)]
    return tolls
end
