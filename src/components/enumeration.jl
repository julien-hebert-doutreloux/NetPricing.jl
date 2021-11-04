struct BFEnumPathInfo
    path::Vector{Int}
    cost::Float64
    spurstart::Int
    excluded::Vector{Int}
end

function get_tolled_list(path, tolledindices)
    tolledlist = Int[]
    for i in 1:(length(path)-1)
        arc = (path[i], path[i+1])
        if haskey(tolledindices, arc)
            push!(tolledlist, tolledindices[arc])
        end
    end
    return tolledlist
end

function enumerate_bilevel_feasible(graph::AbstractGraph, orig, dest, prob, numpaths;
    purge=true,     # Purge all non-bilevel-feasible paths (quadratic time)
    prune1=true)    # Prune rule 1: if the spur part is toll-free, any subsequent path is dominated if it is longer

    reset!(graph, prob)
    init_ws = copy(graph.weights)

    # First shortest path
    spath, scost = shortest_path(graph, orig, dest)
    queue = PriorityQueue(BFEnumPathInfo(spath, scost, orig, []) => scost)
    bfpaths = Vector{Int}[]

    tolledindices = tolled_srcdst_to_index(prob)
    arccosts = srcdst_to_cost(prob)

    purgelist = BitSet[]

    while length(bfpaths) < numpaths && !isempty(queue)
        pathinfo = dequeue!(queue)
        path = pathinfo.path
        spurstart = pathinfo.spurstart
        excluded = pathinfo.excluded

        tolledlist = get_tolled_list(path, tolledindices)

        # Purging non-bilevel-feasible paths
        if purge && !isempty(tolledlist)
            tolledset = BitSet(tolledlist)
            # If no prior set is a subset of this set, it is non-dominated
            if all(set -> set ⊈ tolledset, Iterators.reverse(purgelist))
                push!(bfpaths, pathinfo.path)
                push!(purgelist, tolledset)
            end
        else
            push!(bfpaths, pathinfo.path)
        end

        # If there's no tolled arc, break
        isempty(tolledlist) && break

        # Get the list of tolled arcs between spur nodes
        spurstartidx = something(findfirst(a -> prob.A[a].dst == spurstart, tolledlist), 0) + 1

        tolledlist = tolledlist[spurstartidx:end]

        # Prune 1 upperbound
        upperbound = Inf
        spurnode = spurstart

        for a in tolledlist
            spurnodeidx = findfirst(isequal(spurnode), path)
            rootpart = path[1:(spurnodeidx-1)]
            newexcluded = [excluded; a]

            disable_nodes!(graph, rootpart)
            disable_arcs!(graph, prob, newexcluded)

            spurpart, spurcost = shortest_path(graph, spurnode, dest)
            
            restore!(graph, init_ws)

            if spurcost != Inf
                # Add the cost of the root part and end part
                rootpartcost = get_path_cost([rootpart; spurnode], arccosts)
                spurcost += rootpartcost

                # If prune 1 and the path is longer than the upperbound, skip
                prune1 && spurcost >= upperbound && continue

                newpath = [rootpart; spurpart]

                enqueue!(queue, BFEnumPathInfo(newpath, spurcost, spurnode, newexcluded), spurcost)

                # Prune rule 1: if the spurpart is toll-free, change the upperbound
                if prune1
                    spurtolledlist = get_tolled_list(spurpart, tolledindices)
                    isempty(spurtolledlist) && (upperbound = min(upperbound, spurcost))
                end
            end

            spurnode = prob.A[a].dst
        end
    end

    reset!(graph, prob)

    return bfpaths
end