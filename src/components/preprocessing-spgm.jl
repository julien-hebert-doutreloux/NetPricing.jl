Base.reverse(graph::SimpleWeightedDiGraph) = SimpleWeightedDiGraph(sparse(transpose(graph.weights)))

function preprocess_spgm(prob::Problem, k; graph = build_graph(prob))
    comm = prob.K[k]
    orig, dest = comm.orig, comm.dest

    # Lower bound distances
    reset_tolls!(graph, prob)
    revgraph = reverse(graph)
    lower_orig = dijkstra_shortest_paths(graph, orig).dists
    lower_dest = dijkstra_shortest_paths(revgraph, dest).dists

    # Upper bound distances
    disable_tolls!(graph, prob)
    revgraph = reverse(graph)
    upper_orig = dijkstra_shortest_paths(graph, orig).dists
    upper_dest = dijkstra_shortest_paths(revgraph, dest).dists

    # Bounds from orig to dest
    lower_total = lower_orig[dest]
    upper_total = upper_orig[dest]

    # Maps
    V = prob.V
    A = ProblemArc[]
    Amap = Int[]
    srcs = Dict{Int,BitSet}()
    dsts = Dict{Int,BitSet}()

    # Add all tolled arcs
    for (idx, arc) in enumerate(prob.A)
        arc.toll || continue
        i, j = arc.src, arc.dst
        cost = arc.cost

        # Elimination rules
        # Rules #1, #2: lower and upper prices are the same => arc has no effect
        lower_orig[j] >= upper_orig[j] && continue
        lower_dest[i] >= upper_dest[i] && continue

		# O-D (rule #5), I-D, O-J (enhanced rules #1, #2)
        upper_total <= lower_orig[i] + cost + lower_dest[j] && continue
        upper_dest[i] <= cost + lower_dest[j] && continue
        upper_orig[j] <= lower_orig[i] + cost && continue

        # Pass, add to graph
        push!(A, arc)
        push!(Amap, idx)

        haskey(srcs, i) || (srcs[i] = BitSet())
        push!(srcs[i], j)

        haskey(dsts, j) || (dsts[j] = BitSet())
        push!(dsts[j], i)
    end

    arcsset = Set((a.src, a.dst) for a in A)

    # Add toll-free arcs helper
    function add_tollfree(src, dst, cost)
        # Perturb the cost a little, so any bilevel feasible path is unique
        cost += rand() * 1e-6

        # If there already exists an arc from src to dst, create a virtual node
        if (src, dst) in arcsset
            V += 1
            push!(A, ProblemArc(V, dst, 0, false))
            push!(Amap, 0)
            push!(arcsset, (V, dst))
            dst = V
        end

        push!(A, ProblemArc(src, dst, cost, false))
        push!(Amap, 0)
        push!(arcsset, (src, dst))
    end

    # O-D arc
    add_tollfree(orig, dest, upper_total)

    # O-srcs arcs
    for i in keys(srcs)
        cost = upper_orig[i]
        isinf(cost) && continue
        
        # Rule #7
        upper_total <= cost + lower_dest[i] && continue

        add_tollfree(orig, i, cost)
    end

    # dsts-D arcs
    for j in keys(dsts)
        cost = upper_dest[j]
        isinf(cost) && continue

        # Rule #8
        upper_total <= lower_orig[j] + cost && continue

        add_tollfree(j, dest, cost)
    end

    # dsts-srcs arcs
    for j in keys(dsts)
        upper_j = dijkstra_shortest_paths(graph, j).dists
        for i in keys(srcs)
            cost = upper_j[i]

            # Self-loop or disconnect pair
            i == j && continue
            isinf(cost) && continue

            # Rules #1, #2
            lower_orig[i] >= upper_orig[i] && continue
            lower_dest[j] >= upper_dest[j] && continue

			# O-D (rule #6), J-D (rule #3), O-I (rule #4)
			upper_total <= lower_orig[j] + cost + lower_dest[i] && continue;
			upper_dest[j] <= cost + lower_dest[i] && continue;
			upper_orig[i] <= lower_orig[j] + cost && continue;

            # Pass, add to graph
            if i in dsts[j]
				# i connects to j by a tolled arc
				# Only connect j to i if j has another i' and i has another j'
                if length(dsts[j]) > 1 && length(srcs[i]) > 1
                    add_tollfree(j, i, cost)
                end
            else
                add_tollfree(j, i, cost)
            end
        end
    end

    reset_tolls!(graph, prob)

    # Prune the nodes
    hasinput = BitSet()
    hasoutput = BitSet()
    for a in A
        push!(hasoutput, a.src)
        push!(hasinput, a.dst)
    end
    used_nodes = BitSet(vcat(orig, dest, (hasinput ∩ hasoutput)...))
    Vmap = collect(used_nodes)

    # Translate and filter the arcs
    Vrevmap = revmap(Vmap, Vmap[end])
    map!(arc -> ProblemArc(Vrevmap[arc.src], Vrevmap[arc.dst], arc.cost, arc.toll), A, A)
    filteridx = findall(arc -> arc.src > 0 && arc.dst > 0, A)
    A = A[filteridx]
    Amap = Amap[filteridx]

    Vrevmap = revmap(Vmap, prob.V)
    Arevmap = revmap(Amap, length(prob.A))
    return PreprocessedProblem(prob, length(Vmap), A, Vmap, Amap, Vrevmap, Arevmap, k,
        Vrevmap[prob.K[k].orig],
        Vrevmap[prob.K[k].dest])
end
