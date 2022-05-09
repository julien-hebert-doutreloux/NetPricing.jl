# Standard basis vector
basis(i, n) = [i == j for j in 1:n]

# Count the number of commodities using each tolled arc
function count_tolled_arcs(prob::Problem, paths, weights=ones(length(paths)))
    a1 = tolled_arcs(prob)

    # Extract the set of tolled arcs for each path
    arcdict = srcdst_to_index(prob)
    a1set = BitSet(a1)
    tolled_sets = [BitSet(path_tolled_arcs(path, arcdict, a1set)) for path in paths]

    counts = Dict(a => sum(a in set ? weights[k] : 0 for (k, set) in enumerate(tolled_sets)) for a in a1)
    return counts
end
