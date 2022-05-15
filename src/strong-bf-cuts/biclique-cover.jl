# Return a list of tuple (L, R), where (l, r) is an edge for each l in L, r in R,
# and the union of all (L, R) covers all edges
function biclique_cover(adjmatrix)
    covering = Tuple[]
    uncovered = copy(adjmatrix)
    
    while any(uncovered)
        _, j = Tuple(findfirst(uncovered))

        # Find maximal R set given L
        function find_R(L)
            local R = findall(all.(eachcol(@view adjmatrix[L, :])))
            local numcovered = count(@view uncovered[L, R])
            return R, numcovered
        end

        # Greedy algorithm to find a maximal biclique while maximizing numcovered
        L = findall(@view adjmatrix[:, j])
        R, numcovered = find_R(L)
        best = (L, R, numcovered)
        
        # Remove each element in L greedily and count the number of newly covered pairs
        while length(L) > 1
            Rl = [(l, find_R(l)...) for l in subsets(L, length(L) - 1)]     # List of (l, r, numcovered)
            _, ind = findmax(last.(Rl))
            L, R, numcovered = Rl[ind]
            # Update the best tuple if numcovered is larger, or if equal, then update if L * R is larger
            if numcovered > last(best) || (numcovered == last(best) && length(L) * length(R) > length(best[1]) * length(best[2]))
                best = (L, R, numcovered)
            end
        end

        # Push the best biclique into the list, modify uncovered
        L, R, numcovered = best
        push!(covering, (L, R))
        uncovered[L, R] .= 0
    end

    return covering
end