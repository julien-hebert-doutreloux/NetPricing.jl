# A global linearization method: represent sumk xk as bit string
struct BinaryLinearization <: AbstractLinearization end

function linearize!(model::Model, ::BinaryLinearization, forms, Ms, N; sdtol=1e-10)
    prob = problem(first(forms))
    a1 = tolled_arcs(prob)

    # Count how many commodities use an arc
    arc_vars = Dict(a => VariableRef[] for a in a1)
    arc_ks = Dict(a => Int[] for a in a1)

    for (k, form) in enumerate(forms)
        kprob = problem(primal(form))
        Amap = arcmap(kprob)
        x = primal(form).x

        for a in tolled_arcs(kprob)
            push!(arc_vars[Amap[a]], x[a])
            push!(arc_ks[Amap[a]], k)
        end
    end

    wmax = Dict(a => length(vars) for (a, vars) in arc_vars)            # Max num of commodities that use each arc
    numbits = Dict(a => ceil(Int, log2(w + 1)) for (a, w) in wmax)      # Number of bits to represent 0 to wmax

    # Extract bigM for each bit from Ms
    Mar = Dict{Int,Vector{Float64}}()
    for (i, a) in enumerate(a1)
        Mak = Float64[ Ms[k][i] for k in arc_ks[a] ]
        sort!(Mak, rev=true)
        Mar[a] = [Mak[2^(r-1)] for r in 1:numbits[a]]
    end
    Na = Dict(a1 .=> N)

    # Variables
    @variable(model, 0 ≤ v[a=a1,r=1:numbits[a]] ≤ 1, Bin)
    @variable(model, tv[a=a1,r=1:numbits[a]] ≥ 0)
    t = model[:t]

    # Transform x to v, tx to tv
    @constraint(model, [a=a1], sum(arc_vars[a]) == sum(2^(r-1) * v[a,r] for r=1:numbits[a]))
    alltv = flatten([2^(r-1) * tv[a,r] for r in 1:numbits[a]] for a in a1)
    sumtv = sum(alltv)

    # Linearization
    @constraint(model, [a=a1,r=1:numbits[a]], tv[a,r] ≤ Mar[a][r] * v[a,r])
    @constraint(model, [a=a1,r=1:numbits[a]], t[a] - tv[a,r] ≥ 0)
    @constraint(model, [a=a1,r=1:numbits[a]], t[a] - tv[a,r] ≤ Na[a] * (1 - v[a,r]))

    # Strong duality
    @constraint(model, sumtv <= sum(unnormalized_objective_term.(forms)) + sdtol)

    # Add the McCormick envelop for each commodity
    linearize!(model, EnvelopOnly(), forms, Ms, N; sdtol=sdtol)

    return
end

