## Build model
function make_exact_model(prob::Problem; sdtol=1e-10, silent=true, bigM=1000)
    model = Model(() -> Gurobi.Optimizer(current_env()))
    set_optimizer_attribute(model, MOI.Silent(), silent)
    
    nv = prob.V
    na = length(prob.A)
    nk = length(prob.K)
    a1 = tolled_arcs(prob)
    model[:prob] = prob
    model[:a1] = a1

    c = cost_vector(prob)
    A = incidence_matrix(prob)
    b = [sourcesink_vector(prob, k) for k in 1:nk]
    d = demand_vector(prob)

    @variables(model,
        begin
            0 ≤ x[k=1:nk, a=1:na] ≤ 1, (binary = a in a1)   # All arcs
            λ[k=1:nk,i=1:nv]                                # Duals
            t[a=a1] ≥ 0                                     # Tolls
            tx[k=1:nk,a=a1] ≥ 0                             # Bilinear
        end)

    # Expressions
    tfull = Array{Any}(zeros(na))
    for a in a1
        tfull[a] = t[a]
    end
    
    primalobj = model[:primalobj] = [c' * x[k,:] + sum(tx[k,:]) for k in 1:nk]
    dualobj = model[:dualobj] = [b[k]' * λ[k,:] for k in 1:nk]

    @objective(model, Max, sum(d .* tx))

    @constraints(model,
        begin
            primalfeas[k=1:nk], A * x[k,:] .== b[k]
            dualfeas[k=1:nk], A' * λ[k,:] .≤ c + tfull
            strongdual[k=1:nk], primalobj[k] ≤ dualobj[k] + sdtol
            bilinear1[k=1:nk], tx[k,:] .≥ 0
            bilinear2[k=1:nk], tx[k,:] .≤ bigM .* x[k,a1]
            bilinear3[k=1:nk], t .- tx[k,:] .≥ 0
            bilinear4[k=1:nk], t .- tx[k,:] .≤ bigM .* (1 .- x[k,a1])
        end)
    
    return model
end

## Set Tmin-max
function set_Tminmax!(model, tmindict, tmaxdict)
    nk = length(model[:prob].K)
    a1 = model[:a1]
    tmin = [tmindict[a] for a in a1]
    tmax = [tmaxdict[a] for a in a1]

    bilinear1, bilinear2, bilinear3, bilinear4 = (model[name] for name in [:bilinear1, :bilinear2, :bilinear3, :bilinear4])

    for k in 1:nk
        x = model[:x][k,a1]
        set_normalized_coefficient.(bilinear1[k], x, -tmin)
        set_normalized_coefficient.(bilinear2[k], x, -tmax)
        set_normalized_coefficient.(bilinear3[k], x, tmin)
        set_normalized_coefficient.(bilinear4[k], x, tmax)
        set_normalized_rhs.(bilinear3[k], tmin)
        set_normalized_rhs.(bilinear4[k], tmax)
    end
    set_lower_bound.(model[:t], tmin)
    set_upper_bound.(model[:t], tmax)
    return model
end
