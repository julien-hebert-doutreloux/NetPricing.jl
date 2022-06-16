@with_kw struct ProblemGenerationArgs
    cost_dist = Uniform(5.0, 35.0)
    max_cost::Float64 = 35.0
    max_cost_prob::Float64 = 0.2
    symmetric_cost::Bool = true
    demand_dist = Uniform(1.0, 100.0)
    tolled_proportion::Float64 = 0.2
    symmetric_tolled::Bool = true
    random_tolled_proportion::Float64 = 1/3
end