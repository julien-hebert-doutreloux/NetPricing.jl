function run_experiment(name, probfile, maxpaths, time; outfile=nothing)
    prob = read_problem(probfile)
    func = Core.eval(NetPricing, Symbol(name))

    println("===== PILOT RUN START =====")
    pilotmodel, _ = func(prob, 10, 60)
    optimize!(pilotmodel)
    println("===== PILOT RUN END =====")
    println()

    println("===== MAIN RUN START =====")
    println("===== PREPROCESS =====")
    model, preproc_time = func(prob, maxpaths, time)
    println("===== OPTIMIZE =====")
    optimize!(model)
    println("===== MAIN RUN END =====")
    println()

    totaltime = solve_time(model) + preproc_time
    obj = result_count(model) > 0 ? objective_value(model) : 0
    bound = objective_bound(model)
    gap = relative_gap(model)
    step = node_count(model)

    if !isnothing(outfile)
        open(outfile, "w") do f
            JSON.print(f, Dict(
                :time => totaltime,
                :obj => obj,
                :bound => bound,
                :gap => gap,
                :step => step,
                :preproc_time => preproc_time
                ))
        end
    end

    println("===== SUMMARY =====")
    println("TIME: $totaltime")
    println("OBJ: $obj")
    println("BOUND: $bound")
    println("GAP: $gap")
    println("STEP: $step")
    println("PREPROC TIME: $preproc_time")
end

function adaptive_experiment(prob, maxpaths, time; nonadaptive, preprocessed)
    pprob = nothing
    preprocess_time = preprocessed ? (@elapsed begin
        pprob = preprocess(prob, maxpaths=maxpaths, noempty=true, nospgm=true)
    end) : 0.0
    return adaptive_model(prob, nonadaptive=nonadaptive, preprocess=pprob, timelimit=(time - preprocess_time), threads=1), preprocess_time
end

adaptive_unprocessed(args...) = adaptive_experiment(args..., nonadaptive=false, preprocessed=false)
adaptive_preprocessed(args...) = adaptive_experiment(args..., nonadaptive=false, preprocessed=true)
nonadaptive_unprocessed(args...) = adaptive_experiment(args..., nonadaptive=true, preprocessed=false)
nonadaptive_preprocessed(args...) = adaptive_experiment(args..., nonadaptive=true, preprocessed=true)
