macro makerefs(model, refs)
    model = esc(model)
    ex = Expr(:block)
    ex.args = [:($model[$(Meta.quot(arg))] = []) for arg in eval(refs)]
    return ex
end

macro pushrefs(model, refs)
    model = esc(model)
    ex = Expr(:block)
    ex.args = [:(push!($model[$(Meta.quot(arg))], $(esc(arg)))) for arg in eval(refs)]
    return ex
end

macro pushemptyrefs(model, refs)
    model = esc(model)
    ex = Expr(:block)
    ex.args = [:(push!($model[$(Meta.quot(arg))], [])) for arg in eval(refs)]
    return ex
end
