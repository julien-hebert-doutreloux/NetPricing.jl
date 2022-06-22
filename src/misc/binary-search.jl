function binary_root_search(f, a, b, fa, fb)
    m = (a + b) / 2
    fm = f(m)
    (fm == 0) && return m
    left = (fm * fa < 0)
    return binary_root_search(f, left ? a : b, m, left ? fa : fb, fm)
end

function binary_root_search(f, a, b)
    fa, fb = f(a), f(b)
    (fa == 0) && return a
    (fb == 0) && return b
    binary_root_search(f, a, b, f(a), f(b))
end
