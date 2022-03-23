struct ConsecutivePairIterator{T<:AbstractVector}
    array::T
end

consecutive_pairs(array) = ConsecutivePairIterator(array)

function Base.iterate(iter::ConsecutivePairIterator, i = 1)
    array = iter.array
    return i < length(array) ? ((array[i],array[i + 1]), i + 1) : nothing
end

function Base.eltype(::Type{ConsecutivePairIterator{T}}) where {T}
    R = eltype(T) 
    return Tuple{R,R}
end
Base.length(iter::ConsecutivePairIterator) = length(iter.array) - 1