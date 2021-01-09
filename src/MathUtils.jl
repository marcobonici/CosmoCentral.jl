"""
    LogSpaced(min::Float64, max::Float64, n::Int64)

This function evaluates ``n`` points, logarithmically spaced between
    ``min`` and ``max``.
"""
function LogSpaced(min::Float64, max::Float64, n::Int64)
    logmin = log10(min)
    logmax = log10(max)
    logarray = Array(LinRange(logmin, logmax, n))
    return exp10.(logarray)
end


"""
    BinSearch(x::Float64, Array::Vector{Float64})

Given a value ``z`` and an Array, determines the couple of array elements where
``z`` lies and returns the index correspondin to the first value.
"""
function BinSearch(x::Float64, Array::Vector{Float64})
    idx = 1
    while !(x >= Array[idx] && x <= Array[idx+1])
        idx += 1
    end
    return idx
end
