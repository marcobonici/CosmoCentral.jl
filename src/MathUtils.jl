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
    if x <= first(Array)
        idx = 1
    elseif x >= last(Array)
        idx = length(Array)-1
    else
        while !(x >= Array[idx] && x <= Array[idx+1])
            idx += 1
        end
    end
    return idx
end

function CustomRegression(x::Vector{Float64}, y::Vector{Float64})
    x_mean = mean(x)
    y_mean = mean(y)
    beta = 0
    den = 0
    for i in 1:length(x)
        beta += (x[i]-x_mean)*(y[i]-y_mean)
        den  += (x[i]-x_mean)^2
    end
    beta /= den
    alpha = y_mean - beta* x_mean
    return alpha, beta
end
