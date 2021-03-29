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
    LnSpaced(min::Float64, max::Float64, n::Int64)

This function evaluates ``n`` points, natural logarithmically spaced between
    ``min`` and ``max``.
"""
function LnSpaced(min::Float64, max::Float64, n::Int64)
    lnmin = log(min)
    lnmax = log(max)
    lnarray = Array(LinRange(lnmin, lnmax, n))
    return exp.(lnarray)
end

"""
    BinSearch(x::Float64, Array::Vector{Float64})

Given a value ``z`` and an Array, determines the couple of array elements where
``z`` lies and returns the index corresponding to the first value.
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

"""
    CustomRegression(x::Vector{Float64}, y::Vector{Float64})

Given two arrays ``x`` and ``y``, performs the linear regression and returns the
coefficients ``c`` and ``m`` of the fitted line.
"""
function CustomRegression(x::Vector{Float64}, y::Vector{Float64})
    x_mean = mean(x)
    y_mean = mean(y)
    m = 0
    den = 0
    for i in 1:length(x)
        m += (x[i]-x_mean)*(y[i]-y_mean)
        den  += (x[i]-x_mean)^2
    end
    m /= den
    c = y_mean - m* x_mean
    return c, m
end

"""
    SimpsonWeightArray(n::Int64)

This function evaluates an array with the Simpson weight, for Simpson
Integration
"""
function SimpsonWeightArray(n::Int64)
    number_intervals = floor((n-1)/2)
    weight_array = zeros(n)
    if n == number_intervals*2+1
        for i in 1:number_intervals
            weight_array[Int((i-1)*2+1)] += 1/3
            weight_array[Int((i-1)*2+2)] += 4/3
            weight_array[Int((i-1)*2+3)] += 1/3
        end
    else
        weight_array[1] += 0.5
        weight_array[2] += 0.5
        for i in 1:number_intervals
            weight_array[Int((i-1)*2+1)+1] += 1/3
            weight_array[Int((i-1)*2+2)+1] += 4/3
            weight_array[Int((i-1)*2+3)+1] += 1/3
        end
        weight_array[length(weight_array)]   += 0.5
        weight_array[length(weight_array)-1] += 0.5
        for i in 1:number_intervals
            weight_array[Int((i-1)*2+1)] += 1/3
            weight_array[Int((i-1)*2+2)] += 4/3
            weight_array[Int((i-1)*2+3)] += 1/3
        end
        weight_array ./= 2
    end
    return weight_array
end

"""
    SimpsonWeightMatrix(n::Int64)

This function evaluates an array with the Simpson weight, for Simpson
Integration of the Lensing Efficiency
"""
function SimpsonWeightMatrix(n::Int64)
    number_intervals = floor((n-1)/2)
    weight_matrix = zeros(n, n)
    for i in 1:n
        weight_matrix[i,i:n] = SimpsonWeightArray(n-i+1)
    end
    return weight_matrix
end

"""
    Difference(InputArray::Vector{Float64})

This function evaluates the n-th discrete difference of a given 1-D array.
"""
function Difference(InputArray::Vector{Float64})
    OutputArray = zeros(length(InputArray))
    for i in range 1:length(InputArray)-1
        OutputArray[i] = InputArray[i+1] - InputArray[i]
    end
    return OutputArray
end
