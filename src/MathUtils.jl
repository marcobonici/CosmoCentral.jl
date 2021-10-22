"""
    LogSpaced(min::Float64, max::Float64, n::Int64)

This function evaluates ``n`` points, logarithmically spaced between
    ``min`` and ``max``.
"""
function LogSpaced(min::T, max::T, n::I) where {T,I}
    logmin = log10(min)
    logmax = log10(max)
    logarray = Array(LinRange(logmin, logmax, n))
    return exp10.(logarray)
end

"""
    ΔLogSpaced(min::Float64, max::Float64, n::Int64)

This function evaluates ``n`` points, logarithmically spaced between
    ``min`` and ``max``.
"""
function ΔLogSpaced(min::T, max::T, n::I) where {T,I}
    λmin = log10(min)
    λmax = log10(max)
    Δℓ = zeros(n)
    for i in 1:n
        Δℓ[i] = 10^(λmin+i*(λmax-λmin)/n)-10^(λmin+(i-1)*(λmax-λmin)/n)
    end
    return Δℓ
end

"""
    BinSearch(x::Float64, Array::Vector{Float64})

Given a value ``z`` and an Array, determines the couple of array elements where
``z`` lies and returns the index corresponding to the first value.
"""
function BinSearch(x::T, Array::Vector{T}) where T
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
function CustomRegression(x::Vector{T}, y::Vector{T}) where T
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
function SimpsonWeightArray(n::I) where I
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
function SimpsonWeightMatrix(n::I) where I
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
function Difference(InputArray::Vector{T}) where T
    OutputArray = zeros(length(InputArray)-1)
    for i in 1:length(InputArray)-1
        OutputArray[i] = InputArray[i+1] - InputArray[i]
    end
    return OutputArray
end


"""
    UnevenTrapzWeightArray(InputArray::Vector{Float64})

This function evaluates the array of weights for the integration with uneven
spaced points and the Trapezoidal rule.
"""
function UnevenTrapzWeightArray(InputArray::Vector{T}) where T
    WeightArray = zeros(length(InputArray))
    DifferenceArray = Difference(InputArray)
    for idx in 2:length(WeightArray)-1
        WeightArray[idx] = (DifferenceArray[idx]+DifferenceArray[idx-1])/2
    end
    WeightArray[1] = DifferenceArray[1]/2
    WeightArray[length(InputArray)] = DifferenceArray[length(InputArray)-1]/2
    return WeightArray
end

function UnevenTrapzWeightMatrix(InputMatrix::Matrix{T}) where T
    WeightMatrix = zeros(size(InputMatrix))
    for lidx in 1:length(WeightMatrix[:,1])
        WeightMatrix[lidx,:] = UnevenTrapzWeightArray(InputMatrix[lidx,:])
    end
    return WeightMatrix
end

function uᵢⱼ(n::I, i::I, j::I) where I
    u = zeros(floor(Int,0.5*n*(n+1)))
    u[floor(Int64, (j-1)*n+i-0.5*j*(j-1))] = 1
    return u
end

function eᵢ(n::I, i::I) where I
    e = zeros(n)
    e[i] = 1
    return e
end

function Eᵢⱼ(n::I, i::I, j::I) where I
    E = zeros((n,n))
    E[i,j] = 1
    return E
end

function Tᵢⱼ(n::I, i::I, j::I) where I
    if i != j
        T = Eᵢⱼ(n,i,j) + Eᵢⱼ(n,j,i)
    else
        T = Eᵢⱼ(n,i,i)
    end
    return T
end

function EliminationMatrix(n::I) where I
    L = zeros((floor(Int64,0.5*n*(n+1)), n*n))
    for i in 1:n
        for j in 1:i
            L += kron(uᵢⱼ(n,i,j), transpose(vec(Eᵢⱼ(n,i,j))))
        end
    end
    return L
end

"""
    DuplicationMatrix(n::Int)

Duplication matrix ``\\boldsymbol{D}_n``.
"""
function DuplicationMatrix(n::I) where I
    D = zeros((n*n, floor(Int64,0.5*n*(n+1))))
    for i in 1:n
        for j in 1:i
            D += transpose(kron(uᵢⱼ(n,i,j), transpose(vec(Tᵢⱼ(n,i,j)))))
        end
    end
    return D
end