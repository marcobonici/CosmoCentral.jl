function SteMDerivative(x::Vector{Float64}, y::Vector{Float64})
    coefficients = CustomRegression(x, y)
    return coefficients[2]
end
