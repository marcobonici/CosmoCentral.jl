function SteMDerivative(x::Vector{Float64}, y::Vector{Float64})
    x_copy = deepcopy(x)
    y_copy = deepcopy(y)
    if minimum(y) == maximum(y)
        der = 0
    else
        nonlinear = true
        while nonlinear
            coefficients = CustomRegression(x_copy, y_copy)
            y_fit = zeros(length(x_copy))
            percent_diff = zeros(length(x_copy))
            y_fit .= coefficients[1] .+ x_copy .* coefficients[2]
            percent_diff = abs.((y_copy .- y_fit) ./ y_copy)
            if all(percent_diff .<= 0.01) || length(x_copy) < 3
                nonlinear = false
                der = coefficients[2]
            else
                pop!(x_copy)
                pop!(y_copy)
                popfirst!(y_copy)
                popfirst!(x_copy)
            end
        end
    end
    return der
end
