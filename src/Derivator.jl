"""
    SteMDerivative(x::Vector{Float64}, y::Vector{Float64})

This function evaluates the numerical derivative according to the SteM
algorithm, [Camera et al. 2017](https://arxiv.org/abs/1606.03451):

- a linear regression over ``x`` and ``y`` is performed

- if the points obtained with the fit are close enough (less than 0.01 relative
difference) the linear ansatz is satisfied and the slope gives the derivative

- if the linear ansatz is not satisfied, the external couple of points and the
linear regression is performed again till the linear ansatz is satisfied
"""
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
