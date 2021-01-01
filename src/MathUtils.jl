function LogSpaced(minarray::Float64, maxarray::Float64, n::Int64)
    array = zeros(n)
    logmin = log10(minarray)
    logmax = log10(maxarray)
    logarray = Array(LinRange(logmin, logmax, n))
    for (idx, mylogarray) in enumerate(logarray)
        array[idx] = 10^mylogarray
    end
    return array
end
