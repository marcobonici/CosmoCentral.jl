function LogSpaced(minarray::Float64, maxarray::Float64, n::Int64)
    logmin = log10(minarray)
    logmax = log10(maxarray)
    logarray = Array(LinRange(logmin, logmax, n))
    array = exp10.(array)
    return array
end
