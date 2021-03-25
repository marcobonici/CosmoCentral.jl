function CWindow(N::Vector{Float64}, NCut::Int64)
    NRight = last(N) - NCut
    NR = filter(x->x>=NRight, N)
    ThetaRight = (last(N).-NR) ./ (last(N) - NRight - 1)
    W = ones(length(N))
    W[findall(x->x>=NRight, N)] = ThetaRight .- 1 ./ (2*π) .* sin.(2 .* π .*
	ThetaRight)
    return W
end

function GetCM!(FFTLog::FFTLog)
    FFTLog.CM = FFTW.rfft(FFTLog.FXArray .* FFTLog.XArray .^ (-FFTLog.ν))
    FFTLog.M = Array(0:length(FFTLog.XArray)/2)
    FFTLog.CM = FFTLog.CM .* CWindow(FFTLog.M, floor(Int,
	FFTLog.CWindowWidth*FFTLog.N/2))
end

function EvaluateηM!(FFTLog::FFTLog)
    FFTLog.ηM = 2 .* π ./ (FFTLog.N .* FFTLog.DLnX) .* FFTLog.M
    return FFTLog.ηM
end

function GMVals(Mu::Float64, Q::Array{Complex{Float64},1})
    if(Mu+1+real(Q[1]) == 0)
        println("Γ(0) encountered. Please change another nu value! Try ν=1.1 .")
    end
    ImagQ = imag(Q)
    GM = zeros(ComplexF64, length(Q))
    cut = 200
    AsymQ = filter(x->abs(imag(x)) + abs(Mu) >cut, Q)
    AsymPlus = (Mu .+ 1 .+ AsymQ) ./2
    AsymMinus = (Mu .+ 1 .- AsymQ) ./2

    QGood = filter(x->abs.(imag(x)) .+ abs(Mu) .<=cut && x != Mu + 1 ,Q)
    AlphaPlus  = (Mu .+1 .+ QGood) ./2
    AlphaMinus = (Mu .+1 .- QGood) ./2
    GM[findall(x->abs.(imag(x)) .+ abs(Mu) .<=cut && x != Mu + 1 , Q)] .=
	SpecialFunctions.gamma.(AlphaPlus) ./ SpecialFunctions.gamma.(AlphaMinus)
    GM[findall(x->abs.(imag(x))+abs(Mu) .> cut && x != Mu + 1 , Q)] = exp.(
	(AsymPlus .- 0.5) .* log.(AsymPlus) .- (AsymMinus .- 0.5) .*
	log.(AsymMinus) .-
	AsymQ .+ 1/12 .* (1 ./AsymPlus .- 1 ./ AsymMinus) .+ 1/360 .* (1 ./
	AsymMinus .^3 .- 1 ./AsymPlus .^3) +1/1260 .* (1 ./AsymPlus .^ 5  .- 1 ./
	AsymMinus .^ 5) )
    return GM
end

function GL(Ell::Float64, ZArray::Vector{ComplexF64})
    GL = 2 .^ZArray .* GMVals(Ell+0.5, ZArray .- 1.5)
    return GL
end

function GL1(Ell::Float64, ZArray::Vector{ComplexF64})
    GL1 = -2 .^(ZArray .- 1) * (z_array -1) .* GMVals(Ell+0.5, ZArray .- 2.5)
    return GL1
end

function GL2(Ell::Float64, ZArray::Vector{ComplexF64})
    GL2 = 2 .^ (ZArray .- 2) .* (ZArray .- 2) .* (ZArray .- 1) .*
	GMVals(Ell+0.5, ZArray .- 3.5)
    return GL
end

function LogExtrap(X::Vector{Float64}, NExtrapLow::Int64, NExtrapHigh::Int64)
    DLnXLow = log(X[2]/X[1])
    DLnXHigh= log(reverse(X)[1]/reverse(X)[2])
    if NExtrapLow != 0
        LowX = X[1] .* exp.(DLnXLow .* Array(-NExtrapLow:-1))
        X = vcat(LowX, X)
    end
    if NExtrapHigh != 0
        HighX = last(X) .* exp.(DLnXLow .* Array(1:NExtrapHigh))
        X = vcat(X,HighX)
    end
    return X
end

function ZeroPad!(FFTLog::FFTLog)
    ZeroArray = zeros(FFTLog.NPad)
    FFTLog.XArray = LogExtrap(FFTLog.XArray, FFTLog.NPad, FFTLog.NPad)
    FFTLog.FXArray = vcat(ZeroArray, FFTLog.FXArray, ZeroArray)
end

function CheckNumberElements!(FFTLog::FFTLog)
    if iseven(FFTLog.N+1)
        deleteat!(FFTLog.XArray, length(FFTLog.XArray))
        deleteat!(FFTLog.FXArray, length(FFTLog.FXArray))
        FFTLog.N -= 1
        if (FFTLog.NExtrapHigh != 0)
            FFTLog.NExtrapHigh -= 1
        end
    end
end

function EvaluateFFTLog(FFTLog::FFTLog, Ell::Vector{Float64})
    FFTLog.XArray = LogExtrap(FFTLog.XArray, FFTLog.NExtrapLow,
	FFTLog.NExtrapHigh)
    FFTLog.FXArray = LogExtrap(FFTLog.FXArray, FFTLog.NExtrapLow,
	FFTLog.NExtrapHigh)
    ZeroPad!(FFTLog)
    GetCM!(FFTLog)
    EvaluateηM!(FFTLog)
    X0 = FFTLog.XArray[1]
    ZAr = FFTLog.ν .+ im .* FFTLog.ηM
    YArray = zeros(length(Ell), length(FFTLog.XArray))
    HMArray = zeros(ComplexF64, length(Ell), length(FFTLog.CM))
    FYArray = zeros(length(Ell), length(FFTLog.XArray))
    @inbounds for myl in 1:length(Ell)
        YArray[myl,:] = (Ell[myl] + 1) ./ reverse(FFTLog.XArray)
        HMArray[myl,:]  = FFTLog.CM .* (FFTLog.XArray[1] .* YArray[myl,1] ) .^
		(-im .*FFTLog.ηM) .* GL(Ell[myl], ZAr)
        FYArray[myl,:] = FFTW.irfft(conj(HMArray[myl,:]),
		length(FFTLog.XArray)) .* YArray[myl,:] .^ (-FFTLog.ν) .* sqrt(π) ./4
    end
    return YArray[:,FFTLog.NExtrapLow+FFTLog.NPad+1:FFTLog.NExtrapLow+
	FFTLog.NPad+FFTLog.OriginalLenght], FYArray[:,FFTLog.NExtrapLow+FFTLog.NPad+
	1:FFTLog.NExtrapLow+FFTLog.NPad+FFTLog.OriginalLenght]
end

function EvaluateFFTLogDJ(FFTLog::FFTLog, Ell::Vector{Float64})
    FFTLog.XArray = LogExtrap(FFTLog.XArray, FFTLog.NExtrapLow,
	FFTLog.NExtrapHigh)
    FFTLog.FXArray = LogExtrap(FFTLog.FXArray, FFTLog.NExtrapLow,
	FFTLog.NExtrapHigh)
    ZeroPad!(FFTLog)
    GetCM!(FFTLog)
    EvaluateηM!(FFTLog)
    X0 = FFTLog.XArray[1]
    ZAr = FFTLog.ν .+ im .* FFTLog.ηM
    YArray = zeros(length(Ell), length(FFTLog.XArray))
    HMArray = zeros(ComplexF64, length(Ell), length(FFTLog.CM))
    FYArray = zeros(length(Ell), length(FFTLog.XArray))
    @inbounds for myl in 1:length(Ell)
        YArray[myl,:] = (Ell[myl] + 1) ./ reverse(FFTLog.XArray)
        HMArray[myl,:]  = FFTLog.CM .* (FFTLog.XArray[1] .* YArray[myl,1] ) .^
		(-im .*FFTLog.ηM) .* GL1(Ell[myl], ZAr)
        FYArray[myl,:] = FFTW.irfft(conj(HMArray[myl,:]),
		length(FFTLog.XArray)) .* YArray[myl,:] .^ (-FFTLog.ν) .* sqrt(π) ./4
    end
    return YArray[:,FFTLog.NExtrapLow+FFTLog.NPad+1:FFTLog.NExtrapLow+
	FFTLog.NPad+FFTLog.OriginalLenght], FYArray[:,FFTLog.NExtrapLow+FFTLog.NPad+
	1:FFTLog.NExtrapLow+FFTLog.NPad+FFTLog.OriginalLenght]
end

function EvaluateFFTLogDDJ(FFTLog::FFTLog, Ell::Vector{Float64})
    FFTLog.XArray = LogExtrap(FFTLog.XArray, FFTLog.NExtrapLow,
	FFTLog.NExtrapHigh)
    FFTLog.FXArray = LogExtrap(FFTLog.FXArray, FFTLog.NExtrapLow,
	FFTLog.NExtrapHigh)
    ZeroPad!(FFTLog)
    GetCM!(FFTLog)
    EvaluateηM!(FFTLog)
    X0 = FFTLog.XArray[1]
    ZAr = FFTLog.ν .+ im .* FFTLog.ηM
    YArray = zeros(length(Ell), length(FFTLog.XArray))
    HMArray = zeros(ComplexF64, length(Ell), length(FFTLog.CM))
    FYArray = zeros(length(Ell), length(FFTLog.XArray))
    @inbounds for myl in 1:length(Ell)
        YArray[myl,:] = (Ell[myl] + 1) ./ reverse(FFTLog.XArray)
        HMArray[myl,:]  = FFTLog.CM .* (FFTLog.XArray[1] .* YArray[myl,1] ) .^
		(-im .*FFTLog.ηM) .* GL2(Ell[myl], ZAr)
        FYArray[myl,:] = FFTW.irfft(conj(HMArray[myl,:]),
		length(FFTLog.XArray)) .* YArray[myl,:] .^ (-FFTLog.ν) .* sqrt(π) ./4
    end
    return YArray[:,FFTLog.NExtrapLow+FFTLog.NPad+1:FFTLog.NExtrapLow+
	FFTLog.NPad+FFTLog.OriginalLenght], FYArray[:,FFTLog.NExtrapLow+FFTLog.NPad+
	1:FFTLog.NExtrapLow+FFTLog.NPad+FFTLog.OriginalLenght]
end
