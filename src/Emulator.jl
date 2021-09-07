@kwdef mutable struct CℓEmulator
    NN
    InMinMax::Matrix{Float64} = zeros(8,2)
    OutMinMax::Array{Float64} = zeros(55,2)
    InputParams::Matrix{Float64} = zeros(8,100)
    DuplicationMatrix::Matrix{Float64} = zeros(10,10)
end

function SetℓGrid!(cosmoemu::CℓEmulator, cosmogrid::CosmologicalGrid)
    cosmoemu.InputParams = zeros(8, length(cosmogrid.ℓBinCenters))
end


function SetInputNN!(cosmoemu, cosmogrid::CosmologicalGrid, cosmology::w0waCDMCosmology)
    cosmoemu.InputParams[1,:] .= (cosmology.ΩM - cosmoemu.InMinMax[1,1])/
    (cosmoemu.InMinMax[1,2]-cosmoemu.InMinMax[1,1])
    cosmoemu.InputParams[2,:] .= (cosmology.ΩB - cosmoemu.InMinMax[2,1])/
    (cosmoemu.InMinMax[2,2]-cosmoemu.InMinMax[2,1])
    cosmoemu.InputParams[3,:] .= (cosmology.H0 - cosmoemu.InMinMax[3,1])/
    (cosmoemu.InMinMax[3,2]-cosmoemu.InMinMax[3,1])
    cosmoemu.InputParams[4,:] .= (cosmology.ns - cosmoemu.InMinMax[4,1])/
    (cosmoemu.InMinMax[4,2]-cosmoemu.InMinMax[4,1])
    cosmoemu.InputParams[5,:] .= (cosmology.σ8 - cosmoemu.InMinMax[5,1])/
    (cosmoemu.InMinMax[5,2]-cosmoemu.InMinMax[5,1])
    cosmoemu.InputParams[6,:] .= (cosmology.w0 - cosmoemu.InMinMax[6,1])/
    (cosmoemu.InMinMax[6,2]-cosmoemu.InMinMax[6,1])
    cosmoemu.InputParams[7,:] .= (cosmology.wa - cosmoemu.InMinMax[7,1])/
    (cosmoemu.InMinMax[7,2]-cosmoemu.InMinMax[7,1])
    cosmoemu.InputParams[8,:]  = (log10.(cosmogrid.ℓBinCenters) .-cosmoemu.InMinMax[8,1]) ./
    (cosmoemu.InMinMax[8,2]-cosmoemu.InMinMax[8,1])
end

function ComputeCℓ(cosmoemulator::CℓEmulator, cosmology::CosmoCentral.w0waCDMCosmology, cosmogrid::CosmoCentral.CosmologicalGrid)
    SetℓGrid!(cosmoemulator, cosmogrid)
    SetInputNN!(cosmoemulator, cosmogrid, cosmology)
    y = cosmoemulator.NN(cosmoemulator.InputParams)
    n_ℓ = length(cosmogrid.ℓBinCenters)
    mmm = zeros(n_ℓ, 10, 10)
    c = zeros(100)
    @avx for i in 1:55
        for l in 1:n_ℓ
            y[i,l] *= (cosmoemulator.OutMinMax[i,2]-cosmoemulator.OutMinMax[i,1])
            y[i,l] += (cosmoemulator.OutMinMax[i,1])
        end
    end
    for l in 1:n_ℓ
        my_exp!(y,l)
        mul!(c, cosmoemulator.DuplicationMatrix, @view(y[:,l]))
        mmm[l,:,:] = reshape(c, (10,10))./(
            cosmogrid.ℓBinCenters[l]*(cosmogrid.ℓBinCenters[l]+1.))

    end
    return mmm
end

function my_exp!(y,l)
    @. @views y[:,l] = 10 .^y[:,l]
end