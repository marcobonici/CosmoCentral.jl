@kwdef mutable struct CℓEmulator
    NN
    InMinMax::Matrix{Float64} = zeros(8,2)
    OutMinMax::Array{Float64} = zeros(55,2)
    InputParams::Matrix{Float64} = zeros(8,100)
    DuplicationMatrix::Matrix{Float64} = zeros(10,10)
end

function SetInputNN!(array_input, inMinMax, cosmogrid::CosmoCentral.CosmologicalGrid, cosmology::CosmoCentral.w0waCDMCosmology)
    array_input[1,:] .= (cosmology.ΩM - inMinMax[1,1])/(inMinMax[1,2]-inMinMax[1,1])
    array_input[2,:] .= (cosmology.ΩB - inMinMax[2,1])/(inMinMax[2,2]-inMinMax[2,1])
    array_input[3,:] .= (cosmology.H0 - inMinMax[3,1])/(inMinMax[3,2]-inMinMax[3,1])
    array_input[4,:] .= (cosmology.ns - inMinMax[4,1])/(inMinMax[4,2]-inMinMax[4,1])
    array_input[5,:] .= (cosmology.σ8 - inMinMax[5,1])/(inMinMax[5,2]-inMinMax[5,1])
    array_input[6,:] .= (cosmology.w0 - inMinMax[6,1])/(inMinMax[6,2]-inMinMax[6,1])
    array_input[7,:] .= (cosmology.wa - inMinMax[7,1])/(inMinMax[7,2]-inMinMax[7,1])
    array_input[8,:]  = (log10.(cosmogrid.ℓBinCenters) .-inMinMax[8,1]) ./ (inMinMax[8,2]-inMinMax[8,1])
end

function ComputeCℓ(cosmoemulator::CℓEmulator, cosmology::CosmoCentral.w0waCDMCosmology, cosmogrid::CosmoCentral.CosmologicalGrid)
    SetInputNN!(cosmoemulator.InputParams, cosmoemulator.InMinMax, cosmogrid, cosmology)
    y = cosmoemulator.NN(cosmoemulator.InputParams)
    mmm = zeros(100, 10, 10)
    c = zeros(100)
    @avx for i in 1:55
        for l in 1:100
            y[i,l] *= (cosmoemulator.OutMinMax[i,2]-cosmoemulator.OutMinMax[i,1])
            y[i,l] += (cosmoemulator.OutMinMax[i,1])
        end
    end
    mmm = zeros(100, 10, 10)
    c = zeros(100)
    for i in 1:100
        mul!(c, cosmoemulator.DuplicationMatrix, 10 .^y[:,i])
        mmm[i,:,:] = reshape(c, (10,10))
    end
    return mmm
end