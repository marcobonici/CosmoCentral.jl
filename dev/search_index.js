var documenterSearchIndex = {"docs":
[{"location":"BackgroundUniverse/#Background-Universe","page":"Background Universe","title":"Background Universe","text":"","category":"section"},{"location":"BackgroundUniverse/","page":"Background Universe","title":"Background Universe","text":"Some of the most basics quantities in cosmology are the Hubble factor H(z) and the comoving distance r(z). In this section are presented the functions which evaluates them.","category":"page"},{"location":"BackgroundUniverse/#Hubble-factor","page":"Background Universe","title":"Hubble factor","text":"","category":"section"},{"location":"BackgroundUniverse/","page":"Background Universe","title":"Background Universe","text":"CosmoCentral.ComputeAdimensionalHubbleFactor\nCosmoCentral.ComputeHubbleFactor","category":"page"},{"location":"BackgroundUniverse/#CosmoCentral.ComputeAdimensionalHubbleFactor","page":"Background Universe","title":"CosmoCentral.ComputeAdimensionalHubbleFactor","text":"ComputeAdimensionalHubbleFactor(z::Float64, params::w0waCDMCosmology)\n\nThis function, given the value of the cosmological parameters, evaluate the Adimensional Hubble Factor for w_0 w_aCDM cosmologies. The analitycal expression is given by:\n\nE(z)=sqrtOmega_M(1+z)^3+Omega_R(1+z)^4+Omega_DE(1+z)^3(1+w_0+w_a)expleft(-3w_a fracz1+zrigth)+Omega_k(1+z)^2\n\nExample\n\nusing CosmoCentral\nparams = CosmoCentral.w0waCDMStruct()\nz = 1.\nCosmoCentral.ComputeAdimensionalHubbleFactor(z, params)\n\nwarning: Warning\nThis expression is valid only for the     CPL parameterization     of the Dark Energy Equation of State.\n\n\n\n\n\n","category":"function"},{"location":"BackgroundUniverse/#CosmoCentral.ComputeHubbleFactor","page":"Background Universe","title":"CosmoCentral.ComputeHubbleFactor","text":"ComputeHubbleFactor(z::Float64, params::w0waCDMCosmology)\n\nThis function, given the value of the cosmological parameters, evaluate the Hubble Factor for w_0 w_aCDM cosmologies, whose expression is given by\n\nH(z)=H_0sqrtOmega_M(1+z)^3+Omega_R(1+z)^4+\nOmega_DE(1+z)^3(1+w_0+w_a)exp(-3w_a fracz1+z)+Omega_k(1+z)^2\n\nExample\n\nusing CosmoCentral\nparams = CosmoCentral.w0waCDMStruct(H0=67.)\nz = 1.\nCosmoCentral.ComputeHubbleFactor(z, params)\n\n\n\n\n\n","category":"function"},{"location":"BackgroundUniverse/#Distances","page":"Background Universe","title":"Distances","text":"","category":"section"},{"location":"BackgroundUniverse/","page":"Background Universe","title":"Background Universe","text":"CosmoCentral.ComputeComovingDistance","category":"page"},{"location":"BackgroundUniverse/#CosmoCentral.ComputeComovingDistance","page":"Background Universe","title":"CosmoCentral.ComputeComovingDistance","text":"ComputeComovingDistance(z::Float64, params::w0waCDMCosmology)\n\nThis function, given the value of the cosmological parameters, evaluate the Comoving Distance. It is evaluated as:\n\nr(z)=fraccH_0int_0^z fracdxE(x)\n\nExample\n\nusing CosmoCentral\nparams = CosmoCentral.w0waCDMStruct(H0=67.)\nz = 1.\nCosmoCentral.ComputeComovingDistance(z, params)\n\n\n\n\n\n","category":"function"},{"location":"MathUtils/#Math-Utils","page":"Math Utils","title":"Math Utils","text":"","category":"section"},{"location":"MathUtils/","page":"Math Utils","title":"Math Utils","text":"In order to perform some calculations, we use some custom function, which are listed here.","category":"page"},{"location":"MathUtils/","page":"Math Utils","title":"Math Utils","text":"CosmoCentral.LogSpaced","category":"page"},{"location":"MathUtils/#CosmoCentral.LogSpaced","page":"Math Utils","title":"CosmoCentral.LogSpaced","text":"LogSpaced(min::Float64, max::Float64, n::Int64)\n\nThis function evaluates n points, logarithmically spaced between     min and max.\n\n\n\n\n\n","category":"function"},{"location":"CosmologicalStructure/#Cosmological-Structure","page":"Cosmological Structure","title":"Cosmological Structure","text":"","category":"section"},{"location":"CosmologicalStructure/","page":"Cosmological Structure","title":"Cosmological Structure","text":"In order to organize the code, some custom types are used. This section is under heavy developement and a lot of things may change.","category":"page"},{"location":"CosmologicalStructure/","page":"Cosmological Structure","title":"Cosmological Structure","text":"CosmoCentral.w0waCDMStruct\nCosmoCentral.CosmoGridStruct","category":"page"},{"location":"CosmologicalStructure/#CosmoCentral.w0waCDMStruct","page":"Cosmological Structure","title":"CosmoCentral.w0waCDMStruct","text":"w0waCDMStruct(w0::Float64 = -1, wa::Float64 = 0, ΩM::Float64 = 0.32,\nΩB::Float64  = 0.05, ΩDE::Float64 = 0.68, Ωk::Float64  = 0.,\nΩr::Float64  = 0., ns::Float64  = 0.96, Mν::Float64  = 0.06,\nσ8::Float64  = 0.816, H0::Float64  = 67.)\n\nThis struct contains the value of the cosmological parameters for w_0 w_aCDM cosmologies:\n\nw_0 and w_a, the parameters in the CPL parameterization\nOmega_M, Omega_B, Omega_DE, Omega_R, and Omega_k the density parameters for matter, baryons, Dark Energy, radiation, curvature\nn_s, the scalar spectral index\nM_nu, the sum of the neutrino mass eigenstates in eV\nsigma_8, the amplitude of the scalar fluctuations\nH_0, the value of the Hubble paramater\n\n\n\n\n\n","category":"type"},{"location":"CosmologicalStructure/#CosmoCentral.CosmoGridStruct","page":"Cosmological Structure","title":"CosmoCentral.CosmoGridStruct","text":"CosmoGridStruct(zgrid::Vector{Float64} = Array(LinRange(0.001, 2.5, 300)),\nkgrid::Vector{Float64} = LogSpaced(1e-5, 50., 1000))\n\nThis struct contains the value of the Cosmological Grid, both in k and z.\n\n\n\n\n\n","category":"type"},{"location":"SourceDensity/#Source-Density","page":"Source Density","title":"Source Density","text":"","category":"section"},{"location":"SourceDensity/","page":"Source Density","title":"Source Density","text":"In this page are presented the structures and functions used to deal with source densities and are listed here.","category":"page"},{"location":"SourceDensity/","page":"Source Density","title":"Source Density","text":"Pages = [\"SourceDensity.md\"]","category":"page"},{"location":"SourceDensity/#Source-Density-2","page":"Source Density","title":"Source Density","text":"","category":"section"},{"location":"SourceDensity/","page":"Source Density","title":"Source Density","text":"In the evaluation of Angular Coefficients, central quantities are the source densities. In this section are presented the custom types and function used to deal with the source densities.","category":"page"},{"location":"SourceDensity/","page":"Source Density","title":"Source Density","text":"CosmoCentral.AnalitycalDensityStruct\nCosmoCentral.ComputeDensityFunction\nCosmoCentral.NormalizeAnalitycalDensityStruct","category":"page"},{"location":"SourceDensity/#CosmoCentral.AnalitycalDensityStruct","page":"Source Density","title":"CosmoCentral.AnalitycalDensityStruct","text":"AnalitycalDensityStruct(z0::Float64 = 0.9/sqrt(2.), zmin::Float64 = 0.001,\nzmax::Float64 = 2.5, surfacedensity::Float64 = 30.,\nnormalization::Float64 = 1.)\n\nThis struct contains the parameters of the source galaxy density as given by the official Euclid forecast, whose expression is given by:\n\nn(z)proptoleft(fraczz_0right)^2\nexpleft(-left(fraczz_0right)^-32right)\n\nThe parameters contained in this struct are\n\nz_min and z_max, the minimum and the maximum redshift considered\nz_0, the parameter present in the galaxy distribution\nsurfacedensity , the value of the galaxy source density integrated between z_min and z_max\nnormalization, the value of parameter which multiplies the source dennsity in order to match the correct surface density\n\n\n\n\n\n","category":"type"},{"location":"SourceDensity/#CosmoCentral.ComputeDensityFunction","page":"Source Density","title":"CosmoCentral.ComputeDensityFunction","text":"ComputeDensityFunction(z::Float64, AnalitycalDensity::AnalitycalDensity = AnalitycalDensityStruct())\n\nThis function returns the source density for a given redshift z.\n\n\n\n\n\n","category":"function"},{"location":"SourceDensity/#CosmoCentral.NormalizeAnalitycalDensityStruct","page":"Source Density","title":"CosmoCentral.NormalizeAnalitycalDensityStruct","text":"NormalizeAnalitycalDensityStruct(densityparameters::AnalitycalDensity)\n\nThis function normalize AnalitycalDensityStruct in order to have the correct value of the surface density once integrated.\n\n\n\n\n\n","category":"function"},{"location":"SourceDensity/#Convolved-Source-Density","page":"Source Density","title":"Convolved Source Density","text":"","category":"section"},{"location":"SourceDensity/","page":"Source Density","title":"Source Density","text":"In real surveys we do not deal with the exactr distributions due to errors in the measurement of the source redshifts. The redshift errors are accounted for convolving the source density with a redshift measurement error.","category":"page"},{"location":"SourceDensity/#Intrument-Response","page":"Source Density","title":"Intrument Response","text":"","category":"section"},{"location":"SourceDensity/","page":"Source Density","title":"Source Density","text":"CosmoCentral.InstrumentResponseStruct\nCosmoCentral.ComputeInstrumentResponse","category":"page"},{"location":"SourceDensity/#CosmoCentral.InstrumentResponseStruct","page":"Source Density","title":"CosmoCentral.InstrumentResponseStruct","text":"InstrumentResponseStruct(cb::Float64 = 1.0, zb::Float64 = 0.0,\nσb::Float64 = 0.05, co::Float64 = 1.0, zo::Float64 = 0.1,\nσo::Float64 = 0.05, fout::Float64 = 0.1)\n\nWhen we measure the redshift of a galaxy with redshit z, we will measure a redshift z_p with a probability given by the following expression:\n\np(z_pz)  = frac1-f_outsqrt2 pi sigma_b(1+z) exp left(\n-frac12left(fracz-c_b z_b-z_bsigma_b(1+z)right)^2\nright) + fracf_outsqrt2 pi sigma_mathrmo(1+z) exp\nleft(-frac12left(fracz-c_o z_p-z_osigma_o(1+z)\nright)^2right)\n\nThis struct contains all these parameters.\n\n\n\n\n\n","category":"type"},{"location":"SourceDensity/#CosmoCentral.ComputeInstrumentResponse","page":"Source Density","title":"CosmoCentral.ComputeInstrumentResponse","text":"ComputeInstrumentResponse(z::Float64, zp::Float64, instrumentresponse::InstrumentResponse)\n\nThis function computes the probability that we actually measure a redshift z_p if the real redshift is z.\n\n\n\n\n\n","category":"function"},{"location":"SourceDensity/#Convolved-Source-Density-2","page":"Source Density","title":"Convolved Source Density","text":"","category":"section"},{"location":"SourceDensity/","page":"Source Density","title":"Source Density","text":"CosmoCentral.ConvolvedDensityStruct\nCosmoCentral.ComputeConvolvedDensityFunction\nCosmoCentral.NormalizeConvolvedDensityStruct\nCosmoCentral.ComputeDensityFunctionConvolvedGrid","category":"page"},{"location":"SourceDensity/#CosmoCentral.ConvolvedDensityStruct","page":"Source Density","title":"CosmoCentral.ConvolvedDensityStruct","text":"ConvolvedDensityStruct(AnalitycalDensity::AnalitycalDensity = AnalitycalDensityStruct(),\nInstrumentResponse::InstrumentResponse = InstrumentResponseStruct()\nzbinarray::Vector{Float64} = Array([0.001, 0.418, 0.560, 0.678, 0.789, 0.900, 1.019, 1.155, 1.324, 1.576, 2.50])\ndensityarraynormalization::Vector{Float64} = ones(length(zbinarray)-1)\ndensitygridarray::AbstractArray{Float64, 2} = ones(length(zbinarray)-1, 300))\n\nIn order to take into account the error in the redshift measurement, the source density is convolved with the InstrumentResponseStruct, according to the following equation\n\nn_i(z)=fracint_z_i^-^z_i^+\nmathrmd z_mathrmp n(z) p left(z_mathrmp\nmid zright)int_z_min ^z_max  mathrmd z\nint_z_i^-^z_i^+ mathrmd z_mathrmp n(z) p\nleft(z_mathrmp mid zright)\n\n\n\n\n\n","category":"type"},{"location":"SourceDensity/#CosmoCentral.ComputeConvolvedDensityFunction","page":"Source Density","title":"CosmoCentral.ComputeConvolvedDensityFunction","text":"ComputeConvolvedDensityFunction(z::Float64, i::Int64, convolveddensity::ConvolvedDensity)\n\nThis function computes the Convolved density function for a single bin at a given redshift z.\n\n\n\n\n\n","category":"function"},{"location":"SourceDensity/#CosmoCentral.NormalizeConvolvedDensityStruct","page":"Source Density","title":"CosmoCentral.NormalizeConvolvedDensityStruct","text":"NormalizeConvolvedDensityStruct(convolveddensity::ConvolvedDensity)\n\nThis function normalizes ConvolvedDensity such that the integrals of the convolved densities are normalized to 1.\n\n\n\n\n\n","category":"function"},{"location":"SourceDensity/#CosmoCentral.ComputeDensityFunctionConvolvedGrid","page":"Source Density","title":"CosmoCentral.ComputeDensityFunctionConvolvedGrid","text":"ComputeDensityFunction(CosmoGrid::CosmoGrid, ConvolvedDensityStruct::ConvolvedDensity)\n\nThis function computes the convolved density function for all tomographic bins on the z-grid provided by CosmoGrid.\n\n\n\n\n\n","category":"function"},{"location":"#CosmoCentral.jl","page":"Home","title":"CosmoCentral.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CosmoCentral is a Julia package to perform cosmological calculations. Actually it can evaluate:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Background quantities for w_0 w_aCDM cosmologies\nSource densities, with an analitycal in input","category":"page"},{"location":"","page":"Home","title":"Home","text":"We aim to include also:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Angular Correlation functions, C_ℓ's, for several probes (e.g., Weak Lensing, Galaxy Clustering)\nFisher Matrix evaluation to perform forecasts","category":"page"},{"location":"#Authors","page":"Home","title":"Authors","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Marco Bonici, Dipartimento di Fisica, Università degli Studi di Genova (UniGe)","category":"page"},{"location":"#Usage","page":"Home","title":"Usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Here is an example of how CosmoCentral can be used to evaluate background quantities.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using CosmoCentral\nparams = CosmoCentral.w0waCDMStruct()\nz = 1.\nCosmoCentral.ComputeAdimensionalHubbleFactor(z, params)\nCosmoCentral.ComputeHubbleFactor(z, params)\nCosmoCentral.ComputeComovingDistance(z, params)","category":"page"},{"location":"#Contributing","page":"Home","title":"Contributing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Please make sure to update tests as appropriate.","category":"page"},{"location":"#License","page":"Home","title":"License","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CosmoCentral is licensed under the MIT \"Expat\" license; see LICENSE for the full license text.","category":"page"}]
}