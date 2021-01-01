var documenterSearchIndex = {"docs":
[{"location":"anotherPage/#The-CosmoCentral-Module","page":"Background","title":"The CosmoCentral Module","text":"","category":"section"},{"location":"anotherPage/#Module-Index","page":"Background","title":"Module Index","text":"","category":"section"},{"location":"anotherPage/","page":"Background","title":"Background","text":"Modules = [CosmoCentral]\nOrder   = [:constant, :type, :function, :macro]","category":"page"},{"location":"anotherPage/#Detailed-API","page":"Background","title":"Detailed API","text":"","category":"section"},{"location":"anotherPage/","page":"Background","title":"Background","text":"Modules = [CosmoCentral]\nOrder   = [:constant, :type, :function, :macro]","category":"page"},{"location":"anotherPage/#CosmoCentral.w0waCDMParameters","page":"Background","title":"CosmoCentral.w0waCDMParameters","text":"w0waCDMParameters\n\nThis struct contains the value of the cosmological parameters for w_0 w_aCDM cosmologies:\n\nw_0 and w_a, the parameters in the CPL parameterization\nOmega_M, Omega_B, Omega_DE, Omega_R, and Omega_k the density parameters for matter, baryons, Dark Energy, radiation, curvature\nn_s, the scalar spectral index\nM_nu, the sum of the neutrino mass eigenstates in eV\nsigma_8, the amplitude of the scalar fluctuations\nH_0, the value of the Hubble paramater\n\n\n\n\n\n","category":"type"},{"location":"anotherPage/#CosmoCentral.ComputeAdimensionalHubbleFactor-Tuple{Float64,CosmoCentral.w0waCDMCosmology}","page":"Background","title":"CosmoCentral.ComputeAdimensionalHubbleFactor","text":"ComputeAdimensionalHubbleFactor(z, params)\n\nThis function, given the value of the cosmological parameters, evaluate the Adimensional Hubble Factor for w_0 w_aCDM cosmologies.\n\n...\n\nArguments\n\nz::Float64 the redshift value at which evaluate the Adimensional Hubble Factor\nparams::w0waCDMCosmology  a collection of cosmological parameters, whose expression is given by\n\nE(z)=sqrtOmega_M(1+z)^3+Omega_R(1+z)^4+Omega_DE(1+z)^3(1+w_0+w_a)exp(-3w_a fracz1+z)+Omega_k(1+z)^2\n\nNotes\n\nThis expression is valid only for the CPL parameterization\n\n\n\n\n\n","category":"method"},{"location":"anotherPage/#CosmoCentral.ComputeComovingDistance-Tuple{Float64,CosmoCentral.w0waCDMCosmology}","page":"Background","title":"CosmoCentral.ComputeComovingDistance","text":"ComputeComovingDistance(z, params)\n\nThis function, given the value of the cosmological parameters, evaluate the Comoving Distance. It is evaluated as\n\nr(z)=fraccH_0int_0^z fracdxE(x)\n\n\n\n\n\n","category":"method"},{"location":"anotherPage/#CosmoCentral.ComputeDensityFunction-Tuple{Float64,CosmoCentral.AnalitycalDensity}","page":"Background","title":"CosmoCentral.ComputeDensityFunction","text":"ComputeDensityFunction(z, params)\n\nThis function returns the source density for a given redshift z\n\n\n\n\n\n","category":"method"},{"location":"anotherPage/#CosmoCentral.ComputeHubbleFactor-Tuple{Float64,CosmoCentral.w0waCDMCosmology}","page":"Background","title":"CosmoCentral.ComputeHubbleFactor","text":"ComputeHubbleFactor(z, params)\n\n...\n\nArguments\n\nz::Float64 the redshift value at which evaluate the Adimensional Hubble Factor\nparams::w0waCDMCosmology  a collection of cosmological parameters\n\nThis function, given the value of the cosmological parameters, evaluate the Hubble Factor for w_0 w_aCDM cosmologies, whose expression is given by\n\nH(z)=H_0sqrtOmega_M(1+z)^3+Omega_R(1+z)^4+Omega_DE(1+z)^3(1+w_0+w_a)exp(-3w_a fracz1+z)+Omega_k(1+z)^2\n\nNotes\n\nThis expression is valid only for the CPL parameterization\n\n\n\n\n\n","category":"method"},{"location":"anotherPage/#CosmoCentral.ComputeInstrumentResponse-Tuple{Float64,Float64,CosmoCentral.InstrumentResponseStruct}","page":"Background","title":"CosmoCentral.ComputeInstrumentResponse","text":"ComputeInstrumentResponse(z, params)\n\nThis function returns the instrument response.\n\n\n\n\n\n","category":"method"},{"location":"anotherPage/#CosmoCentral.NormalizeAnalitycalDensityStruct-Tuple{CosmoCentral.AnalitycalDensity}","page":"Background","title":"CosmoCentral.NormalizeAnalitycalDensityStruct","text":"ComputeDensityFunction(params)\n\nThis function modifies the normalization constant in the AnalitycalDensityStruct in order to have the same value of the surface density once integrated.\n\n\n\n\n\n","category":"method"},{"location":"anotherPage/#CosmoCentral.NormalizeConvolvedDensityStruct-Tuple{CosmoCentral.ConvolvedDensity}","page":"Background","title":"CosmoCentral.NormalizeConvolvedDensityStruct","text":"ComputeDensityFunction(params)\n\nThis function modifies the normalization constant in the AnalitycalDensityStruct in order to have the same value of the surface density once integrated.\n\n\n\n\n\n","category":"method"},{"location":"#CosmoCentral.jl","page":"Home","title":"CosmoCentral.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for CosmoCentral.jl. This is a Julia package to perform cosmological calculations.","category":"page"},{"location":"#Authors","page":"Home","title":"Authors","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Marco Bonici, Dipartimento di Fisica, Università degli Studi di Genova (UniGe)","category":"page"},{"location":"#Usage","page":"Home","title":"Usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"import CosmoCentral\n\nparams = CosmoCentral.w0waCDMParameters()\n\nCosmoCentral.ComputeAdimensionalHubbleFactor(z, params) # returns the adimensional Hubble factor\nCosmoCentral.ComputeHubbleFactor(z, params) # returns the Hubble factor\nCosmoCentral.ComputeComovingDistance(z, params) # returns the comoving distance","category":"page"},{"location":"#Contributing","page":"Home","title":"Contributing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Please make sure to update tests as appropriate.","category":"page"},{"location":"#License","page":"Home","title":"License","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"CosmoCentral is licensed under the MIT \"Expat\" license; see LICENSE for the full license text.","category":"page"}]
}
