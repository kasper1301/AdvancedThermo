#=
All the data is downloaded from DIPPR
    http://dippr.byu.edu/students/chemsearch.asp
=#
# Read the data
data        = readdlm("../data/componentData.csv",',','\n',skipstart=2)
# Component names
const global components  = convert(Array{String,1},data[:,1])[Ncomp]
# Critical temperatures [K]
const global tc          = convert(Array{Float64,1},data[:,2])[Ncomp]
# Critical pressures [Pa]
const global pc          = convert(Array{Float64,1},data[:,3])[Ncomp]
# Ideal gas reference enthalpies [J mol^{-1}]
const global h_ref       = convert(Array{Float64,1},data[:,4])[Ncomp]/1000
# Ideal gas reference entropies [J mol^{-1} K^{-1}]
const global s_ref       = convert(Array{Float64,1},data[:,5])[Ncomp]/1000
# Ideal gas heat capacities (assumed constant for now)[J mol^{-1} K^{-1}]
const global cpCoeffs    = convert(Array{Float64,2},data[:,6:end-1])[Ncomp,:]./
                            [1e3 1e3 1 1e3 1 1e3 1]
const global DIPPRnumber = convert(Array{Int64,1},data[:,end])[Ncomp]
# Redlich-Kwong constants b [m^3 mol^{-1}]
const global b           = (2^(1/3)-1)/3*R*tc./pc
# Redlich-Kwong constants sqrt(a) [J^{1/2} K^{1/4} m^{-3/2} mol^{-1}]
const global a           = sqrt(1/(9*(2^(1/3)-1))*R^2*tc.^(5/2)./pc)