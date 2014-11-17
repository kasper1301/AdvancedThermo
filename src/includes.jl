# Constants
const global R          = 8.3145              # Gas constant [J K^{-1} mol^{-1}]
const global T_ref      = 298.15                     # Reference temperature [K]
const global p_ref      = 1.0e5                        # Reference pressure [Pa]
const global Ncomp      = [2, 3, 4, 5, 6]                    # Component indices
const global Debug      = false                                     # Debug flag
const global maxIter    = 1000               # Max iterations for Newton-Raphson
const global tolerance  = 1e-6                    # Tolerance for Newton-Raphson

# Includes
require("componentData.jl")                            # Data for each component
require("functions/idealGas.jl")                # Set of functions for ideal gas
require("functions/redlichKwong.jl")    # Set of functions for the Redlich-Kwong 
                                        # equation of state
require("functions/equilibrium.jl")  # Newton-Raphson iteration for equilibrium 
                                     # calculation
require("functions/saveState.jl")        # Function for updating solution arrays

const global e_n = ones(Int64,size(components))                    # Unit vector
