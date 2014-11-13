function equilibrium(T::Float64, x::Array{Float64,1}, x_tot::Array{Float64,1},
                     stepTolerance = 1e-14, attempt = 1)
    #=
    Function for solving a multiphase equilibrium.
    Needs the following global variable:
        maxIter     :   maximum allowed iterations
        tolerance   :   tolerance for convergence check
    Input:
        T:      temperature [K]
                Float64
        x:      initial state guess [V_v; n_v] [[m^3]; [mol]]
                Array{Float64,1}
        x_tot:  total volume and mole numbers [V_tot; n_tot] [[m^3]; [mol]]
                Array{Float64,1}
    Output:
        x:              equilibrium state [V_v; n_v] [[m^3]; [mol]]
                        Array{Float64,1}
    =#
    iter      = 1                                  # Initialize iteration number
    converged = false                                         # Convergence flag
    ndxOld    = 1.0e9                                      # Initialize residual
    res       = 1.0e9                                          # Initialize norm
    f = 0
    p_v = 0
    mu_v = 0
    while !converged
        # Hessians
        H_v   = Hessian(T, x[1], x[2:end])
        H_l   = Hessian(T, x_tot[1] - x[1], x_tot[2:end] - x[2:end])
        
        # Chemical potentials
        mu_v  = chemicalPotential(T, x[1], x[2:end])
        mu_l  = chemicalPotential(T, x_tot[1] - x[1], x_tot[2:end] - x[2:end])

        # Pressure
        p_v   = pressure(T, x[1], x[2:end])
        p_l   = pressure(T, x_tot[1] - x[1], x_tot[2:end] - x[2:end])

        # Equations to be sovled
        f     = [p_l - p_v;
                 mu_l - mu_v]
        H     = H_v + H_l                                        # Total Hessian
        dx    = -H\f                                     # Solve for change in x
        x    += dx                                           # Apply change in x
        ndx   = norm(dx)                          # Check the size of the change

        # Check for invalid volumes and moles
        while reduceStep(x, x_tot) && maximum(dx) > eps()
            dx /= 2                                         # Reduce step length
            x  -= dx                                     # Apply new step length 
            if Debug
                println("Reducing step length")
            end
        end

        # Break if max iterations is reached
        if iter > maxIter
            if Debug
                warn("Maximum number of iterations reached: $maxIter
                     Current residual: $(norm(f))")
            end
            break
        end

        # Display info if debugging
        if Debug
            println("Residual:  $(norm(f))")
            println("Iteration: $iter")
            println("Vapor volume: $(x[1])")
            println("Min vapor volume: $(dot(x[2:end],b))")
            println("Liquid volume: $(x_tot[1]-x[1])")
            println("Min liquid volume: $(dot(x_tot[2:end]-x[2:end],b))")
            println("Volume fraction: $(x[1]/x_tot[1])")
            println("Mole fractions: $(x[2:end]./x_tot[2:end])")
        end

        
        iter += 1                                    # Increase iteration number
        converged = ndx < stepTolerance && ndx >= ndxOld            # Converged?
        ndxOld  = ndx                                      # Update old residual
        res = norm(f./([p_v;mu_v]+eps()))
        # Check for division by zero error
        if sum(isnan(x)) > 0
            error("Not a number")
        end
    end

    # End of iterations
    # Check for trivial solution
    if trivialSolution(x, x_tot)
        error("Trivial solution")
    end
    # Check if solution has converged
    res = norm(f./([p_v;mu_v]+eps()))
    if res > tolerance && attempt < 2
        # Try to solve for the other phase
        attempt += 1
        return x_tot - equilibrium(T, x_tot - x, x_tot, stepTolerance, attempt)
    elseif res > tolerance
        error("Not converged, norm: $res")
    end

    return x
end

function reduceStep(x::Array{Float64,1}, x_tot::Array{Float64,1})
    #=
    Check for negative mole numbers and volumes less than B
    =#

    n_v = x[2:end]                                                 # Vapor moles
    n_l = x_tot[2:end] - x[2:end]                                 # Liquid moles 
    V_v = x[1]                                                    # Vapor volume
    V_l = x_tot[1] - x[1]                                        # Liquid volume
    minV_v = dot(n_v, b)                                  # Minimal vapor volume
    minV_l = dot(n_l, b)                                 # Minimal liquid volume
    return (
        minimum(n_v) < 0 
     || minimum(n_l) < 0
     || V_v <= minV_v 
     || V_l <= minV_l
    )

end

function trivialSolution(x, x_tot)
    #=
    Check if a trivial solution has been found
    =#
    meanFrac    = mean(x./x_tot)                         # Mean of the fractions 
    numEquals   = 0                     # 
    for i = 1:length(x)
        numEquals += abs(x[i]/x_tot[i] - meanFrac) < tolerance
    end
    return numEquals == length(x)
end