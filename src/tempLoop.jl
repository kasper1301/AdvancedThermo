    println("T: $T")
    # Solve for increasing volume
    while true
        try
            V    = V*dV                                                 # Update V
            x[1] = V
            println("V: $V")
            x_v  = equilibrium(T, x_v, x)                    # Solve equilibrium
            saveState!(T, V, x_v, resultT, resultV, resultX, trace) # Save result
        catch
            warn("Vapor volume fraction: $(x_v[1]/V)")
            break
        end
    end
    # Reset state
    V       = initV
    x[1]    = V
    x_v     = initx_v
    # Solve for decreasing volume
    while true
        try
            V    = V/dV                                               # Update V
            x[1] = V
            println("V: $V")
            x_v  = equilibrium(T, x_v, x)                    # Solve equilibrium
            saveState!(T, V, x_v, resultT, resultV, resultX, trace) # Save result
        catch
            warn("Vapor volume fraction: $(x_v[1]/V)")
            break
        end
    end
    # Calculate new initial state
    V       = initV                                               # Reset volume
    T      += dT                                            # Update temperature
    x[1]    = V                                             # Update total state
    n_v     = initx_v[2:end]                           # Vapor composition guess
    n_l     = n - n_v                                       # Liquid composition
    V_l     = 1.1dot(n_l,b)                                # Liquid volume guess
    V_v     = V - V_l                                             # Vapor volume
    initx_v = equilibrium(T, [V_v, n_v], x)         # Update initial vapor state
    x_v     = initx_v                                      # Update vapor states
    trace  += 1                                                 # Increase trace
