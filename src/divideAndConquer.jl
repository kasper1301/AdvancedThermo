include("includes.jl")

const global n_tot = [1.60; 94.50; 2.60; 0.81; 0.52]/100.03

V_tot   = [7.83238e-5, 1e-2]
T       = [60.0, 170.0]

resultT = Float64[]
resultV = Float64[]
resultX = {Float64[] for i in n_tot}

function divideAndConquer(T, V, resultT, resultV, resultX, depth = 1)
    #=
    A divide and conquer algorithm for solving phase equilibrium problems.
    The domain is a triangle in temperature and volume. If the calculation in 
    all or none of the corners are successful, the calculation is stopped.
    If there are 1 or 2 successful calculations, the triangle is divided into
    three smaller triangles. 

    All successful calculations are stored.

    Input:
        T: two T-coordinates
        V: two V-coordinates
        resultT: A vector with temperature for successful calculations
        resultV: A vector with total volume for successful calculations
        resultX: An array of vectors with resulting vapor solution
        depth: The recursion depth
    =#
    # Guess the vapor composition
    #n_v = [0.7, 0.6, 0.5, 0.3, 0.2].*n_tot
    n_v =  0.5.*n_tot
    # Calculate the liquid composition
    n_l = n_tot - n_v
    # Guess a small liquid volume
    V_l = 1.1dot(n_l,b)

    # Create a vector with temperature for the three corners
    T = [T[1], T[1], T[2]]
    # Create a vector with the total volume for the three corners
    V = [V[1], V[2], 0.5*(V[1] + V[2])]

    # Check for minimum differences
    if abs(T[3] - T[1]) < 0.1 || abs(V[2] - V[1]) < 1e-6
        return
    end

    nSuccess = 0                  # Initialize number of successful calculations

    # Try to solve equilibrium condition in all corners
    for i = 1:3
        try
            V_v = V[i] - V_l                              # Guess a vapor volume
            x_tot = [V[i]; n_tot]                                  # Total state
            
            x_v = equilibrium(T[i], [V_v, n_v], x_tot)   # Calculate equilibrium

            # Store the results
            push!(resultT, T[i]) 
            push!(resultV, V[i])
            for r in 1:length(resultX)
                push!(resultX[r],x_v[r])
            end
            nSuccess += 1            # Increase number of successful calculatons
        catch
            # Vapor phase did not converge, try to solve the liquid phase
            try
                x_tot = [V[i]; n_tot]
                x_l = equilibrium(T[i], [V_l, n_l], x_tot)
                push!(resultT, T[i]) 
                push!(resultV, V[i])
                for r in 1:length(resultX)
                    push!(resultX[r],x_tot[r] - x_l[r])
                end
                nSuccess += 1        # Increase number of successful calculatons
            end
        end
    end

    # Check if triangle is to be split
    if nSuccess < 3 && nSuccess > 0
        depth += 1                                    # Increase recursion depth
        println("Depth: $depth")
        divideAndConquer([T[1],0.5(T[1]+T[3])],
                         [V[1],V[3]],
                         resultT,
                         resultV,
                         resultX,
                         depth)
        divideAndConquer([T[1],0.5(T[1]+T[3])],
                         [V[3],V[2]],
                         resultT,
                         resultV,
                         resultX,
                         depth)
        divideAndConquer([0.5(T[1]+T[3]),T[3]],
                         [0.5(V[1]+V[3]),0.5(V[2]+V[3])], 
                         resultT,
                         resultV,
                         resultX,
                         depth)
    end
    if nSuccess == 3
        meanVaporFraction = mean(resultX[1][end-2:end]./resultV[end-2:end])
        color = "$meanVaporFraction"
        fill(V,T,color)
    end
    return
end

divideAndConquer(T, V_tot, resultT, resultV, resultX)