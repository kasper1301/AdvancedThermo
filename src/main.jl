include("includes.jl")

#=
The solution strategy:
    0.  Start in the middle of the phase diagram (T = 180, V = 3e-4)
    1.  Loop through the volumes at a given temperature
        1.1  Increase volume until the fraction of vapor volume is almost 1
        1.2  Jump back to the initial volume
        1.3  Decrease volume until the fraction of vapor volume is almost 0
        1.4  Jump back to the inital volume
    2.  Increase the temperature and jump to 1.
    3.  Jump back to the initial temperature
    4.  Decrease the temperature and jump to 1.
=#

# Obtain initial solution
const n = [1.60; 94.50; 2.60; 0.81; 0.52]/100.03             # Total composition
initT   = 180.0                                            # Initial temperature
initV   = 3e-4                                                  # Initial volume
n_v     = 0.5n                                     # Guess for vapor composition
n_l     = n - n_v                                           # Liquid composition
V_l     = 1.1dot(n_l, b)                               # Guess for liquid volume
V_v     = initV - V_l                                             # Vapor volume
x_v     = [V_v, n_v]                                               # Vapor state
x       = [initV, n]                                               # Total state
initx_v = equilibrium(initT, x_v, x)         # Solve equilibrium for vapor state

# Start from initial solution
T       = initT
V       = initV
x_v     = initx_v
x_v0    = copy(x_v) # Keep a copy of the initial state

# Array for storing results
resultT = {Float64[]}                           # Array for storing temperatures
resultV = {Float64[]}                                # Array for storing volumes
resultX = {{Float64[]} for i in x}              # Arrays for storing vapor state

trace   = 1                                                       # Trace number
dV      = 1.0001        # The volume is updated by multiplication/division by dV
dT      = 0.1          # The temperature is update by addition/subtraction of dT

while true
    try
        include("tempLoop.jl")
    catch
        println("Done with increasing temperature")
        break
    end
end

# Reset state and to the calculations for decreasing temperature
T   = initT
V   = initV
x_v = x_v0
dT  = -dT

while true
    try
        include("tempLoop.jl")
    catch
        println("Done with decreasing temperature")
        break
    end
end

# Sort the volumes
for v in resultV
    sort!(v)
end

# Prepare array for storing the entropy
resultS = deepcopy(resultV)

# Create phase diagram and store entropy calculations
fig = figure()
for i in 1:length(resultV)
    if length(resultV[i]) > 0
        Vmax = maximum(resultV[i])
        idx  = indmax(resultV[i])
        # Check if Vmax is a bubble or a dew point
        if resultX[1][i][idx]/Vmax > 0.5
            # Bubble point
            semilogx(Vmax,resultT[i][idx],"bo",fillstyle="none")
        else
            # Dew point
            semilogx(Vmax,resultT[i][idx],"rx")
        end
        Vmin = minimum(resultV[i])
        idx  = indmin(resultV[i])
        # Check if Vmin is a bubble or a dew point
        if resultX[1][i][idx]/Vmin > 0.5
            semilogx(Vmin,resultT[i][idx],"bo",fillstyle="none")
        else
            semilogx(Vmin,resultT[i][idx],"rx")
        end
        # Calculate entropy
        for j in 1:length(resultV[i])
            V = resultV[i][j]
            T = resultT[i][j]
            V_v = resultX[1][i][j]
            n_v = Float64[resultX[c][i][j] for c in 2:length(resultX)]
            V_l = V - V_v
            n_l = n - n_v
            resultS[i][j] = Entropy(T,V_v,n_v) + Entropy(T,V_l,n_l)
        end
    end
end
xlabel(L"$V\quad\left[\SI{}{\cubic\meter}\right]$")
ylabel(L"$T\quad\left[\SI{}{\kelvin}\right]$")
legend(["Dew points", "Bubble points"])
xlim([1e-5,1e-2])
printfig(fig, filename="../fig/TV.tex")
close()

# Plot the volume-entropy diagram for different isotherms
fig = figure()
isotherms = Float64[170,180,190,200,210]

for i = 1:length(resultT)
    if in(resultT[i][1],isotherms)
        semilogx(resultV[i],resultS[i])
    end
end

xlabel(L"$V\quad\left[\SI{}{\cubic\meter}\right]$")
ylabel(L"$S\quad\left[\SI{}{\joule\per\kelvin\per\mole}\right]$")
printfig(fig, filename="../fig/VS.tex")