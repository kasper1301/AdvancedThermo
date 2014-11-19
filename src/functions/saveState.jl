function saveState!(T, V, x, resultT, resultV, resultX, index=1)
    #=
    Save the current T, V and x into the result arrays
    =#
    if length(resultT) == index
        push!(resultT[index], T)
        push!(resultV[index], V)
        for i = 1:length(x)
            push!(resultX[i][index], x[i])
        end
    else
        push!(resultT, Float64[])
        push!(resultV, Float64[])
        for i = 1:length(x)
            push!(resultX[i], Float64[])
        end
        saveState!(T, V, x, resultT, resultV, resultX, index)
    end
end
