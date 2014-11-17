#=
    Functions for calculating properties of an ideal gas. The variables are
    temperature, volume and number of moles. 

    Needs the following global variables:
        R       :   The ideal gas constant [J K^{-1} mol^{-1}]  Float64
        T_ref   :   The reference temperature                   Float64
        p_ref   :   The reference pressure                      Float64
        h_ref   :   The reference enthalpies                    Array{Float64,1}
        s_ref   :   The reference entropies                     Array{Float64,1}
        cp      :   The heat capacities                         Array{Float64,1}
        e_n     :   A unit vector with the same size as n       Array{Int64,1}

    Functions:
        calculateIdealEnthalpy  :   The enthalpy of an ideal gas
        calculateIdealEntropy   :   The entropy of an ideal gas
        calculateIdealPotential :   The chemical potential of an ideal gas
        idealHelmholtz          :   The Helmholtz' energy of an ideal gas
        Aig_T                   :   dA^ig/dT = -S^ig
        Aig_V                   :   dA^ig/dV = -p^ig
        Aig_n                   :   dA^ig/dn = mu^ig
        Aig_T_T                 :   d2A^ig/dTdT
        Aig_T_V                 :   d2A^ig/dTdV
        Aig_T_n                 :   d2A^ig/dTdn
        Aig_V_V                 :   d2A^ig/dVdV
        Aig_n_V                 :   d2A^ig/dndV
        Aig_n_n                 :   d2A^ig/dndn
        idealHessian            :   The Hessian matrix for Helmholtz' energy of 
                                    an ideal gas. Based on volume and mole 
                                    vector only, temperature is given.
=#

function integrateCp(T::Float64)
    #=
    Integrate the cp from T_ref to T
    =#
    intCp = zeros(h_ref)
    for i = 1:length(Ncomp)
        if DIPPRnumber[i] == 107
            # Use DIPPR equation number 107
            intCp[i] = integrateCp107(T, i)
        elseif DIPPRnumber[i] == 127
            # Use DIPPR equation number 127
            intCp[i] = integrateCp127(T, i)
        else
            # Assume constant cp
            intCp[i] = cpCoeffs[i,1]*(T - T_ref)
        end
    end
    return intCp
end

function integrateCp107(T::Float64, i::Int64)
    #=
    Integrate cp from T_ref to T for component number i using DIPPR equation
    number 107
    =#
    A = cpCoeffs[i,1]
    B = cpCoeffs[i,2]
    C = cpCoeffs[i,3]
    D = cpCoeffs[i,4]
    E = cpCoeffs[i,5]
    A*(T-T_ref) + 
    2B*C*(1/(exp(2C/T) - 1) - 1/(exp(2C/T_ref) - 1)) + 
    2D*E*(1/(exp(2E/T) + 1) - 1/(exp(2E/T_ref) + 1))
end

function integrateCp127(T::Float64, i::Int64)
    #=
    Integrate cp from T_ref to T for component number i using DIPPR equation
    number 127
    =#
    A = cpCoeffs[i,1]
    B = cpCoeffs[i,2]
    C = cpCoeffs[i,3]
    D = cpCoeffs[i,4]
    E = cpCoeffs[i,5]
    F = cpCoeffs[i,6]
    G = cpCoeffs[i,7]
    A*(T - T_ref) +
    B*C*(1/(exp(C/T) - 1) - 1/(exp(C/T_ref) - 1)) + 
    D*E*(1/(exp(E/T) - 1) - 1/(exp(E/T_ref) - 1)) + 
    F*G*(1/(exp(G/T) - 1) - 1/(exp(G/T_ref) - 1))
end

function integrateCpByT(T::Float64)
    #=
    Integrate cp/T from T_ref to T
    =#
    intCp = zeros(s_ref)
    for i = 1:length(Ncomp)
        if DIPPRnumber[i] == 107
            # Use DIPPR equation number 107
            intCp[i] = integrateCpByT107(T, i)
        elseif DIPPRnumber[i] == 127
            # Use DIPPR equation number 127
            intCp[i] = integrateCpByT127(T, i)
        else
            # Assume constant cp
            intCp[i] = cpCoeffs[i,1]*log(T/T_ref)
        end
    end
    return intCp
end

function integrateCpByT107(T::Float64, i::Int64)
    #=
    Integrate cp/Tfrom T_ref to T for component number i using DIPPR equation
    number 107
    =#
    A = cpCoeffs[i,1]
    B = cpCoeffs[i,2]
    C = cpCoeffs[i,3]
    D = cpCoeffs[i,4]
    E = cpCoeffs[i,5]
    A*log(T/T_ref) - 
    B*log((exp((2C)/T) - 1)/(exp((2C)/T_ref) - 1)) +
    2B*C*
    (
        1/(T*(exp(2C/T) - 1)) - 1/(T_ref*(exp(2C/T_ref) - 1)) + 
        (1/T - 1/T_ref)
    ) +
    D*log((exp((2E)/T) - 1)/(exp((2E)/T_ref) - 1)) +
    2D*E*
    (
        1/(T*(exp(2E/T) - 1)) - 1/(T_ref*(exp(2E/T_ref) - 1)) - 
        (1/T - 1/T_ref)
    )
end

function integrateCpByT127(T::Float64, i::Int64)
    #=
    Integrate cp/T from T_ref to T for component number i using DIPPR equation
    number 127
    =#
    A = cpCoeffs[i,1]
    B = cpCoeffs[i,2]
    C = cpCoeffs[i,3]
    D = cpCoeffs[i,4]
    E = cpCoeffs[i,5]
    F = cpCoeffs[i,6]
    G = cpCoeffs[i,7]
    A*log(T/T_ref) - 
    B*log((exp(C/T) - 1)/(exp(C/T_ref) - 1)) -
    D*log((exp(E/T) - 1)/(exp(E/T_ref) - 1)) -
    F*log((exp(G/T) - 1)/(exp(G/T_ref) - 1)) -
    B*C*
    (
        exp(C/T)/(T*(1 - exp(C/T))) - 
        exp(C/T_ref)/(T_ref*(1 - exp(C/T_ref)))
    ) - 
    D*E*
    (
        exp(E/T)/(T*(1 - exp(E/T))) - 
        exp(E/T_ref)/(T_ref*(1 - exp(E/T_ref)))
    ) - 
    F*G*
    (
        exp(G/T)/(T*(1 - exp(G/T))) - 
        exp(G/T_ref)/(T_ref*(1 - exp(G/T_ref)))
    )

end

function calculateIdealEnthalpy(T::Float64)
    #=
    Calculate the ideal gas enthalpy of each component
    h^ig(T) = h_ref(T_ref) + int_Tref^T cp dt
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]    
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        h:      enthalpy of each component [J mol^{-1}]
                Array{Float64,1}
    =#
    h = h_ref + integrateCp(T)
end

function calculateIdealEntropy(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    Calculate the ideal gas entropy of each component.
    s^ig(T,V,n) = s_ref(T_ref; pref) + int_Tref^T cp/t dt - R*log(nRT/Vpref)
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]    
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        s:      entropy of each component [J K^{-1} mol^{-1}]
                Array{Float64,1}
    =#
    s = s_ref + integrateCpByT(T) - R*log(n*R*T/(V*p_ref))
end

function calculateIdealPotential(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    Calculate ideal gas chemical potential of each component
    mu^ig(T,V,n) = h^ig(T) - T*s^ig(T,V,n)
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        mu:     chemical potential of each component [J mol^{-1}]
                Array{Float64,1}
    =#
    h   = calculateIdealEnthalpy(T)
    s   = calculateIdealEntropy(T, V, n)
    mu  = h - T*s
end

function idealHelmholtz(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    Calculate the ideal gas Helmholtz' energy.
    A^ig(T,V,n) = -NRT + n'*mu^ig(T,V,n)
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        A:      Helmholtz' energy [J]
                Float64
    =#
    mu  = calculateIdealPotential(T, V, n)
    A   = -sum(n)*R*T + dot(n, mu)
end

# Derivatives

function Aig_T(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    The derivative of Helmholtz' energy with respect to temperature
    dA^ig/dT = -S^ig
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        A_T:    dA^ig/dT [J K^{-1}]
                Float64
    =#
    s   = calculateIdealEntropy(T, V, n)
    A_T = -dot(n, s)
end

function Aig_V(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    The derivative of Helmholtz' energy with respect to volume
    dA^ig/dV = -p^ig
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        A_V:    dA^ig/dV [J m^{-3}]
                Float64
    =#
    A_V = -sum(n)*R*T/V
end

function Aig_n(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    The derivative of Helmholtz' energy with respect to mole vector
    dA^ig/dn = mu^ig
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        A_n:    dA^ig/dn [J mol^{-1}]
                Array{Float64,1}
    =#
    A_n = calculateIdealPotential(T, V, n)
end

# Second derivatives
function Aig_T_T(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    The second derivative of Helmholtz' energy with respect to temperature
    d2A^ig/dTdT = -n'*(cp/R - R/T)
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        A_T_T:  d2A^ig/dTdT [J K^{-2}]
                Float64
    =#
    A_T_T = -dot(n,(cp/R - R/T))
end

function Aig_T_V(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    The second derivative of Helmholtz' energy with respect to temperature and
    volume
    d2A^ig/dTdV = -NR/V
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        A_T_V:  d2A^ig/dTdV [J K^{-1} m^{-3}]
                Float64
    =#
    A_T_V = -sum(n)*R/V
end

function Aig_T_n(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    The second derivative of Helmholtz' energy with respect to temperature and
    mole vector
    d2A^ig/dTdn = -s + R
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        A_T_n:  d2A^ig/dTdn [J K^{-1} mol^{-1}]
                Array{Float64,1}
    =#
    s       = calculateIdealEntropy(T, V, n)
    A_T_n   = -s + R
end

function Aig_V_V(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    The second derivative of Helmholtz' energy with respect to volume
    d2A^ig/dVdV = NRT/V^2
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        A_V_V:  d2A^ig/dVdV [J m^{-6}]
                Float64
    =#
    A_V_V =sum(n)*R*T/V^2
end

function Aig_n_V(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    The second derivative of Helmholtz' energy with respect to volume and mole
    vector
    d2A^ig/dVdn = -eRT/V
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        A_n_V:  d2A^ig/dVdV [J m^{-3} mol^{-1}]
                Array{Float64,1}
    =#
    A_n_V = -e_n*R*T/V
end

function Aig_n_n(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    The second derivative of Helmholtz' energy with respect to mole vector
    d2A^ig/dndn = diag(n)^{-1}RT
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        A_n_n:  d2A^ig/dndn [J mol^{-2}]
                Diagonal{Float64}
    =#
    A_n_n = Diagonal(-1./n)*R*T
end

# Return the Hessian matrix of Aig with respect to volume and moles

function idealHessian(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    The Hessian of the ideal gas Helmholtz' energy.
    H^ig = [d2A^ig/dVdV     d2A^ig/dVdn'
            d2A^ig/dndV     d2A^ig/dndn]
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        H:      H^ig
                Array{Float64,2}
    =#
    A_V_V   = Aig_V_V(T, V, n)
    A_n_V   = Aig_n_V(T, V, n)
    A_n_n   = Aig_n_n(T, V, n)
    H       = [A_V_V A_n_V'
               A_n_V A_n_n]
end
