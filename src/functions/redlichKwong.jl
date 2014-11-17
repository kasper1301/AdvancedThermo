#=
    Functions for calculating properties for the Redlich-Kwong equation
    of state. The variables aretemperature, volume and number of moles. 

    Needs the following global variables:
        R       :   The ideal gas constant [J K^{-1} mol^{-1}]  Float64
        T_ref   :   The reference temperature                   Float64
        p_ref   :   The reference pressure                      Float64
        h_ref   :   The reference enthalpies                    Array{Float64,1}
        s_ref   :   The reference entropies                     Array{Float64,1}
        cp      :   The heat capacities                         Array{Float64,1}
        e_n     :   A unit vector with the same size as n       Array{Int64,1}
        a       :   Redlich-Kwong constant sqrt(a)              Array{Float64,1}
        b       :   Redlich-Kwong constant b                    Array{Float64,1}

    Requires the idealGas functions to calculate the properties

    Functions:
        redlichKwongB               :   Hard sphere volume B
        redlichKwongA               :   Redlich-Kwong parameter A
        residualHelmholtz           :   Residual Helmholtz energy
        residualEntropy             :   Residual entropy
        residualChemicalPotential   :   Residual chemical potential
        Helmholtz                   :   Helmholtz energy
        Enthalpy                    :   Enthalpy
        Entropy                     :   Entropy
        chemicalPotential           :   Chemical potential
        pressure                    :   pressure
        residualHessian             :   d2A^rv/d(V,n)d(V,n)
        Hessian                     :   d2A/d(V,n)d(V,n)
        Arv_T_T                     :   d2A^rv/dT2
        Arv_T_V                     :   d2A^rv/dTdV
        Arv_T_n                     :   d2A^rv/dTdn
        Arv_V_V                     :   d2A^rv/dVdV
        Arv_n_V                     :   d2A^rv/dndV
        Arv_n_V                     :   d2A^rv/dndV
        Arv_n_n                     :   d2A^rv/dndn
=#

# Parameters A and B

function redlichKwongB(n::Array{Float64,1})
    #=
    Calculate the parameter B in the Redlich-Kwong equation of state.
    B(n) = b'*n
    Input:
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        B:      parameter B [m^3]
                Float64
    =#
    B = dot(b,n)
end

function redlichKwongA(T::Float64, n::Array{Float64,1})
    #=
    Calculate the parameter A in the Redlich-Kwong equation of state
    A(T,n) = n'*(a*a')*n/sqrt(T)
    Input:
        T:      temperature [K]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        A:      parameter A [J m^3]
                Float64
    =#
    A = dot(n,(a*a')*n)/sqrt(T)
end

# Residuals

function residualHelmholtz(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    Calculate the residual Helmholtz' energy
    A^rv(T,V,n) = NRTlog(V/(V-B)) + A/B*log(V/(V+B)
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        Arv:    residual Helmholtz' energy [J]
                Float64
    =#
    A   = redlichKwongA(T, n)
    B   = redlichKwongB(n)
    Arv = sum(n)*R*T*log(V/(V-B)) + A/B*log(V/(V+B))
end

function residualEntropy(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    Calculate residual entropy
    S^rv(T,V,n) = NRlog((V-B)/B)+A_T/Blog((V+B)/V)
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        Srv:    residual entropy [J K^{-1}]
                Float64
    =#
    A   = redlichKwongA(T, n)
    A_T = -A/(2T)
    B   = redlichKwongB(n)
    Srv = sum(n)*R*log((V-B)/B) + A_T/B*log((V+B)/V)
end

function residualChemicalPotential(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    Calculate residual chemical potential
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        mu:     residual chemical potential [J mol^{-1}]
                Array{Float64,1}
    =#
    A   = redlichKwongA(T, n)
    B   = redlichKwongB(n)
    A_n = 2(a*a')*n/sqrt(T)
    mu  = e_n*R*T*log(V/(V-B)) + sum(n)*R*T*b/(V-B) +
          (A_n*B-A*b)/B^2*log(V/(V+B)) - A*b/(B*(V+B))
end

# Properties

function Helmholtz(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    Calculate the Helmholtz' energy
    A(T,V,n) = A^rv + A^ig
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
    A = residualHelmholtz(T, V, n) + idealHelmholtz(T, V, n)
end

function Entropy(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    Calculate the entropy
    S(T,V,n) = S^rv+ S^ig
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        S:      entropy [J K^{-1}]
                Float64
    =#
    S = residualEntropy(T, V, n) - Aig_T(T, V, n)
end

function Enthalpy(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    Calculate the enthalpy
    H(T,V,n) = A + TS + pV
    =#
    H = Helmholtz(T, V, n) + T*Entropy(T, V, n) + V*pressure(T, V, n)
end

function chemicalPotential(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    Calculate the chemical potential
    mu(T,V,n) = mu^rv+ mu^ig
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        mu:     chemical potential [J mol^{-1}]
                Array{Float64,1}
    =#
    mu = residualChemicalPotential(T, V, n) + calculateIdealPotential(T, V, n)
end 

function pressure(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    Calculate the pressure
    p^RK(T,V,n) = NRT/(V-B) - A/(V(V+B))
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        p:      pressure [Pa]
                Float64
    =#
    A = redlichKwongA(T, n)
    B = redlichKwongB(n)
    p = sum(n)*R*T/(V-B) - A/(V*(V+B))
end

# Second derivatives

function Arv_T_T(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    Calculate the second derivative of the Helmholtz' energy residual with 
    respect to temperature
    d2A^rv/dTdT = 3A/(4T^2*B)*log(V/(V+B))
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        A_T_T:  d2A^rv/dTdT [J K^{-2}]
                Float64
    =#
    A     = redlichKwongA(T, n)
    B     = redlichKwongB(n)
    A_T_T = 3A/(4T^2*B)*log(V/(V+B))
end

function Arv_T_V(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    Calculate the second derivative of the Helmholtz' energy residual with 
    respect to temperature and volume
    d2A^rv/dTdV = -(A/2T)/(V(V+B)) - NRB/(V(V-B))
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        A_T_V:  d2A^rv/dTdV [J K^{-1} m^{-3}]
                Float64
    =#
    A     = redlichKwongA(T, n)
    B     = redlichKwongB(n)
    A_T_V = -(A/2T)/(V*(V+B)) - sum(n)*R*B/(V*(V-B))
end

function Arv_T_n(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    Calculate the second derivative of the Helmholtz' energy residual with 
    respect to temperature and mole vector
    d2A^rv/dTdn = eRlog(V/(V-B)) + NRb/(V-B) + (A_nT B-A_T B_n)/B^2log(V/(V+B))-
                  A_T B_n/(B(V+B))
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
    A     = redlichKwongA(T, n)
    B     = redlichKwongB(n)
    A_n   = 2(a*a')*n/sqrt(T)
    A_T   = -A/2T
    A_nT  = -A_n/2T
    A_T_n = e_n*R*log(V/(V-B)) + sum(n)*R*b/(V-B) +
            (A_nT*B-A_T*b)/B^2*log(V/(V+B)) - A_T*b/(B*(B+V))
end

function Arv_V_V(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    Calculate the second derivative of the Helmholtz' energy residual with 
    respect to volume
    d2A^rv/dVdV = -NRT/V^2 + NRT/(V-B)^2 - A/(V^2(V+B)) - A/(V(V+B)^2)
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
    A     = redlichKwongA(T,n)
    B     = redlichKwongB(n)
    A_V_V = sum(n)*(-R*T/V^2 + R*T/(V-B)^2) - A/(V^2*(V+B)) - A/(V*(V+B)^2)
end

function Arv_n_V(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    Calculate the second derivative of the Helmholtz' energy residual with 
    respect to volume and mole vector
    d2A^rv/dndV = eRT/V - (eRT(V-B)+NRTb)/(V-B)^2 + (A_n(V+B)-Ab)/(V(V+B)^2)
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        A_n_V:  d2A^ig/dndV [J m^{-3} mol^{-1}]
                Array{Float64,1}
    =#
    A     = redlichKwongA(T,n)
    B     = redlichKwongB(n)
    A_n   = 2(a*a')*n/sqrt(T)
    A_n_V = e_n*R*T/V - (e_n*R*T*(V-B)+sum(n)*R*T*b)/(V-B)^2 + 
            (A_n*(V+B) - A*b)/(V*(V+B)^2)
end

function Arv_n_n(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    Calculate the second derivative of the Helmholtz' energy residual with 
    respect to mole vector
    d2A^rv/dndn = RT/(V-B)*(eb'+be') + NRTbb'/(V-B)^2 + 
                  (B^2A_nn-(A_n b'+bA_n')B+2Abb')/B^3log(V/(V+B)) - 
                  (A_nB-Ab)/B^2*b'/(V+B) - b/(V+B)*(A_n'B-b'A)/B^2 + 
                  A/B*(bb')/(V+B)^2
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        A_n_n:  d2A^ig/dndn [J mol^{-2}]
                Symmetric{Float64}
    =#
    A     = redlichKwongA(T,n)
    B     = redlichKwongB(n)
    A_n   = 2(a*a')*n/sqrt(T)
    A_nn  = 2*(a*a')/sqrt(T)
    A_n_n = R*T/(V-B)*(e_n*b'+b*e_n') + sum(n)*R*T*b*b'/(V-B)^2 +
            (B^2*A_nn - (A_n*b'+b*A_n')*B+2A*b*b')/B^3*log(V/(V+B)) - 
            (A_n*B-A*b)*b'/(B^2*(V+B)) - b*(B*A_n' - A*b')/(B^2*(V+B)) +
            A*(b*b')/(B*(V+B)^2)
end

# Hessian matrices

function residualHessian(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    The Hessian of the residual Helmholtz' energy.
    H^rv = [d2A^rv/dVdV     d2A^rv/dVdn'
            d2A^rv/dndV     d2A^rv/dndn]
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        H:      H^rv
                Array{Float64,2}
    =#
    A_V_V   = Arv_V_V(T, V, n)
    A_n_V   = Arv_n_V(T, V, n)
    A_n_n   = Arv_n_n(T, V, n)
    H       = [A_V_V A_n_V'
               A_n_V A_n_n]
end

function Hessian(T::Float64, V::Float64, n::Array{Float64,1})
    #=
    The Hessian of the Helmholtz' energy.
    H = H^rv + H^ig
    Input:
        T:      temperature [K]
                Float64
        V:      volume [m^3]
                Float64
        n:      mole vector [mol]
                Array{Float64,1}
    Output:
        H:      H
                Array{Float64,2}
    =#
    H = idealHessian(T, V, n) + residualHessian(T, V, n)
end
