include("includes.jl")

n_tot = [1.60; 94.50; 2.60; 0.81; 0.52]/100.03                     # Composition

n_v = 0.5n_tot
n_l = n_tot - n_v

V_v = logspace(-1.0, log10(1.01dot(n_v,b)), 100)
V_l = logspace(-1.0, log10(1.01dot(n_l,b)), 100)
T = 180.0
p_v = zeros(length(V_v))
for i in 1:length(V_v)
    p_v[i] = pressure(T, V_v[i], n_v)
end

fig = plt.figure()

semilogx(V_l, p_v, "g")

V_l = 1.1dot(n_l, b)
V_v = 3e-4 - V_l
x_v = equilibrium(T, [V_v, n_v], [3e-4, n_tot])
n_v = x_v[2:end]
n_l = n_tot - n_v
V_v = logspace(-1.0, log10(1.01dot(n_v,b)), 100)
V_l = logspace(-1.0, log10(1.01dot(n_l,b)), 100)
p_v = zeros(length(V_v))
p_l = copy(p_v)

for i in 1:length(V_v)
    p_v[i] = pressure(T, V_v[i], n_v)
    p_l[i] = pressure(T, V_l[i], n_l)
end

semilogx(V_l, p_l, "b")
semilogx(V_v, p_v, "r")
semilogx(V_l, ones(length(V_l))*pressure(T, x_v[1], x_v[2:end]),"k:")

legend([L"$p^\mathrm{l,v},\quad\mathbf{n}^\mathrm{l,v} = \frac{1}{2}\mathbf{n}$",
        L"$p^\mathrm{l},\quad\mathbf{n}^\mathrm{l} = \mathbf{n}^\mathrm{l,eq}$",
        L"$p^\mathrm{v},\quad\mathbf{n}^\mathrm{v} =  \mathbf{n}^\mathrm{v,eq}$",
        L"$p^\mathrm{eq}=\SI{2.78e6}{\pascal}$"])
xlabel(L"$V\quad\left[\SI{}{\cubic\meter}\right]$")
ylabel(L"$p\quad\left[\SI{}{\pascal}\right]$")
printfig(fig, filename="../fig/pV.tex")
ylim([0, 1e7])
printfig(fig, filename="../fig/pVzoom.tex")