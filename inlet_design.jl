using PyPlot

plot_stuffz = true

γ = 1.4 # calorically perfect air
# Lecture 20, Equation 1
Γ = √(γ) * (2.0 / (γ + 1.0)) ^ ((γ + 1.0) / (2.0 * (γ - 1.0)))
println("Γ = ", Γ)

M0 = 1.5 # 3.0
M2 = 0.5

δ0 = (1.0 + (γ - 1.0) / 2.0 * M0 ^ 2.0) ^ (γ / (γ - 1.0))
println("δ0 = ", δ0)

θ0 = 1.0 + (γ - 1.0) / 2.0 * M0 ^ 2.0
println("θ0 = ", θ0)

################################################################################
######################### (a) Calc Area Ratios #################################
################################################################################
println()

function Astar_A_ratio_fn(γ, M)
    exp = (γ + 1.0) / (2.0 * (γ - 1.0))
    num = (γ + 1.0) / 2.0
    den = 1.0 + (γ - 1.0) / 2.0 * M ^ 2.0
    return M * (num / den) ^ exp
end

At_Astar = 1 / Astar_A_ratio_fn(γ, 1.0)
println("At_Astar = ", At_Astar)

A1_Astar = 1 / Astar_A_ratio_fn(γ, M0)
println("A1_Astar = ", A1_Astar)

A2_Astar = 1 / Astar_A_ratio_fn(γ, M2)
println("A2_Astar = ", A2_Astar)

At_A2 = At_Astar / A2_Astar
A1_A2 = A1_Astar / A2_Astar

println("At_A2 = ", At_A2)
println("A1_A2 = ", A1_A2)

################################################################################
######################### (b) Normal Shock Calcs ###############################
################################################################################
println()
M0 = 1.25

num = 1.0 + ((γ - 1.0) / 2.0) * M0 ^ 2.0
den = γ * M0 ^ 2.0 - (γ - 1.0) / 2.0
M1 = √(num / den)

println("M1 = ", M1)

trm1 = 2.0 * γ * M0 ^ 2.0 / (γ + 1.0)
trm2 = (γ - 1.0) / (γ + 1.0)
trm3 = (trm1 - trm2) ^ (1.0 / (γ - 1.0))
trm4 = (γ - 1.0) * M0 ^ 2.0 + 2.0
trm5 = (γ + 1.0) * M0 ^ 2.0
trm6 = (trm4 / trm5) ^ (γ / (γ - 1.0))

Pt0_Pt1 = trm3 * trm6
P_loss = 1 / Pt0_Pt1
println("P_loss = ", P_loss)

A1_Astar_schock = 1 / Astar_A_ratio_fn(γ, M1)
println("A1_Astar_schock = ", A1_Astar_schock)

f_spill = 1.0 - A1_Astar_schock / A1_Astar
println("f_spill = ", f_spill)

################################################################################
######################### (c) Shock Re-Position ################################
################################################################################
println()

# not really Atnew / A* anymore, but compared to the original ratio
Atnew_Astar = A1_Astar * At_Astar / A1_Astar_schock

Atnew_A2 = Atnew_Astar / A2_Astar
println("Atnew_A2 = ", Atnew_A2)

# A2_Atnew = 1.0 / Atnew_A2
# println("A2_Atnew = ", A2_Atnew)

# Would be sick af to do bisection search here (on section that is monotonically decreasing...)

M2new = .61
Atnew_A2_check = Astar_A_ratio_fn(γ, M2new)
println("Atnew_A2_check = ", Atnew_A2_check)
println("M2new = ", M2new)

################################################################################
######################### (d) Shock "Pops In" ##################################
################################################################################
println()

Mtnew = 1.126

Atnew_A1 = Atnew_Astar / A1_Astar

# Forrealz this time
Atnew_A1_check = Astar_A_ratio_fn(γ, M0) / Astar_A_ratio_fn(γ, Mtnew)

println("Mtnew = ", Mtnew)
println("Atnew_A1 = ", Atnew_A1)
println("Atnew_A1_check = ", Atnew_A1_check)

################################################################################
######################### (e) Thrust Calculation ###############################
################################################################################
println()

M0 = Mtnew

trm1 = 2.0 * γ * M0 ^ 2.0 / (γ + 1.0)
trm2 = (γ - 1.0) / (γ + 1.0)
trm3 = (trm1 - trm2) ^ (1.0 / (γ - 1.0))
trm4 = (γ - 1.0) * M0 ^ 2.0 + 2.0
trm5 = (γ + 1.0) * M0 ^ 2.0
trm6 = (trm4 / trm5) ^ (γ / (γ - 1.0))

Pt0_Pt1 = trm3 * trm6
println("Pt0_Pt1 = ", Pt0_Pt1)
P_loss = 1 / Pt0_Pt1
println("P_loss = ", P_loss)

Astar_Atnew = Astar_A_ratio_fn(γ, M1)
println("Astar_Atnew = ", Astar_Atnew)
# Atnewnew_A2 = 

Astar_A2new = Astar_Atnew * Atnew_A2
println("Astar_A2new = ", Astar_A2new)

M2new = 0.57771
Astar_A2new_check = Astar_A_ratio_fn(γ, M2new)
println("Astar_A2new_check = ", Astar_A2new_check)
println("M2new = ", M2new)



################################################################################
######################### Nozzle and Diffuser Plots ############################
################################################################################

close("all")

# Re-create plot from Lecture 6

if plot_stuffz

    γv = [1.3, 1.4, 1.66]
    Mv = Vector(0.0:0.1:6.0)

    Astar_A_ratio = zeros(length(γv), length(Mv))

    for (i, γi) in enumerate(γv)
        for (j, Mj) in enumerate(Mv)
            Astar_A_ratio[i, j] = Astar_A_ratio_fn(γi, Mj)
        end
    end

    plt.figure()
    plt.plot(Mv, Astar_A_ratio[1, :])
    plt.plot(Mv, Astar_A_ratio[2, :])
    plt.plot(Mv, Astar_A_ratio[3, :])

    # Re-create plot from Lecture 22

    Mv = Vector(0.2:0.05:3.0)
    A_Astar_ratio = zeros(length(Mv))

    for (i, Mi) in enumerate(Mv)
        A_Astar_ratio[i] = 1.0 / Astar_A_ratio_fn(γ, Mi)
    end

    plt.figure()
    plt.plot(Mv, A_Astar_ratio)
    plt.show()

end