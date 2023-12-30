using PyPlot

# Standard constants
g = 9.81 # m/s^2 (TODO: replace with standard gravity)
R = 286.0 #  8.31446261815324?
γ = 1.4 # calorically perfect air

# See Lecture 19, Page 3 for numerical example with these numbers
cp = 1004.0 # J/kg*K
# specific enthalpy of fuel (air mixture?), I think
h = 4.3e7 # [J/kg]

# Lecture 20, Equation 1
Γ = √(γ) * (2.0 / (γ + 1.0)) ^ ((γ + 1.0) / (2.0 * (γ - 1.0)))
println("Γ = ", Γ)

# Climb conditions
θ = 20.0 # deg 
M0 = 0.6
alt = 6.0e3 # m

# Ambient conditions
T0 = 261.0 # K
# 486 mb = 48,600 Pa
P0 = 4.86e4 # Pa

# Aircraft params
L_over_D = 15.0
m = 2.e3 # kg

################################################################################
######################### (a) Thrust Calculation ###############################
################################################################################

# m * dh/dt = L - mgcos(θ) = 0
L = m * g * cosd(θ)
println("L = ", L)
D = L / L_over_D
println("D = ", D)

# m * du/dt = F - D - mgsinθ = 0
F = D + m * g * sind(θ)
println("F = ", F)
println()

################################################################################
####################### (b) Temperature & Pressure #############################
################################################################################

# Speed of sound
a0 = √(γ * R * T0)
println("a0 = ", a0)

# Engine design points
M2 = 0.4
Tt4 = 1200 # K

θ0 = (1.0 + (γ - 1.0) / 2.0 * M0 ^ 2.0)
δ0 = (1.0 + (γ - 1.0) / 2.0 * M0 ^ 2.0) ^ (γ / (γ - 1.0))
println("θ0 = ", θ0)
println("δ0 = ", δ0)

Tt0 = θ0 * T0
v = Tt4 / T0
θt = v
println("θt = ", θt)
println()

# Temperatures and pressures ratios for each component

# Compressor design (optimal τc)
τc = √(θt) / θ0 # Lecture 19, Page 2
println("τc = ", τc)

πc = τc ^ (γ / (γ - 1.0))
println("πc = ", πc)

# Shaft balance to get τt
τt = 1 - (τc - 1) / v
println("τt = ", τt)

πt = τt ^ (γ / (γ - 1.0))
println("πt = ", πt)
println()

# Total temperatures and pressures at each station

M2 = 0.4
θ2 = (1.0 + (γ - 1.0) / 2.0 * M2 ^ 2.0)
δ2 = (1.0 + (γ - 1.0) / 2.0 * M2 ^ 2.0) ^ (γ / (γ - 1.0))

Tt2 = θ2 * T0
Pt2 = δ2 * P0
println("Tt2 = ", Tt2, "\tPt2 = ", Pt2)

Tt3 = √(T0 * Tt4) # for optimal compressor temperature ratio
Pt3 = πc * Pt2
println("Tt3 = ", Tt3, "\t\tPt3 = ", Pt3)

# Tt4 = designed to above
Pt4 = Pt3 # ideally
println("Tt4 = ", Tt4, "\t\t\tPt4 = ", Pt4)

Tt5 = τt * Tt4
Pt5 = πt * Pt4
println("Tt5 = ", Tt5, "\t\tPt5 = ", Pt5)

# "For non-afterburning turbojet..." (Lecture 20, Equation 4)
Pt7 = Pt5
Tt7 = Tt5
println("Tt7 = ", Tt7, "\t\tPt7 = ", Pt7)

# Lecture 18, Eqautions 4 & 5
Tte = T0 * θt * τt
Pte = P0 * δ0 * πc * πt # assume πb ≈ 1
println("Tte = ", Tte, "\t\tPte = ", Pte)
println()

################################################################################
######################### (c) Air Flow Calculation #############################
################################################################################

# Lecture 18, Equation 10
trm1 = 2.0 / (γ - 1.0)
trm2 = θt / (θ0 * τc)
trm3 = θt - θ0 * (τc - 1.0) - trm2
eqn10_rhs = √(trm1 * trm3) - M0
mdot = F / (a0 * eqn10_rhs) # no bleeds, so same for each station
println("mdot =  ", mdot)

# Double check with max thrust per unit airflow calc
F_per_unit_airflow = F / (mdot * a0)
println("F_per_unit_airflow = ", F_per_unit_airflow)

# Lecture 19, Equation 11
trm4 = (√(θt) - 1.0) ^ 2.0
F_per_unit_airflow_2 = √(trm1 * trm4 + M0 ^ 2.0) - M0
println("F_per_unit_airflow_2 = ", F_per_unit_airflow_2)
println()

################################################################################
############################ (d) Key Flow Areas ################################
################################################################################

# mdot = Γ * Pti * Ai / √(R * Tti)

# Turbine inlet area
A4 = mdot * √(R * Tt4)/ (Γ * Pt4)
println("A4 = ", A4)

num = (γ + 1.0) / 2.0
den = 1.0 + (γ - 1.0) / 2.0 * M2 ^ 2.0
exp = (γ + 1.0) / (2.0 * (γ - 1.0))
mbar2 = M2 * (num / den) ^ exp
println("mbar2 = ", mbar2)

# Lecture 20, Equation 9
# mbar2 = τc ^ (γ / (γ - 1.0)) / √(v) * A4 / A2
A2 = τc ^ (γ / (γ - 1.0)) / √(v) * A4 / mbar2
println("A2 = ", A2)

# Double check A2
# Lecture 20, Equation 8
# mdot = Γ * Pt2 * A2 / √(R * Tt2) * mbar2
A2_2 = mdot * √(R * Tt2) / (Γ * Pt2 * mbar2)
println("A2_2 = ", A2_2)

# Nozzle throat area (Lecture 20, Equation 5)
# τt = (A4 / A7) ^ (2.0 * (γ - 1.0) / (γ + 1.0))
A7 = A4 * τt ^ (-1.0 * (γ + 1.0) / (2.0 * (γ - 1.0)))
println("A7 = ", A7)

# # Double check exponent math
A7_2 = √(τt) * A4 / πt
println("A7_2 = ", A7_2)

# A = πr^2 = πD^2/4
d7 = √(4.0 * A7 / π)
println("d7 = ", d7)
println()

################################################################################
############################ (e) Fuel Flow & ISP ###############################
################################################################################

mdotf = mdot * cp * T0 / h * (θt - θ0 * τc)
println("mdotf = ", mdotf)

mdotf_2 = mdot * cp * Tt2 * (v - τc) / h
println("mdotf_2 = ", mdotf_2)

mdotf_3 = mdot * cp * (Tt4 - Tt3) / h
println("mdotf_3 = ", mdotf_3)
println()

f = mdotf / mdot
println("f = ", f)

trm5 = h * a0 / (g * cp * T0)
println("trm5 = ", trm5)

Isp = trm5 * (F / (mdot * a0)) / (θt - θ0 * τc)
println("Isp = ", Isp)

Isp_2 = trm5 * (F / (mdot * a0)) / (θt - √(θt))
println("Isp_2 = ", Isp_2)

println()

################################################################################
################## (f) Exhaust Velocity & Efficiencies #########################
################################################################################

u0 = M0 * a0
println("u0 = ", u0)

# Lecture 18, Equation 6
Me = √((2.0 / (γ - 1.0)) * (θ0 * τc * τt - 1.0))
println("Me = ", Me)

ue = Me * a0
println("ue = ", ue)

η_prop = F * u0 / (mdot * (ue ^ 2.0 / 2.0 - u0 ^ 2.0 / 2.0))
println("η_prop = ", η_prop)

η_prop_2 = 2.0 / (F / (mdot * a0 * M0) + 2.0)
println("η_prop_2 = ", η_prop_2)

η_therm = mdot * (ue ^ 2.0 / 2.0 - u0 ^ 2.0 / 2.0) / (mdotf * h)
println("η_therm = ", η_therm)

η_all = F * u0 / (mdotf * h)
println("η_all = ", η_all)

η_all_2 = η_therm * η_prop
println("η_all_2 = ", η_all_2)

η_all_3 = g * u0 * Isp / h
println("η_all_3 = ", η_all_3)

println()

################################################################################
###################### Plot compressor ratio sweeps ############################
################################################################################

# find mbar2, M2, θ, f for sweep of πc (also do f for M0 = 0.0)

close("all")

figure()

Πc = Vector{Float64}(4.0:0.5:14.0)
println(Πc)
println("Πc = ", Πc)

Τc = Πc .^ ((γ - 1.0) / γ)
println("Τc = ", Τc)

plt.plot(Πc, Τc)
# plt.label("")

Τt = 1.0 .- (Τc .- 1.0) ./ v
plt.plot(Πc, Τt)

mbar2v = Πc ./ v .* (A4 / A2)
plt.plot(Πc, mbar2v)

figure()
fv = cp * Tt0 / h .* (v .- Τc)
plt.plot(Πc, fv)

figure()
plt.plot(mbar2v ./ mbar2, Πc)
plt.plot(mbar2v ./ mbar2v[end], Πc)

return nothing
