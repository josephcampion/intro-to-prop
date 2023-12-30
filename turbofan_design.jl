
# Standard constants
g = 9.81 # m/s^2 (TODO: replace with standard gravity)
R = 287.0 #  8.31446261815324?
γ = 1.4 # calorically perfect air

# See Lecture 19, Page 3 for numerical example with these numbers
cp = 1004.0 # J/kg*K
# specific enthalpy of fuel (air mixture?), I think
h = 4.3e7 # [J/kg]

# Design conditions
M0 = 0.9
T0 = 220.0 # K
θt = 6.25
α = 6.0

# Speed of sound
a0 = √(γ * R * T0)
println("a0 = ", a0)

u0 = M0 * a0
println("u0 = ", u0)

θ0 = (1.0 + (γ - 1.0) / 2.0 * M0 ^ 2.0)

Tt4 = θt * T0
println("Tt4 = ", Tt4)

println()

################################################################################
############################ Assuming u6 = u8 ##################################
################################################################################
println("Assuming u6 = u8\n")

# Lecture 21, Page 4: (τc)_Fmax
τc = √(θt) / θ0
println("τc = ", τc)

πc = τc ^ (γ / (γ - 1.0))
println("πc = ", πc)

# Lecture 21, Equation 17b
# for maximum thrust
τf = (√θt - 1.0) ^ 2.0 / (θ0 * (1.0 + α)) + 1
println("τf = ", τf)

πf = τf ^ (γ / (γ - 1.0))

println("πf = ", πf)

# Lecture 21, Equation 15
τt = 1.0 - (θ0 / θt) * ((τc - 1.0) + α * (τf - 1.0))
println("τt = ", τt)

πt = τt ^ (γ / (γ - 1.0))
println("πt = ", πt)
println()

trm1 = 2.0 / (γ - 1.0)

F_mdota = (1.0 + α) * (√(trm1 * (θ0 * τf - 1.0)) - M0)

# trm2 = θt / (θ0 * τc)
# trm3 = θ0 * τc * τt - 1.0
# Fc_mdota = √(trm1 * trm2 * trm3) - M0

trm4 = θ0 * τf - 1.0
Fbp_mdota = α * (√(trm1 * trm4) - M0)

println("F_mdot0 = ", F_mdota)
# println("Fc_mdot0 = ", Fc_mdota)
println("Fc_mdot0 = ", F_mdota - Fbp_mdota)
println("Fbp_mdot0 = ", Fbp_mdota)
# println("F_mdot0_2 = ", Fc_mdota + Fbp_mdota)
println()

# Assuming P8 = P0 (Lecture 21, Page 2)
M8 = √(trm1 * (θ0 * τf - 1.0))
println("M8 = ", M8)

u8 = a0 * M8
println("u8 = ", u8)

u8_2 = u0 * (Fbp_mdota / (α * M0) + 1.0)
println("u8_2 = ", u8_2)
println()

u6 = u8

"""
Energy balance in combustor: 
mdotf * h = mdot * cp * (Tt4 - Tt3)

Divide by T0, and multiply by Tt2 / Tt2 
Remember (Tt0 = Tt2, assuming ideal inlet)

mdot * cp (Tt4 / T0 - Tt3 / T0 * Tt2 / Tt2)
mdot * cp (θt - τc * θ0)

f == mdotf / mdot 
"""
f = T0 / h * cp * (θt - τc * θ0)
println("f = ", f)

# Lecture 19, Equation 11 (still applies?)
trm5 = a0 * h / (g * cp * T0)
Isp = trm5 * F_mdota / (θt - √θt)
println("Isp = ", Isp)

Isp_2 = η_overall * h / (g * u0)
println("Isp_2 = ", Isp_2)
println()

η_prop = 2.0 / (u6 / u0 + 1.0)
println("η_prop = ", η_prop)
# ^ this seems to match graph on Page 5 of Lecture 21

η_overall = g * u0 * Isp / h
println("η_overall = ", η_overall)

η_therm = η_overall / η_prop
println("η_therm = ", η_therm)

# Need to include total mass for "power-in-jet" term:
η_therm_2 = (1.0 + α) * (u6 ^ 2.0 / 2.0 - u0 ^ 2.0 / 2.0) / (f * h)
println("η_therm_2 = ", η_therm_2)

println()

################################################################################
########################### Assuming u6 != u8 ##################################
################################################################################
println("Assuming u6 != u8")

πf = 2.0
println("πf = ", πf)
println()
τf = πf ^ ((γ - 1.0) / γ)
println("τf = ", τf)

# Assuming P8 = P0 (Lecture 21, Page 2)
M8 = √(trm1 * (θ0 * τf - 1.0))
println("M8 = ", M8)

u8 = a0 * M8
println("u8 = ", u8)
println()

# Assume same τc as before:
# Lecture 21, Equation 15
τt = 1.0 - (θ0 / θt) * ((τc - 1.0) + α * (τf - 1.0))
println("τt = ", τt)

# Lecture 18, Equation 7
trm6 = θ0 * τc * τt - 1.0
trm7 = θt / (θ0 * τc)
u6 = u0 / M0 * √(trm1 * trm6 * trm7)
println("u6 = ", u6)

trm4 = θ0 * τf - 1.0
Fbp_mdota = α * (√(trm1 * trm4) - M0)
println("Fbp_mdot0 = ", Fbp_mdota)

trm2 = θt / (θ0 * τc)
trm3 = θ0 * τc * τt - 1.0
Fc_mdota = √(trm1 * trm2 * trm3) - M0
println("Fc_mdota = ", Fc_mdota)

F_mdota = Fbp_mdota + Fc_mdota
println("F_mdota = ", F_mdota)
println()

# Calculate Efficiencies

# f stays the same because πc unchanged
# need to split up mass flows though
# Need to include total mass for "power-in-jet" term:
num1 = (u6 ^ 2.0 / 2.0 - u0 ^ 2.0 / 2.0)
num2 = α * (u8 ^ 2.0 / 2.0 - u0 ^ 2.0 / 2.0)

η_therm = (num1 + num2) / (f * h)
println("η_therm = ", η_therm)

"""
Thermal efficiency unchanged because enthalpy balance unchanged.
The same power coming from fuel (mdotf * h) is converted to change in temperature of air flow in combustor.
"""

"""
η_prop = F * u0 / ddt_sum_of_Ke (see above)

Divide both num and den by u0 * a0 * mdot:
"""

num = F_mdota * u0

trm1 = (u6 ^ 2.0 / 2.0 - u0 ^ 2.0 / 2.0)
trm2 = α * (u8 ^ 2.0 / 2.0 - u0 ^ 2.0 / 2.0)
den = (trm1 + trm2) / a0

η_prop = num / den
println("η_prop = ", η_prop)

η_overall = η_prop *  η_therm
println("η_overall = ", η_overall)

# η_overall = g * u0 * Isp / h
Isp = h * η_overall / (g * u0)
println("Isp = ", Isp)
