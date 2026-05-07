using Unitful

# Boltzmann constant
const k_B = 8.617333e-5u"eV/K"

# ── Unit conversions ───────────────────────────────────────────────────────────

"""
    ea_to_ta(E_a) → T_A

Convert activation energy in eV to Arrhenius temperature in K: `T_A = E_a / k_B`.
Accepts Unitful quantities or bare Float64 (assumed eV).
"""
ea_to_ta(E_a) = E_a / k_B

"""
    ta_to_ea(T_A) → E_a

Convert Arrhenius temperature in K to activation energy in eV: `E_a = T_A × k_B`.
Accepts Unitful quantities or bare Float64 (assumed K).
"""
ta_to_ea(T_A) = T_A * k_B

# ── Temperature helpers (internal) ────────────────────────────────────────────

# Convert anything to Kelvin (Unitful)
_to_kelvin(T::Unitful.Temperature) = uconvert(u"K", T)
_to_kelvin(T::Real)                = (T + 273.15) * u"K"   # assume °C

# Bare Float64 in Kelvin for arithmetic inside model equations
_K(T::Unitful.Temperature) = ustrip(u"K", uconvert(u"K", T))
_K(T::Real)                = T + 273.15

# Bare Float64 in °C
_C(T::Unitful.Temperature) = ustrip(u"°C", uconvert(u"°C", T))
_C(T::Real)                = T   # already °C
