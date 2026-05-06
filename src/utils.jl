using Unitful

# Boltzmann constant
const k_B = 8.617333e-5u"eV/K"

# ── Unit conversions ───────────────────────────────────────────────────────────

ea_to_ta(E_a) = E_a / k_B   # activation energy (eV) → Arrhenius temperature (K)
ta_to_ea(T_A) = T_A * k_B   # Arrhenius temperature (K) → activation energy (eV)

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
