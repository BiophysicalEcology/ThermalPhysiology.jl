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
# _kelvin_param: struct field init — normalises to Quantity(K), units required.
# _to_kelvin:    Quantity(K) — for arithmetic against Quantity(K) struct fields
#                (same-unit ratios auto-cancel to bare Float64 under Unitful).
# _K:            bare Float64 K — for arithmetic with non-Quantity fields.
# All three treat bare Real input as °C (user-facing convention).

# Canonical Kelvin quantity type used for all Arrhenius-family struct fields.
const _KelvinQuantity = typeof(1.0u"K")

_to_kelvin(T::Unitful.Temperature) = uconvert(u"K", T)
_to_kelvin(T::Real)                = (T + 273.15) * u"K"   # assume °C

_K(T::Unitful.Temperature) = ustrip(u"K", uconvert(u"K", T))
_K(T::Real)                = T + 273.15

# Struct constructor parameter: Unitful temperature → Quantity(K). Requires units
# (no bare-Real fallback) to rule out silent K/°C/eV mix-ups.
_kelvin_param(x::Unitful.Temperature) = _K(x) * u"K"
_kelvin_param(x::Real) = throw(ArgumentError(
    "Arrhenius model temperature parameters must carry units, e.g. `9000.0u\"K\"` " *
    "(got bare value $x)."
))

# Bare Float64 in °C
_C(T::Unitful.Temperature) = ustrip(u"°C", uconvert(u"°C", T))
_C(T::Real)                = T   # already °C

# Activation parameter (T_A, T_AL, T_AH): temperature or energy, either way with units.
_arrhenius_temperature(x::Unitful.Temperature) = _kelvin_param(x)
_arrhenius_temperature(x::Unitful.Energy) = _kelvin_param(ea_to_ta(x))
_arrhenius_temperature(x::Real) = throw(ArgumentError(
    "Arrhenius activation parameters must carry units — a temperature " *
    "(e.g. `9000.0u\"K\"`) or an energy (e.g. `0.65u\"eV\"`); got bare value $x."
))
