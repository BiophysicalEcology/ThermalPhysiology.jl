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

# Struct constructor parameter: Unitful temperature → Quantity(K). Units required.
_kelvin_param(x::Unitful.Temperature) = _K(x) * u"K"
_kelvin_param(x::Real) = throw(ArgumentError(
    "Arrhenius model temperature parameters must carry units, e.g. `9000.0u\"K\"` " *
    "(got bare value $x)."
))

# Bare Float64 in °C
_C(T::Unitful.Temperature) = ustrip(u"°C", uconvert(u"°C", T))
_C(T::Real)                = T   # already °C

# Validates units without forcing K or °C -- preserves whichever the caller gave
# (mirrors _time_param; DryAirProperties etc. in FluidProperties.jl do the same).
_temperature_param(x::Unitful.Temperature) = x
_temperature_param(x::Real) = throw(ArgumentError(
    "Thermal death time temperature parameters must carry units, e.g. `39.0u\"°C\"` " *
    "(got bare value $x)."
))

# z_value: a difference, stored in K (Unitful rejects arithmetic on °C-typed
# quantities). °C input accepted too, reinterpreted directly (no offset).
_z_value_param(x::Unitful.Quantity) = Unitful.unit(x) in (u"K", u"°C") ? ustrip(x) * u"K" :
    throw(ArgumentError("z_value must be given in K or °C (a temperature " *
        "difference); got $x."))
_z_value_param(x::Real) = throw(ArgumentError(
    "z_value must carry units, e.g. `4.0u\"°C\"` or `4.0u\"K\"` (got bare value $x)."
))

# Activation parameter (T_A, T_AL, T_AH): temperature or energy, either way with units.
_arrhenius_temperature(x::Unitful.Temperature) = _kelvin_param(x)
_arrhenius_temperature(x::Unitful.Energy) = _kelvin_param(ea_to_ta(x))
_arrhenius_temperature(x::Real) = throw(ArgumentError(
    "Arrhenius activation parameters must carry units — a temperature " *
    "(e.g. `9000.0u\"K\"`) or an energy (e.g. `0.65u\"eV\"`); got bare value $x."
))

# ── Time helpers (internal, TDT) ──────────────────────────────────────────────

# Bare Float64 in minutes, for the handful of formulas (Rezende/Jørgensen) that
# are defined against a specific minutes-numeric convention.
_min(t::Unitful.Time) = ustrip(u"minute", uconvert(u"minute", t))

# survival_time returns Unitful time for LogLinearTDTModel, bare Float64 for
# ToleranceLandscape; generic AbstractTDTModel code normalises via this.
_min_or_bare(t::Unitful.Time) = _min(t)
_min_or_bare(t::Real) = t

# Validates units without forcing a specific one -- callers keep whichever time
# unit they passed (e.g. reference_duration in hours stays in hours), since
# Unitful correctly handles ratios/products across different time units.
_time_param(x::Unitful.Time) = x
_time_param(x::Real) = throw(ArgumentError(
    "Thermal death time parameters must carry units, e.g. `60.0u\"minute\"` " *
    "(got bare value $x)."
))

# Ramp-rate validator: temperature/time dimension (e.g. K/minute). °C/time is
# itself rejected by Unitful (AffineError) since °C can't be used in a rate.
const _RAMP_RATE_DIM = Unitful.dimension(u"K/minute")
_ramp_rate(x::Unitful.Quantity) = Unitful.dimension(x) == _RAMP_RATE_DIM ?
    ustrip(u"K/minute", uconvert(u"K/minute", x)) :
    throw(ArgumentError("ramp_rate must have dimensions of temperature/time, " *
        "e.g. `0.5u\"K/minute\"`; got $x."))
_ramp_rate(x::Real) = throw(ArgumentError(
    "ramp_rate must carry units, e.g. `0.5u\"K/minute\"` (got bare value $x)."
))
