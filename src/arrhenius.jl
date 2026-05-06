# ── Arrhenius-family temperature-correction models ────────────────────────────
#
# These models return a dimensionless correction factor (= 1 at T_ref).
# All temperatures stored as bare Float64 in Kelvin (DEBtool convention).
# Kooijman (2000 §2.6): calling Arrhenius "mechanistic" for organisms overstates
# the case — T_A is preferred over E_a to avoid false mechanistic connotations.

# ── ArrheniusModel ────────────────────────────────────────────────────────────

"""
    ArrheniusModel(; T_A, T_ref)

One-parameter temperature-correction model. Returns a dimensionless correction
factor equal to 1 at `T_ref`.

    temperature_correction(m, T) = exp(T_A/T_ref - T_A/T)

Parameters stored in Kelvin. Convenience constructor `arrhenius` accepts Unitful
quantities or bare values.
"""
@kwdef struct ArrheniusModel <: AbstractArrheniusModel
    T_A::Float64   = 8000.0    # Arrhenius temperature (K); T_A = E_a / k_B
    T_ref::Float64 = 293.15    # reference temperature (K); correction = 1 here
end

temperature_correction(m::ArrheniusModel, T) = exp(m.T_A / m.T_ref - m.T_A / _K(T))
(m::ArrheniusModel)(T) = temperature_correction(m, T)

"""
    arrhenius(; activation_energy=0.65u"eV", reference_temperature=293.15u"K")

Named constructor for [`ArrheniusModel`](@ref). Accepts Unitful eV/K quantities
or bare Float64 values (assumed eV and K respectively).
"""
function arrhenius(;
    activation_energy    = 0.65u"eV",
    reference_temperature = 293.15u"K",
)
    ArrheniusModel(T_A=_K(ea_to_ta(activation_energy)), T_ref=_K(reference_temperature))
end

# ── SharpSchoolHighModel ──────────────────────────────────────────────────────

"""
    SharpSchoolHighModel(; T_A, T_ref, T_H, T_AH, rate_at_reference)

Sharpe-Schoolfield model with high-temperature enzyme deactivation.

    rate(T) = rate_at_reference * exp(T_A/T_ref - T_A/T) /
              (1 + exp(T_AH/T_H - T_AH/T))

Formula from Schoolfield, Sharpe & Magnuson (1981).
"""
@kwdef struct SharpSchoolHighModel <: AbstractArrheniusModel
    T_A::Float64             = 8000.0
    T_ref::Float64           = 293.15
    T_H::Float64             = 318.15   # high-temperature transition (K)
    T_AH::Float64            = 90000.0  # Arrhenius temperature for high deactivation (K)
    rate_at_reference::Float64 = 1.0
end

function temperature_correction(m::SharpSchoolHighModel, T)
    Tk = _K(T)
    m.rate_at_reference *
        exp(m.T_A / m.T_ref - m.T_A / Tk) /
        (1 + exp(m.T_AH / m.T_H - m.T_AH / Tk))
end
(m::SharpSchoolHighModel)(T) = temperature_correction(m, T)

function sharpe_schoolfield_high(;
    activation_energy       = 0.65u"eV",
    reference_temperature    = 293.15u"K",
    high_temperature         = 318.15u"K",
    high_deactivation_energy = 10.0u"eV",
    rate_at_reference        = 1.0,
)
    SharpSchoolHighModel(
        T_A  = _K(ea_to_ta(activation_energy)),
        T_ref = _K(reference_temperature),
        T_H  = _K(high_temperature),
        T_AH = _K(ea_to_ta(high_deactivation_energy)),
        rate_at_reference = rate_at_reference,
    )
end

# ── SharpSchoolLowModel ───────────────────────────────────────────────────────

"""
    SharpSchoolLowModel(; T_A, T_ref, T_L, T_AL, rate_at_reference)

Sharpe-Schoolfield model with low-temperature enzyme suppression.
"""
@kwdef struct SharpSchoolLowModel <: AbstractArrheniusModel
    T_A::Float64             = 8000.0
    T_ref::Float64           = 293.15
    T_L::Float64             = 277.15   # low-temperature transition (K)
    T_AL::Float64            = 50000.0  # Arrhenius temperature for low deactivation (K)
    rate_at_reference::Float64 = 1.0
end

function temperature_correction(m::SharpSchoolLowModel, T)
    Tk = _K(T)
    m.rate_at_reference *
        exp(m.T_A / m.T_ref - m.T_A / Tk) /
        (1 + exp(m.T_AL / Tk - m.T_AL / m.T_L))
end
(m::SharpSchoolLowModel)(T) = temperature_correction(m, T)

function sharpe_schoolfield_low(;
    activation_energy      = 0.65u"eV",
    reference_temperature   = 293.15u"K",
    low_temperature         = 277.15u"K",
    low_deactivation_energy = 5.0u"eV",
    rate_at_reference       = 1.0,
)
    SharpSchoolLowModel(
        T_A  = _K(ea_to_ta(activation_energy)),
        T_ref = _K(reference_temperature),
        T_L  = _K(low_temperature),
        T_AL = _K(ea_to_ta(low_deactivation_energy)),
        rate_at_reference = rate_at_reference,
    )
end

# ── SharpSchoolFullModel ──────────────────────────────────────────────────────

"""
    SharpSchoolFullModel(; T_A, T_ref, T_L, T_AL, T_H, T_AH, rate_at_reference)

Schoolfield (1981) full model with both low- and high-temperature enzyme deactivation.
`rate_at_reference` is the rate in the **absence** of enzyme inactivation at T_ref
(i.e. assumes T_ref is in the central Arrhenius zone).

    rate(T) = rate_at_reference * exp(T_A/T_ref - T_A/T) /
              (1 + exp(T_AL*(1/T - 1/T_L)) + exp(T_AH*(1/T_H - 1/T)))

Low term `exp(T_AL*(1/T - 1/T_L))` is large at T < T_L (cold suppression);
high term `exp(T_AH*(1/T_H - 1/T))` is large at T > T_H (heat denaturation).

For the DEBtool-normalised variant where `rate_at_reference` is the **actual observed**
rate at T_ref (correction applied in numerator too), use [`SharpSchoolDEBModel`](@ref).

References: Schoolfield, Sharpe & Magnuson (1981) J Theor Biol 88:719;
`ArrFunc5` in Rezende `dynamic.landscape1.R`; `sharpeschoolfull_1981.R` in rTPC.
"""
@kwdef struct SharpSchoolFullModel <: AbstractArrheniusModel
    T_A::Float64             = 8000.0
    T_ref::Float64           = 293.15
    T_L::Float64             = 277.15
    T_AL::Float64            = 50000.0
    T_H::Float64             = 318.15
    T_AH::Float64            = 90000.0
    rate_at_reference::Float64 = 1.0
end

function temperature_correction(m::SharpSchoolFullModel, T)
    Tk = _K(T)
    boltzmann    = exp(m.T_A / m.T_ref - m.T_A / Tk)
    low_term     = exp(m.T_AL / Tk - m.T_AL / m.T_L)    # large at T < T_L (cold suppression)
    high_term    = exp(m.T_AH / m.T_H - m.T_AH / Tk)   # large at T > T_H (heat denaturation)
    inactivation = 1 / (1 + low_term + high_term)
    m.rate_at_reference * boltzmann * inactivation
end
(m::SharpSchoolFullModel)(T) = temperature_correction(m, T)

# ── SharpSchoolDEBModel ───────────────────────────────────────────────────────

"""
    SharpSchoolDEBModel(; T_A, T_ref, T_L, T_AL, T_H, T_AH, rate_at_reference)

DEBtool-normalised Sharpe-Schoolfield model. `rate_at_reference` is the **actual
observed rate** at T_ref — the correction is applied in both numerator and denominator,
so `temperature_correction(m, T_ref) == rate_at_reference` for any T_ref.

    rate(T) = rate_at_reference * exp(T_A/T_ref - T_A/T) *
              (1 + exp(T_AL/T_ref - T_AL/T_L) + exp(T_AH/T_H - T_AH/T_ref)) /
              (1 + exp(T_AL/T   - T_AL/T_L)   + exp(T_AH/T_H - T_AH/T))

This is `tempcorr` in DEBtool_J and the form used in NicheMapR for CTE calculations.
Use this variant when fitting to rate data collected at a known reference temperature
or when computing constant temperature equivalents within DEB.

Reference: Kooijman (2010) DEB Theory §2.6; DEBtool_J `tempcorr.m`.
"""
@kwdef struct SharpSchoolDEBModel <: AbstractArrheniusModel
    T_A::Float64             = 8000.0
    T_ref::Float64           = 293.15
    T_L::Float64             = 277.15
    T_AL::Float64            = 50000.0
    T_H::Float64             = 318.15
    T_AH::Float64            = 90000.0
    rate_at_reference::Float64 = 1.0
end

function temperature_correction(m::SharpSchoolDEBModel, T)
    Tk   = _K(T)
    Tref = m.T_ref
    boltzmann    = exp(m.T_A / Tref - m.T_A / Tk)
    low_ref      = exp(m.T_AL / Tref - m.T_AL / m.T_L)
    high_ref     = exp(m.T_AH / m.T_H - m.T_AH / Tref)
    low_T        = exp(m.T_AL / Tk   - m.T_AL / m.T_L)
    high_T       = exp(m.T_AH / m.T_H - m.T_AH / Tk)
    norm         = (1 + low_ref + high_ref) / (1 + low_T + high_T)
    m.rate_at_reference * boltzmann * norm
end
(m::SharpSchoolDEBModel)(T) = temperature_correction(m, T)

"""
    sharpe_schoolfield_deb(; activation_energy, reference_temperature,
                             low_temperature, low_deactivation_energy,
                             high_temperature, high_deactivation_energy,
                             rate_at_reference)

Named constructor for [`SharpSchoolDEBModel`](@ref).
`rate_at_reference` is the actual observed rate at `reference_temperature`.
"""
function sharpe_schoolfield_deb(;
    activation_energy       = 0.65u"eV",
    reference_temperature    = 293.15u"K",
    low_temperature         = 277.15u"K",
    low_deactivation_energy  = 5.0u"eV",
    high_temperature        = 318.15u"K",
    high_deactivation_energy = 10.0u"eV",
    rate_at_reference        = 1.0,
)
    SharpSchoolDEBModel(
        T_A  = _K(ea_to_ta(activation_energy)),
        T_ref = _K(reference_temperature),
        T_L  = _K(low_temperature),
        T_AL = _K(ea_to_ta(low_deactivation_energy)),
        T_H  = _K(high_temperature),
        T_AH = _K(ea_to_ta(high_deactivation_energy)),
        rate_at_reference = rate_at_reference,
    )
end

function sharpe_schoolfield(;
    activation_energy       = 0.65u"eV",
    reference_temperature    = 293.15u"K",
    low_temperature         = 277.15u"K",
    low_deactivation_energy  = 5.0u"eV",
    high_temperature        = 318.15u"K",
    high_deactivation_energy = 10.0u"eV",
    rate_at_reference        = 1.0,
)
    SharpSchoolFullModel(
        T_A  = _K(ea_to_ta(activation_energy)),
        T_ref = _K(reference_temperature),
        T_L  = _K(low_temperature),
        T_AL = _K(ea_to_ta(low_deactivation_energy)),
        T_H  = _K(high_temperature),
        T_AH = _K(ea_to_ta(high_deactivation_energy)),
        rate_at_reference = rate_at_reference,
    )
end

# ── JohnsonLewinModel ─────────────────────────────────────────────────────────

"""
    JohnsonLewinModel(; T_A, T_ref, T_H, T_AH, rate_at_reference)

Original Johnson-Lewin (1946) enzyme-kinetics model. Equivalent to the
high-deactivation Sharpe-Schoolfield form but historically distinct.
"""
@kwdef struct JohnsonLewinModel <: AbstractArrheniusModel
    T_A::Float64             = 8000.0
    T_ref::Float64           = 293.15
    T_H::Float64             = 318.15
    T_AH::Float64            = 90000.0
    rate_at_reference::Float64 = 1.0
end

function temperature_correction(m::JohnsonLewinModel, T)
    Tk = _K(T)
    m.rate_at_reference *
        exp(m.T_A / m.T_ref - m.T_A / Tk) /
        (1 + exp(m.T_AH / m.T_H - m.T_AH / Tk))
end
(m::JohnsonLewinModel)(T) = temperature_correction(m, T)

function johnson_lewin(;
    activation_energy       = 0.65u"eV",
    reference_temperature    = 293.15u"K",
    high_temperature        = 318.15u"K",
    high_deactivation_energy = 10.0u"eV",
    rate_at_reference        = 1.0,
)
    JohnsonLewinModel(
        T_A  = _K(ea_to_ta(activation_energy)),
        T_ref = _K(reference_temperature),
        T_H  = _K(high_temperature),
        T_AH = _K(ea_to_ta(high_deactivation_energy)),
        rate_at_reference = rate_at_reference,
    )
end
