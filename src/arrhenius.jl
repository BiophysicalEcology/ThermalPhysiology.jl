# ── Arrhenius-family temperature-correction models ────────────────────────────
#
# These models return a dimensionless correction factor (= 1 at T_ref).
# Temperature fields are Quantity(K) (_KelvinQuantity); same-unit ratios auto-cancel
# to bare Float64 under Unitful, so Tk must come from _to_kelvin (not _K) to match.
# Kooijman (2000 §2.6): calling Arrhenius "mechanistic" for organisms overstates
# the case — T_A is preferred over E_a to avoid false mechanistic connotations.

"""
    temperature_correction(model, T) → Float64

Return the dimensionless temperature-correction factor at temperature `T` for
an [`AbstractArrheniusModel`](@ref). The factor equals 1.0 at the model's
reference temperature `T_ref`.

`T` may be a Unitful temperature quantity or a bare `Float64` (°C assumed).
"""
function temperature_correction end

# ── ArrheniusModel ────────────────────────────────────────────────────────────

"""
    ArrheniusModel(; T_A, T_ref)

One-parameter temperature-correction model. Returns a dimensionless correction
factor equal to 1 at `T_ref`.

    temperature_correction(m, T) = exp(T_A/T_ref - T_A/T)

`T_A`, `T_ref` must be given as Unitful temperatures. No defaults — every parameter
is biologically meaningful and must be supplied deliberately. Example:

    ArrheniusModel(T_A=9000.0u"K", T_ref=293.15u"K")

`ArrheniusModel(E_A::Unitful.Energy; T_ref)` is also available, converting an
activation energy (e.g. `0.65u"eV"`) via [`ea_to_ta`](@ref).
"""
struct ArrheniusModel <: AbstractArrheniusModel
    T_A::_KelvinQuantity   # Arrhenius temperature; T_A = E_a / k_B
    T_ref::_KelvinQuantity # reference temperature; correction = 1 here
end

temperature_correction(m::ArrheniusModel, T) = exp(m.T_A / m.T_ref - m.T_A / _to_kelvin(T))
(m::ArrheniusModel)(T) = temperature_correction(m, T)

ArrheniusModel(E_A::Unitful.Energy; T_ref) = ArrheniusModel(T_A=ea_to_ta(E_A), T_ref=T_ref)
ArrheniusModel(x::Real; kwargs...) = throw(ArgumentError(
    "ArrheniusModel requires a Unitful quantity — a temperature (T_A, e.g. `9000.0u\"K\"`) " *
    "or an energy (E_A, e.g. `0.65u\"eV\"`); got bare value $x."
))

# ── SharpSchoolHighModel ──────────────────────────────────────────────────────

"""
    SharpSchoolHighModel(; T_A, T_ref, T_H, T_AH, rate_at_reference)

Sharpe-Schoolfield model with high-temperature enzyme deactivation.

    rate(T) = rate_at_reference * exp(T_A/T_ref - T_A/T) /
              (1 + exp(T_AH/T_H - T_AH/T))

Formula from Schoolfield, Sharpe & Magnuson (1981). Example:

    SharpSchoolHighModel(T_A=9000.0u"K", T_ref=293.15u"K", T_H=318.15u"K",
                          T_AH=90000.0u"K", rate_at_reference=1.0u"d^-1")
"""
struct SharpSchoolHighModel{R} <: AbstractArrheniusModel
    T_A::_KelvinQuantity
    T_ref::_KelvinQuantity
    T_H::_KelvinQuantity     # high-temperature transition
    T_AH::_KelvinQuantity    # Arrhenius temperature for high deactivation
    rate_at_reference::R     # Real or Unitful rate, preserved as-is
end

function temperature_correction(m::SharpSchoolHighModel, T)
    Tk = _to_kelvin(T)
    m.rate_at_reference *
        exp(m.T_A / m.T_ref - m.T_A / Tk) /
        (1 + exp(m.T_AH / m.T_H - m.T_AH / Tk))
end
(m::SharpSchoolHighModel)(T) = temperature_correction(m, T)

"""
    sharpe_schoolfield_high(; activation, reference_temperature, high_temperature,
                              high_deactivation, rate_at_reference)

Named constructor for [`SharpSchoolHighModel`](@ref). `activation` and
`high_deactivation` each accept either a Unitful temperature (`9000.0u"K"`) or
energy (`0.65u"eV"`, converted via [`ea_to_ta`](@ref)). Example:

    sharpe_schoolfield_high(activation=0.65u"eV", reference_temperature=293.15u"K",
                             high_temperature=318.15u"K", high_deactivation=10.0u"eV",
                             rate_at_reference=1.0u"d^-1")
"""
function sharpe_schoolfield_high(;
    activation,
    reference_temperature,
    high_temperature,
    high_deactivation,
    rate_at_reference,
)
    SharpSchoolHighModel(
        T_A  = _arrhenius_temperature(activation),
        T_ref = reference_temperature,
        T_H  = high_temperature,
        T_AH = _arrhenius_temperature(high_deactivation),
        rate_at_reference = rate_at_reference,
    )
end

# ── SharpSchoolLowModel ───────────────────────────────────────────────────────

"""
    SharpSchoolLowModel(; T_A, T_ref, T_L, T_AL, rate_at_reference)

Sharpe-Schoolfield model with low-temperature enzyme suppression. Example:

    SharpSchoolLowModel(T_A=9000.0u"K", T_ref=293.15u"K", T_L=273.15u"K",
                         T_AL=50000.0u"K", rate_at_reference=1.0u"d^-1")
"""
struct SharpSchoolLowModel{R} <: AbstractArrheniusModel
    T_A::_KelvinQuantity
    T_ref::_KelvinQuantity
    T_L::_KelvinQuantity     # low-temperature transition
    T_AL::_KelvinQuantity    # Arrhenius temperature for low deactivation
    rate_at_reference::R     # Real or Unitful rate, preserved as-is
end

function temperature_correction(m::SharpSchoolLowModel, T)
    Tk = _to_kelvin(T)
    m.rate_at_reference *
        exp(m.T_A / m.T_ref - m.T_A / Tk) /
        (1 + exp(m.T_AL / Tk - m.T_AL / m.T_L))
end
(m::SharpSchoolLowModel)(T) = temperature_correction(m, T)

"""
    sharpe_schoolfield_low(; activation, reference_temperature, low_temperature,
                             low_deactivation, rate_at_reference)

Named constructor for [`SharpSchoolLowModel`](@ref). `activation` and
`low_deactivation` each accept either a Unitful temperature or energy (see
[`sharpe_schoolfield_high`](@ref)). Example:

    sharpe_schoolfield_low(activation=0.65u"eV", reference_temperature=293.15u"K",
                            low_temperature=277.15u"K", low_deactivation=5.0u"eV",
                            rate_at_reference=1.0u"d^-1")
"""
function sharpe_schoolfield_low(;
    activation,
    reference_temperature,
    low_temperature,
    low_deactivation,
    rate_at_reference,
)
    SharpSchoolLowModel(
        T_A  = _arrhenius_temperature(activation),
        T_ref = reference_temperature,
        T_L  = low_temperature,
        T_AL = _arrhenius_temperature(low_deactivation),
        rate_at_reference = rate_at_reference,
    )
end

# ── SharpSchoolFullModel ──────────────────────────────────────────────────────

"""
    SharpSchoolFullModel(; T_A, T_ref, T_L, T_AL, T_H, T_AH, rate_at_reference)

Schoolfield (1981) full model with both low- and high-temperature enzyme deactivation.
`rate_at_reference` is the rate in the **absence** of enzyme inactivation at T_ref
(i.e. assumes T_ref is in the central Arrhenius zone).

    rate(T) = rate_at_reference * (T/T_ref) * exp(T_A/T_ref - T_A/T) /
              (1 + exp(T_AL*(1/T - 1/T_L)) + exp(T_AH*(1/T_H - 1/T)))

The `(T/T_ref)` pre-factor comes from collision-frequency theory and is part of the
original Schoolfield et al. (1981) Eq. (4); it is also present in rTPC's
`sharpeschoolfull_1981.R`. At biological temperatures it accounts for ~5–10% variation.

Low term `exp(T_AL*(1/T - 1/T_L))` is large at T < T_L (cold suppression);
high term `exp(T_AH*(1/T_H - 1/T))` is large at T > T_H (heat denaturation).

For the DEBtool-normalised variant (no T/T_ref, actual rate at T_ref guaranteed),
use [`SharpSchoolDEBModel`](@ref).

References: Schoolfield, Sharpe & Magnuson (1981) J Theor Biol 88:719;
`sharpeschoolfull_1981.R` in rTPC; `ArrFunc5` in Rezende `dynamic.landscape1.R`.

Example:

    SharpSchoolFullModel(T_A=9000.0u"K", T_ref=293.15u"K", T_L=273.15u"K",
                          T_AL=50000.0u"K", T_H=303.15u"K", T_AH=90000.0u"K",
                          rate_at_reference=1.0u"d^-1")
"""
struct SharpSchoolFullModel{R} <: AbstractArrheniusModel
    T_A::_KelvinQuantity
    T_ref::_KelvinQuantity
    T_L::_KelvinQuantity
    T_AL::_KelvinQuantity
    T_H::_KelvinQuantity
    T_AH::_KelvinQuantity
    rate_at_reference::R     # Real or Unitful rate, preserved as-is
end

function temperature_correction(m::SharpSchoolFullModel, T)
    Tk = _to_kelvin(T)
    boltzmann    = (Tk / m.T_ref) * exp(m.T_A / m.T_ref - m.T_A / Tk)
    low_term     = exp(m.T_AL / Tk - m.T_AL / m.T_L)    # large at T < T_L (cold suppression)
    high_term    = exp(m.T_AH / m.T_H - m.T_AH / Tk)    # large at T > T_H (heat denaturation)
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

Use this variant when fitting to rate data collected at a known reference temperature
or when computing constant temperature equivalents within DEB.

Reference: Kooijman (2010) DEB Theory §2.6; DEBtool_J `tempcorr.m`.

Example:

    SharpSchoolDEBModel(T_A=9000.0u"K", T_ref=293.15u"K", T_L=273.15u"K",
                         T_AL=50000.0u"K", T_H=303.15u"K", T_AH=90000.0u"K",
                         rate_at_reference=1.0u"d^-1")
"""
struct SharpSchoolDEBModel{R} <: AbstractArrheniusModel
    T_A::_KelvinQuantity
    T_ref::_KelvinQuantity
    T_L::_KelvinQuantity
    T_AL::_KelvinQuantity
    T_H::_KelvinQuantity
    T_AH::_KelvinQuantity
    rate_at_reference::R     # Real or Unitful rate, preserved as-is
end

function temperature_correction(m::SharpSchoolDEBModel, T)
    Tk = _to_kelvin(T)
    boltzmann    = exp(m.T_A / m.T_ref - m.T_A / Tk)
    low_ref      = exp(m.T_AL / m.T_ref - m.T_AL / m.T_L)
    high_ref     = exp(m.T_AH / m.T_H - m.T_AH / m.T_ref)
    low_T        = exp(m.T_AL / Tk     - m.T_AL / m.T_L)
    high_T       = exp(m.T_AH / m.T_H  - m.T_AH / Tk)
    norm         = (1 + low_ref + high_ref) / (1 + low_T + high_T)
    m.rate_at_reference * boltzmann * norm
end
(m::SharpSchoolDEBModel)(T) = temperature_correction(m, T)

"""
    sharpe_schoolfield_deb(; activation, reference_temperature,
                             low_temperature, low_deactivation,
                             high_temperature, high_deactivation,
                             rate_at_reference)

Named constructor for [`SharpSchoolDEBModel`](@ref). `activation`, `low_deactivation`
and `high_deactivation` each accept either a Unitful temperature or energy (see
[`sharpe_schoolfield_high`](@ref)). `rate_at_reference` is the actual observed rate
at `reference_temperature`. Example:

    sharpe_schoolfield_deb(activation=0.65u"eV", reference_temperature=293.15u"K",
                            low_temperature=277.15u"K", low_deactivation=5.0u"eV",
                            high_temperature=318.15u"K", high_deactivation=10.0u"eV",
                            rate_at_reference=1.0u"d^-1")
"""
function sharpe_schoolfield_deb(;
    activation,
    reference_temperature,
    low_temperature,
    low_deactivation,
    high_temperature,
    high_deactivation,
    rate_at_reference,
)
    SharpSchoolDEBModel(
        T_A  = _arrhenius_temperature(activation),
        T_ref = reference_temperature,
        T_L  = low_temperature,
        T_AL = _arrhenius_temperature(low_deactivation),
        T_H  = high_temperature,
        T_AH = _arrhenius_temperature(high_deactivation),
        rate_at_reference = rate_at_reference,
    )
end

"""
    sharpe_schoolfield(; activation, reference_temperature,
                         low_temperature, low_deactivation,
                         high_temperature, high_deactivation,
                         rate_at_reference)

Named constructor for [`SharpSchoolFullModel`](@ref). `activation`, `low_deactivation`
and `high_deactivation` each accept either a Unitful temperature or energy (see
[`sharpe_schoolfield_high`](@ref)). Example:

    sharpe_schoolfield(activation=0.65u"eV", reference_temperature=293.15u"K",
                        low_temperature=277.15u"K", low_deactivation=5.0u"eV",
                        high_temperature=318.15u"K", high_deactivation=10.0u"eV",
                        rate_at_reference=1.0u"d^-1")
"""
function sharpe_schoolfield(;
    activation,
    reference_temperature,
    low_temperature,
    low_deactivation,
    high_temperature,
    high_deactivation,
    rate_at_reference,
)
    SharpSchoolFullModel(
        T_A  = _arrhenius_temperature(activation),
        T_ref = reference_temperature,
        T_L  = low_temperature,
        T_AL = _arrhenius_temperature(low_deactivation),
        T_H  = high_temperature,
        T_AH = _arrhenius_temperature(high_deactivation),
        rate_at_reference = rate_at_reference,
    )
end

# ── JohnsonLewinModel ─────────────────────────────────────────────────────────

"""
    JohnsonLewinModel(; T_A, T_ref, T_H, T_AH, rate_at_reference)

Original Johnson-Lewin (1946) enzyme-kinetics model. Equivalent to the
high-deactivation Sharpe-Schoolfield form but historically distinct. Example:

    JohnsonLewinModel(T_A=9000.0u"K", T_ref=293.15u"K", T_H=318.15u"K",
                       T_AH=90000.0u"K", rate_at_reference=1.0u"d^-1")
"""
struct JohnsonLewinModel{R} <: AbstractArrheniusModel
    T_A::_KelvinQuantity
    T_ref::_KelvinQuantity
    T_H::_KelvinQuantity
    T_AH::_KelvinQuantity
    rate_at_reference::R     # Real or Unitful rate, preserved as-is
end

function temperature_correction(m::JohnsonLewinModel, T)
    Tk = _to_kelvin(T)
    m.rate_at_reference *
        exp(m.T_A / m.T_ref - m.T_A / Tk) /
        (1 + exp(m.T_AH / m.T_H - m.T_AH / Tk))
end
(m::JohnsonLewinModel)(T) = temperature_correction(m, T)

"""
    johnson_lewin(; activation, reference_temperature, high_temperature,
                    high_deactivation, rate_at_reference)

Named constructor for [`JohnsonLewinModel`](@ref). `activation` and
`high_deactivation` each accept either a Unitful temperature or energy (see
[`sharpe_schoolfield_high`](@ref)). Example:

    johnson_lewin(activation=0.65u"eV", reference_temperature=293.15u"K",
                  high_temperature=318.15u"K", high_deactivation=10.0u"eV",
                  rate_at_reference=1.0u"d^-1")
"""
function johnson_lewin(;
    activation,
    reference_temperature,
    high_temperature,
    high_deactivation,
    rate_at_reference,
)
    JohnsonLewinModel(
        T_A  = _arrhenius_temperature(activation),
        T_ref = reference_temperature,
        T_H  = high_temperature,
        T_AH = _arrhenius_temperature(high_deactivation),
        rate_at_reference = rate_at_reference,
    )
end

# ── Unitful-accepting keyword constructors ────────────────────────────────────
# Temperature params require units (_kelvin_param throws on bare numbers).
# rate_at_reference stays flexible (Real or Unitful rate), passed through as-is.
# No defaults — every parameter is biologically meaningful (see docstring examples).

ArrheniusModel(; T_A, T_ref) =
    ArrheniusModel(_kelvin_param(T_A), _kelvin_param(T_ref))

SharpSchoolHighModel(; T_A, T_ref, T_H, T_AH, rate_at_reference) = SharpSchoolHighModel(
    _kelvin_param(T_A), _kelvin_param(T_ref),
    _kelvin_param(T_H), _kelvin_param(T_AH), rate_at_reference,
)

SharpSchoolLowModel(; T_A, T_ref, T_L, T_AL, rate_at_reference) = SharpSchoolLowModel(
    _kelvin_param(T_A), _kelvin_param(T_ref),
    _kelvin_param(T_L), _kelvin_param(T_AL), rate_at_reference,
)

SharpSchoolFullModel(; T_A, T_ref, T_L, T_AL, T_H, T_AH, rate_at_reference) = SharpSchoolFullModel(
    _kelvin_param(T_A), _kelvin_param(T_ref),
    _kelvin_param(T_L), _kelvin_param(T_AL),
    _kelvin_param(T_H), _kelvin_param(T_AH), rate_at_reference,
)

SharpSchoolDEBModel(; T_A, T_ref, T_L, T_AL, T_H, T_AH, rate_at_reference) = SharpSchoolDEBModel(
    _kelvin_param(T_A), _kelvin_param(T_ref),
    _kelvin_param(T_L), _kelvin_param(T_AL),
    _kelvin_param(T_H), _kelvin_param(T_AH), rate_at_reference,
)

JohnsonLewinModel(; T_A, T_ref, T_H, T_AH, rate_at_reference) = JohnsonLewinModel(
    _kelvin_param(T_A), _kelvin_param(T_ref),
    _kelvin_param(T_H), _kelvin_param(T_AH), rate_at_reference,
)
