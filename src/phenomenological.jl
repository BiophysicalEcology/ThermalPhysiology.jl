# ── Phenomenological TPC models ───────────────────────────────────────────────
#
# Temperature parameters require Unitful units (K or °C); internal formulas work
# in bare °C via _C, or bare K via _K, on already-validated fields. "Rate" fields
# (maximum_rate, rate_constant, rate_at_reference, maximum_performance) are plain
# multiplicative prefactors and accept either a bare Real (dimensionless/relative
# scale) or a Unitful rate of whatever process the model represents — preserved
# as given, never normalised. Internal symbols (T_opt, CTmax, rmax) follow
# literature convention.

"""
    thermal_performance(model, T) → Float64 or corrected rate

Return the thermal performance (absolute or relative rate) at temperature `T`
for an [`AbstractPhenomenologicalModel`](@ref). `T` must be a Unitful temperature
quantity (e.g. `25.0u"°C"`); bare numbers are rejected.
"""
function thermal_performance end

# ── Universal TPC (Arnoldi, Jackson, Peralta-Maraver & Payne 2025, PNAS) ──────

"""
    UniversalTPCModel(; T_opt, E, maximum_performance)

Universal thermal performance curve (Arnoldi et al. 2025, PNAS).
Two-parameter model: `y(x) = exp(x) * (1 - x)` where `x = (T - T_opt) / E`.

- `T_opt`: optimal temperature — Unitful temperature; bare rejected
- `E`: thermal breadth — Unitful K or °C difference; bare rejected. Related to
  Arrhenius temperature by `E = T_ref² / T_A`
- `maximum_performance`: peak rate at T_opt

The UTPC thermal breadth connects to TDT z-value: `z = E * log(10)`. Example:

    UniversalTPCModel(T_opt=303.15u"K", E=10.0u"K", maximum_performance=1.0)
"""
struct UniversalTPCModel{C,P} <: AbstractPhenomenologicalModel
    T_opt::C
    E::_KelvinQuantity
    maximum_performance::P
end

UniversalTPCModel(; T_opt, E, maximum_performance) =
    UniversalTPCModel(_temperature_param(T_opt), _z_value_param(E), maximum_performance)

function thermal_performance(m::UniversalTPCModel, T)
    x = (_strict_C(T) - _C(m.T_opt)) / ustrip(m.E)
    m.maximum_performance * exp(x) * (1 - x)
end
(m::UniversalTPCModel)(T) = thermal_performance(m, T)

"""
    utpc(; optimal_temperature, thermal_breadth, maximum_performance)

Named constructor for [`UniversalTPCModel`](@ref) with friendlier keyword names.
Example:

    utpc(optimal_temperature=30.0u"°C", thermal_breadth=10.0u"K", maximum_performance=1.0)
"""
utpc(; optimal_temperature, thermal_breadth, maximum_performance) =
    UniversalTPCModel(T_opt=optimal_temperature, E=thermal_breadth, maximum_performance=maximum_performance)

# ── Deutsch 2008 ──────────────────────────────────────────────────────────────

"""
    DeutschModel(; maximum_rate, optimal_temperature, critical_thermal_maximum, width_parameter)

Modified Deutsch (2008, PNAS) thermal performance curve.
Gaussian rise below T_opt; quadratic decline above T_opt to CTmax.

    T < T_opt: rate = maximum_rate * exp(-((T - T_opt) / (2*width_parameter))^2)
    T > T_opt: rate = maximum_rate * (1 - ((T - T_opt) / (T_opt - CTmax))^2)

Temperature parameters require Unitful units. Example:

    DeutschModel(maximum_rate=1.0, optimal_temperature=25.0u"°C",
                 critical_thermal_maximum=40.0u"°C", width_parameter=5.0u"K")
"""
struct DeutschModel{R,C1,C2} <: AbstractPhenomenologicalModel
    maximum_rate::R
    optimal_temperature::C1
    critical_thermal_maximum::C2
    width_parameter::_KelvinQuantity
end

DeutschModel(; maximum_rate, optimal_temperature, critical_thermal_maximum, width_parameter) =
    DeutschModel(maximum_rate, _temperature_param(optimal_temperature),
                 _temperature_param(critical_thermal_maximum), _z_value_param(width_parameter))

function thermal_performance(m::DeutschModel, T)
    Tc, Topt = _strict_C(T), _C(m.optimal_temperature)
    if Tc < Topt
        m.maximum_rate * exp(-((Tc - Topt) / (2 * ustrip(m.width_parameter)))^2)
    else
        m.maximum_rate * (1 - ((Tc - Topt) / (Topt - _C(m.critical_thermal_maximum)))^2)
    end
end
(m::DeutschModel)(T) = thermal_performance(m, T)

deutsch(; kwargs...) = DeutschModel(; kwargs...)

# ── Briere 1 (1999) ───────────────────────────────────────────────────────────

"""
    Briere1Model(; rate_constant, minimum_temperature, maximum_temperature)

Briere 1 (1999) thermal performance curve.
`rate = rate_constant * T * (T - T_min) * sqrt(T_max - T)` for T_min < T < T_max.
Temperature parameters require Unitful units. Example:

    Briere1Model(rate_constant=0.01, minimum_temperature=10.0u"°C", maximum_temperature=40.0u"°C")
"""
struct Briere1Model{R,C1,C2} <: AbstractPhenomenologicalModel
    rate_constant::R
    minimum_temperature::C1
    maximum_temperature::C2
end

Briere1Model(; rate_constant, minimum_temperature, maximum_temperature) =
    Briere1Model(rate_constant, _temperature_param(minimum_temperature), _temperature_param(maximum_temperature))

function thermal_performance(m::Briere1Model, T)
    Tc = _strict_C(T)
    Tmin, Tmax = _C(m.minimum_temperature), _C(m.maximum_temperature)
    Tc <= Tmin || Tc >= Tmax && return zero(m.rate_constant)
    val = m.rate_constant * Tc * (Tc - Tmin) * sqrt(Tmax - Tc)
    max(zero(m.rate_constant), val)
end
(m::Briere1Model)(T) = thermal_performance(m, T)

briere(; kwargs...) = Briere1Model(; kwargs...)

# ── Briere 2 (1999) ───────────────────────────────────────────────────────────

"""
    Briere2Model(; rate_constant, minimum_temperature, maximum_temperature, shape_parameter)

Briere 2 (1999) — generalised form with shape exponent.
`rate = rate_constant * T * (T - T_min) * (T_max - T)^(1/shape_parameter)`
Temperature parameters require Unitful units. Example:

    Briere2Model(rate_constant=0.01, minimum_temperature=10.0u"°C",
                 maximum_temperature=40.0u"°C", shape_parameter=2.0)
"""
struct Briere2Model{R,C1,C2} <: AbstractPhenomenologicalModel
    rate_constant::R
    minimum_temperature::C1
    maximum_temperature::C2
    shape_parameter::Float64
end

Briere2Model(; rate_constant, minimum_temperature, maximum_temperature, shape_parameter) =
    Briere2Model(rate_constant, _temperature_param(minimum_temperature),
                 _temperature_param(maximum_temperature), Float64(shape_parameter))

function thermal_performance(m::Briere2Model, T)
    Tc = _strict_C(T)
    Tmin, Tmax = _C(m.minimum_temperature), _C(m.maximum_temperature)
    Tc <= Tmin || Tc >= Tmax && return zero(m.rate_constant)
    val = m.rate_constant * Tc * (Tc - Tmin) *
          (Tmax - Tc)^(1 / m.shape_parameter)
    max(zero(m.rate_constant), val)
end
(m::Briere2Model)(T) = thermal_performance(m, T)

# ── Gaussian ──────────────────────────────────────────────────────────────────

"""
    GaussianModel(; maximum_rate, optimal_temperature, width_parameter)

Symmetric Gaussian TPC.
`rate = maximum_rate * exp(-0.5 * ((T - T_opt) / width_parameter)^2)`
Temperature parameters require Unitful units. Example:

    GaussianModel(maximum_rate=1.0, optimal_temperature=25.0u"°C", width_parameter=5.0u"K")
"""
struct GaussianModel{R,C} <: AbstractPhenomenologicalModel
    maximum_rate::R
    optimal_temperature::C
    width_parameter::_KelvinQuantity
end

GaussianModel(; maximum_rate, optimal_temperature, width_parameter) =
    GaussianModel(maximum_rate, _temperature_param(optimal_temperature), _z_value_param(width_parameter))

function thermal_performance(m::GaussianModel, T)
    Tc = _strict_C(T)
    m.maximum_rate * exp(-0.5 * ((Tc - _C(m.optimal_temperature)) / ustrip(m.width_parameter))^2)
end
(m::GaussianModel)(T) = thermal_performance(m, T)

gaussian(; kwargs...) = GaussianModel(; kwargs...)

# ── Thomas 2012 ───────────────────────────────────────────────────────────────

"""
    Thomas2012Model(; rate_constant, shape_parameter, optimal_temperature)

Thomas et al. (2012) thermal performance curve.
`rate = rate_constant * (T - T_opt + b) * (T - T_opt - b) * (-1)`
where `b` (`shape_parameter`) is the half-width of the performance curve.
Returns 0 outside the performance range. Temperature parameters require Unitful
units. Example:

    Thomas2012Model(rate_constant=0.5, shape_parameter=15.0u"K", optimal_temperature=25.0u"°C")
"""
struct Thomas2012Model{R,C} <: AbstractPhenomenologicalModel
    rate_constant::R
    shape_parameter::_KelvinQuantity
    optimal_temperature::C
end

Thomas2012Model(; rate_constant, shape_parameter, optimal_temperature) =
    Thomas2012Model(rate_constant, _z_value_param(shape_parameter), _temperature_param(optimal_temperature))

function thermal_performance(m::Thomas2012Model, T)
    Tc = _strict_C(T)
    b  = ustrip(m.shape_parameter)
    T0 = _C(m.optimal_temperature)
    val = m.rate_constant * (Tc - T0 + b) * (b - (Tc - T0))
    max(zero(m.rate_constant), val)
end
(m::Thomas2012Model)(T) = thermal_performance(m, T)

thomas(; kwargs...) = Thomas2012Model(; kwargs...)

# ── Thomas 2017 ───────────────────────────────────────────────────────────────

"""
    Thomas2017Model(; maximum_rate, optimal_temperature, width_parameter, skewness)

Thomas et al. (2017) asymmetric TPC: skewed-Gaussian form. Temperature
parameters require Unitful units. Example:

    Thomas2017Model(maximum_rate=1.0, optimal_temperature=25.0u"°C",
                     width_parameter=5.0u"K", skewness=0.0)
"""
struct Thomas2017Model{R,C} <: AbstractPhenomenologicalModel
    maximum_rate::R
    optimal_temperature::C
    width_parameter::_KelvinQuantity
    skewness::Float64            # positive → right skew
end

Thomas2017Model(; maximum_rate, optimal_temperature, width_parameter, skewness) =
    Thomas2017Model(maximum_rate, _temperature_param(optimal_temperature),
                     _z_value_param(width_parameter), Float64(skewness))

function thermal_performance(m::Thomas2017Model, T)
    Tc = _strict_C(T)
    d  = Tc - _C(m.optimal_temperature)
    σ  = ustrip(m.width_parameter) * (1 + m.skewness * sign(d))
    m.maximum_rate * exp(-0.5 * (d / σ)^2)
end
(m::Thomas2017Model)(T) = thermal_performance(m, T)

# ── Pawar 2018 ────────────────────────────────────────────────────────────────

"""
    PawarModel(; rate_at_reference, activation_energy, deactivation_energy,
                 peak_temperature, reference_temperature)

Pawar et al. (2018) metabolic TPC: Arrhenius rise with high-T deactivation,
parameterised in biologically meaningful terms.

`activation_energy`/`deactivation_energy` each accept either a Unitful energy
(e.g. `0.65u"eV"`) or an Arrhenius temperature (e.g. `9000.0u"K"`, via
[`ea_to_ta`](@ref)); `peak_temperature`/`reference_temperature` require a
Unitful temperature. Example:

    PawarModel(rate_at_reference=1.0, activation_energy=0.65u"eV", deactivation_energy=1.15u"eV",
               peak_temperature=30.0u"°C", reference_temperature=20.0u"°C")
"""
struct PawarModel{R} <: AbstractPhenomenologicalModel
    rate_at_reference::R
    activation_energy::_KelvinQuantity     # Arrhenius-temperature equivalent (Ea/k_B)
    deactivation_energy::_KelvinQuantity
    peak_temperature::_KelvinQuantity
    reference_temperature::_KelvinQuantity
end

PawarModel(; rate_at_reference, activation_energy, deactivation_energy, peak_temperature, reference_temperature) =
    PawarModel(rate_at_reference, _arrhenius_temperature(activation_energy), _arrhenius_temperature(deactivation_energy),
               _kelvin_param(peak_temperature), _kelvin_param(reference_temperature))

function thermal_performance(m::PawarModel, T)
    Tk = _to_kelvin(_temperature_param(T))
    boltzmann  = exp(m.activation_energy / m.reference_temperature - m.activation_energy / Tk)
    correction = 1 / (1 + exp(m.deactivation_energy / m.peak_temperature - m.deactivation_energy / Tk))
    m.rate_at_reference * boltzmann * correction
end
(m::PawarModel)(T) = thermal_performance(m, T)

pawar(; kwargs...) = PawarModel(; kwargs...)

# ── Lactin 2 ──────────────────────────────────────────────────────────────────

"""
    Lactin2Model(; rate_constant, maximum_temperature, delta_temperature, intercept)

Lactin 2 (1995) thermal performance curve.
`rate = exp(rate_constant * T) - exp(rate_constant * T_max - (T_max - T) / delta_T) + intercept`

`rate_constant` and `intercept` are dimensionless curve-shape coefficients
(the Lactin equation is defined against the numeric magnitude of T in °C, not
a portable physical rate) and stay bare `Float64`; `maximum_temperature` and
`delta_temperature` require Unitful units. Example:

    Lactin2Model(rate_constant=0.1, maximum_temperature=40.0u"°C",
                 delta_temperature=2.0u"K", intercept=-1.0)
"""
struct Lactin2Model{C} <: AbstractPhenomenologicalModel
    rate_constant::Float64
    maximum_temperature::C
    delta_temperature::_KelvinQuantity
    intercept::Float64
end

Lactin2Model(; rate_constant, maximum_temperature, delta_temperature, intercept) =
    Lactin2Model(Float64(rate_constant), _temperature_param(maximum_temperature),
                 _z_value_param(delta_temperature), Float64(intercept))

function thermal_performance(m::Lactin2Model, T)
    Tc = _strict_C(T)
    Tmax = _C(m.maximum_temperature)
    exp(m.rate_constant * Tc) -
    exp(m.rate_constant * Tmax - (Tmax - Tc) / ustrip(m.delta_temperature)) +
    m.intercept
end
(m::Lactin2Model)(T) = thermal_performance(m, T)

lactin2(; kwargs...) = Lactin2Model(; kwargs...)
