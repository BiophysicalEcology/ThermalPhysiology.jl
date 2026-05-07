# ── Phenomenological TPC models ───────────────────────────────────────────────
#
# Parameters in °C; return absolute or relative performance.
# Internal symbols (T_opt, CTmax, rmax) follow literature convention.

"""
    thermal_performance(model, T) → Float64

Return the thermal performance (absolute or relative rate) at temperature `T`
for an [`AbstractPhenomenologicalModel`](@ref).

`T` may be a Unitful temperature quantity or a bare `Float64` (°C assumed).
"""
function thermal_performance end

# ── Universal TPC (Arnoldi, Jackson, Peralta-Maraver & Payne 2025, PNAS) ──────

"""
    UniversalTPCModel(; T_opt, E, maximum_performance)

Universal thermal performance curve (Arnoldi et al. 2025, PNAS).
Two-parameter model: `y(x) = exp(x) * (1 - x)` where `x = (T - T_opt) / E`.

- `T_opt`: optimal temperature (stored in K internally)
- `E`: thermal breadth (K); related to Arrhenius temperature by `E = T_ref² / T_A`
- `maximum_performance`: peak rate at T_opt

The UTPC thermal breadth connects to TDT z-value: `z = E * log(10)`.
"""
@kwdef struct UniversalTPCModel <: AbstractPhenomenologicalModel
    T_opt::Float64           = 303.15   # K
    E::Float64               = 10.0     # K
    maximum_performance::Float64 = 1.0
end

function thermal_performance(m::UniversalTPCModel, T)
    x = (_K(T) - m.T_opt) / m.E
    m.maximum_performance * exp(x) * (1 - x)
end
(m::UniversalTPCModel)(T) = thermal_performance(m, T)

"""
    utpc(; optimal_temperature, thermal_breadth, maximum_performance=1.0)

Named constructor for [`UniversalTPCModel`](@ref).
Accepts Unitful K/°C quantities or bare values (assumed °C for temperature, K for breadth).
"""
function utpc(;
    optimal_temperature  = 30.0,
    thermal_breadth      = 10.0u"K",
    maximum_performance  = 1.0,
)
    T_opt = _K(optimal_temperature)
    E     = thermal_breadth isa Unitful.Quantity ? ustrip(u"K", uconvert(u"K", thermal_breadth)) :
            Float64(thermal_breadth)
    UniversalTPCModel(T_opt=T_opt, E=E, maximum_performance=maximum_performance)
end

# ── Deutsch 2008 ──────────────────────────────────────────────────────────────

"""
    DeutschModel(; maximum_rate, optimal_temperature, critical_thermal_maximum, width_parameter)

Modified Deutsch (2008, PNAS) thermal performance curve.
Gaussian rise below T_opt; quadratic decline above T_opt to CTmax.

    T < T_opt: rate = maximum_rate * exp(-((T - T_opt) / (2*width_parameter))^2)
    T > T_opt: rate = maximum_rate * (1 - ((T - T_opt) / (T_opt - CTmax))^2)
"""
@kwdef struct DeutschModel <: AbstractPhenomenologicalModel
    maximum_rate::Float64             = 1.0
    optimal_temperature::Float64      = 25.0   # °C
    critical_thermal_maximum::Float64 = 40.0   # °C
    width_parameter::Float64          = 5.0    # °C (related to full curve width)
end

function thermal_performance(m::DeutschModel, T)
    Tc = _C(T)
    if Tc < m.optimal_temperature
        m.maximum_rate * exp(-((Tc - m.optimal_temperature) / (2 * m.width_parameter))^2)
    else
        m.maximum_rate * (1 - ((Tc - m.optimal_temperature) / (m.optimal_temperature - m.critical_thermal_maximum))^2)
    end
end
(m::DeutschModel)(T) = thermal_performance(m, T)

function deutsch(;
    maximum_rate             = 1.0,
    optimal_temperature      = 25.0,
    critical_thermal_maximum = 40.0,
    width_parameter          = 5.0,
)
    DeutschModel(
        maximum_rate=Float64(maximum_rate),
        optimal_temperature=_C(optimal_temperature),
        critical_thermal_maximum=_C(critical_thermal_maximum),
        width_parameter=Float64(width_parameter),
    )
end

# ── Briere 1 (1999) ───────────────────────────────────────────────────────────

"""
    Briere1Model(; rate_constant, minimum_temperature, maximum_temperature)

Briere 1 (1999) thermal performance curve.
`rate = rate_constant * T * (T - T_min) * sqrt(T_max - T)` for T_min < T < T_max.
"""
@kwdef struct Briere1Model <: AbstractPhenomenologicalModel
    rate_constant::Float64        = 0.01
    minimum_temperature::Float64  = 10.0   # °C
    maximum_temperature::Float64  = 40.0   # °C
end

function thermal_performance(m::Briere1Model, T)
    Tc = _C(T)
    Tc <= m.minimum_temperature || Tc >= m.maximum_temperature && return 0.0
    val = m.rate_constant * Tc * (Tc - m.minimum_temperature) * sqrt(m.maximum_temperature - Tc)
    max(0.0, val)
end
(m::Briere1Model)(T) = thermal_performance(m, T)

function briere(;
    rate_constant       = 0.01,
    minimum_temperature = 10.0,
    maximum_temperature = 40.0,
)
    Briere1Model(
        rate_constant=Float64(rate_constant),
        minimum_temperature=_C(minimum_temperature),
        maximum_temperature=_C(maximum_temperature),
    )
end

# ── Briere 2 (1999) ───────────────────────────────────────────────────────────

"""
    Briere2Model(; rate_constant, minimum_temperature, maximum_temperature, shape_parameter)

Briere 2 (1999) — generalised form with shape exponent.
`rate = rate_constant * T * (T - T_min) * (T_max - T)^(1/shape_parameter)`
"""
@kwdef struct Briere2Model <: AbstractPhenomenologicalModel
    rate_constant::Float64        = 0.01
    minimum_temperature::Float64  = 10.0
    maximum_temperature::Float64  = 40.0
    shape_parameter::Float64      = 2.0
end

function thermal_performance(m::Briere2Model, T)
    Tc = _C(T)
    Tc <= m.minimum_temperature || Tc >= m.maximum_temperature && return 0.0
    val = m.rate_constant * Tc * (Tc - m.minimum_temperature) *
          (m.maximum_temperature - Tc)^(1 / m.shape_parameter)
    max(0.0, val)
end
(m::Briere2Model)(T) = thermal_performance(m, T)

# ── Gaussian ──────────────────────────────────────────────────────────────────

"""
    GaussianModel(; maximum_rate, optimal_temperature, width_parameter)

Symmetric Gaussian TPC.
`rate = maximum_rate * exp(-0.5 * ((T - T_opt) / width_parameter)^2)`
"""
@kwdef struct GaussianModel <: AbstractPhenomenologicalModel
    maximum_rate::Float64        = 1.0
    optimal_temperature::Float64 = 25.0
    width_parameter::Float64     = 5.0
end

function thermal_performance(m::GaussianModel, T)
    Tc = _C(T)
    m.maximum_rate * exp(-0.5 * ((Tc - m.optimal_temperature) / m.width_parameter)^2)
end
(m::GaussianModel)(T) = thermal_performance(m, T)

function gaussian(;
    maximum_rate        = 1.0,
    optimal_temperature = 25.0,
    width_parameter     = 5.0,
)
    GaussianModel(
        maximum_rate=Float64(maximum_rate),
        optimal_temperature=_C(optimal_temperature),
        width_parameter=Float64(width_parameter),
    )
end

# ── Thomas 2012 ───────────────────────────────────────────────────────────────

"""
    Thomas2012Model(; rate_constant, shape_parameter, optimal_temperature)

Thomas et al. (2012) thermal performance curve.
`rate = rate_constant * (T - T_opt + b) * (T - T_opt - b) * (-1)`
where b is a shape parameter. Returns 0 outside the performance range.
"""
@kwdef struct Thomas2012Model <: AbstractPhenomenologicalModel
    rate_constant::Float64       = 0.5
    shape_parameter::Float64     = 15.0   # °C; half-width of performance curve
    optimal_temperature::Float64 = 25.0
end

function thermal_performance(m::Thomas2012Model, T)
    Tc = _C(T)
    b  = m.shape_parameter
    T0 = m.optimal_temperature
    val = m.rate_constant * (Tc - T0 + b) * (b - (Tc - T0))
    max(0.0, val)
end
(m::Thomas2012Model)(T) = thermal_performance(m, T)

function thomas(;
    rate_constant       = 0.5,
    shape_parameter     = 15.0,
    optimal_temperature = 25.0,
)
    Thomas2012Model(
        rate_constant=Float64(rate_constant),
        shape_parameter=Float64(shape_parameter),
        optimal_temperature=_C(optimal_temperature),
    )
end

# ── Thomas 2017 ───────────────────────────────────────────────────────────────

"""
    Thomas2017Model(; maximum_rate, optimal_temperature, width_parameter, skewness)

Thomas et al. (2017) asymmetric TPC: skewed-Gaussian form.
"""
@kwdef struct Thomas2017Model <: AbstractPhenomenologicalModel
    maximum_rate::Float64        = 1.0
    optimal_temperature::Float64 = 25.0
    width_parameter::Float64     = 5.0
    skewness::Float64            = 0.0   # positive → right skew
end

function thermal_performance(m::Thomas2017Model, T)
    Tc = _C(T)
    d  = Tc - m.optimal_temperature
    σ  = m.width_parameter * (1 + m.skewness * sign(d))
    m.maximum_rate * exp(-0.5 * (d / σ)^2)
end
(m::Thomas2017Model)(T) = thermal_performance(m, T)

# ── Pawar 2018 ────────────────────────────────────────────────────────────────

"""
    PawarModel(; rate_at_reference, activation_energy, deactivation_energy,
                 peak_temperature, reference_temperature)

Pawar et al. (2018) metabolic TPC: Arrhenius rise with high-T deactivation,
parameterised in biologically meaningful terms.
"""
@kwdef struct PawarModel <: AbstractPhenomenologicalModel
    rate_at_reference::Float64    = 1.0
    activation_energy::Float64    = 0.65   # eV
    deactivation_energy::Float64  = 1.15   # eV
    peak_temperature::Float64     = 30.0   # °C
    reference_temperature::Float64 = 20.0  # °C
end

function thermal_performance(m::PawarModel, T)
    Tk   = _K(T)
    Tref = _K(m.reference_temperature)
    Tpk  = _K(m.peak_temperature)
    k    = 8.617333e-5  # eV/K
    Ea   = m.activation_energy
    Ed   = m.deactivation_energy
    boltzmann  = exp(Ea / k * (1 / Tref - 1 / Tk))
    correction = 1 / (1 + exp((Ed / k) * (1 / Tpk - 1 / Tk)))
    m.rate_at_reference * boltzmann * correction
end
(m::PawarModel)(T) = thermal_performance(m, T)

function pawar(;
    rate_at_reference    = 1.0,
    activation_energy    = 0.65,
    deactivation_energy  = 1.15,
    peak_temperature     = 30.0,
    reference_temperature = 20.0,
)
    PawarModel(
        rate_at_reference=Float64(rate_at_reference),
        activation_energy=Float64(activation_energy),
        deactivation_energy=Float64(deactivation_energy),
        peak_temperature=_C(peak_temperature),
        reference_temperature=_C(reference_temperature),
    )
end

# ── Lactin 2 ──────────────────────────────────────────────────────────────────

"""
    Lactin2Model(; rate_constant, maximum_temperature, delta_temperature, intercept)

Lactin 2 (1995) thermal performance curve.
`rate = exp(rate_constant * T) - exp(rate_constant * T_max - (T_max - T) / delta_T) + intercept`
"""
@kwdef struct Lactin2Model <: AbstractPhenomenologicalModel
    rate_constant::Float64     = 0.1
    maximum_temperature::Float64 = 40.0
    delta_temperature::Float64  = 2.0
    intercept::Float64          = -1.0
end

function thermal_performance(m::Lactin2Model, T)
    Tc = _C(T)
    exp(m.rate_constant * Tc) -
    exp(m.rate_constant * m.maximum_temperature - (m.maximum_temperature - Tc) / m.delta_temperature) +
    m.intercept
end
(m::Lactin2Model)(T) = thermal_performance(m, T)

function lactin2(;
    rate_constant       = 0.1,
    maximum_temperature = 40.0,
    delta_temperature   = 2.0,
    intercept           = -1.0,
)
    Lactin2Model(
        rate_constant=Float64(rate_constant),
        maximum_temperature=_C(maximum_temperature),
        delta_temperature=Float64(delta_temperature),
        intercept=Float64(intercept),
    )
end
