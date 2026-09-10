using Roots: find_zero, Brent
using Statistics: mean

# ── TPC properties ─────────────────────────────────────────────────────────────

"""
    optimal_temperature(m) → Float64 (K)

Temperature at which `thermal_performance` is maximised.
Analytic for `UniversalTPCModel`; numerical via Roots.jl for others.
Returns `Inf` for monotone models (e.g. `ArrheniusModel`).
"""
optimal_temperature(m::UniversalTPCModel)  = m.T_opt
critical_thermal_minimum(::UniversalTPCModel; kwargs...) = -Inf * u"K"  # asymptotic approach to 0
optimal_temperature(m::DeutschModel)       = m.optimal_temperature
optimal_temperature(m::GaussianModel)      = m.optimal_temperature
optimal_temperature(m::Thomas2017Model)    = m.optimal_temperature
optimal_temperature(::ArrheniusModel)      = Inf * u"K"   # monotone — no finite peak

function optimal_temperature(m::AbstractTPCModel)
    # Numerical search over a reasonable °C range (bare Float64 domain; wrapped
    # in u"°C" at each call into thermal_performance, per the root-finder convention).
    T_search = (1.0, 80.0)
    try
        T_opt_C = find_zero(T -> _tpc_derivative(m, T), T_search, Brent())
        T_opt_C * u"°C"
    catch
        T_opt_C = _grid_maximum(m, T_search...)
        T_opt_C * u"°C"
    end
end

# Finite-difference derivative of thermal_performance
_tpc_derivative(m, T; h=0.01) =
    (thermal_performance(m, (T + h) * u"°C") - thermal_performance(m, (T - h) * u"°C")) / (2h)

function _grid_maximum(m, T_lo, T_hi; n=500)
    Ts = range(T_lo, T_hi, length=n)
    Ts[argmax(thermal_performance.(Ref(m), Ts .* u"°C"))]
end

"""
    maximum_rate(m) → Float64

Peak thermal performance (at `optimal_temperature`).
"""
maximum_rate(m::AbstractTPCModel) = thermal_performance(m, optimal_temperature(m))

"""
    critical_thermal_maximum(m; threshold=0.0) → Unitful temperature

Temperature above the optimum where `thermal_performance` falls to `threshold`.
"""
critical_thermal_maximum(m::DeutschModel)      = m.critical_thermal_maximum
critical_thermal_maximum(m::UniversalTPCModel) = m.T_opt + m.E   # root of y(x)=0 → x=1 → T=T_opt+E

function critical_thermal_maximum(m::AbstractTPCModel; threshold=0.0)
    T_opt_C = _C(optimal_temperature(m))   # search in °C — wrapped in u"°C" before each thermal_performance call
    try
        T_hi = T_opt_C + 50.0
        T_C  = find_zero(T -> thermal_performance(m, T * u"°C") - threshold, (T_opt_C, T_hi), Brent())
        T_C * u"°C"
    catch
        NaN * u"°C"
    end
end

"""
    critical_thermal_minimum(m; threshold=0.0) → Unitful temperature

Temperature below the optimum where `thermal_performance` falls to `threshold`.
"""
function critical_thermal_minimum(m::AbstractTPCModel; threshold=0.0)
    T_opt_C = _C(optimal_temperature(m))   # search in °C
    try
        T_lo = T_opt_C - 50.0
        T_C  = find_zero(T -> thermal_performance(m, T * u"°C") - threshold, (T_lo, T_opt_C), Brent())
        T_C * u"°C"
    catch
        NaN * u"°C"
    end
end

"""
    thermal_breadth(m; threshold=0.0) → Unitful temperature difference

Temperature range over which `thermal_performance` exceeds `threshold`.
= CTmax - CTmin.
"""
function thermal_breadth(m::AbstractTPCModel; threshold=0.0)
    ctmax = critical_thermal_maximum(m; threshold)
    ctmin = critical_thermal_minimum(m; threshold)
    uconvert(u"K", ctmax - ctmin)
end

"""
    q10(m, T; delta=10.0) → Float64

Temperature coefficient: ratio of performance at T+delta to T.
`delta` is in °C (or K — equivalent for differences).
"""
q10(m::AbstractTPCModel, T; delta=10.0) =
    thermal_performance(m, (_strict_K(T) + delta) * u"K") / thermal_performance(m, _strict_K(T) * u"K")

# ── TDT properties ─────────────────────────────────────────────────────────────

"""
    lethal_temperature(m, duration; T_search=(273.0, 373.0)) → Unitful temperature

Temperature at which `survival_time` equals `duration`.
Inverse of `survival_time`.
"""
function lethal_temperature(m::AbstractTDTModel, duration;
                            T_search=(20.0, 80.0))
    d = _min(_time_param(duration))
    find_zero(Tc -> _min_or_bare(survival_time(m, Tc * u"°C")) - d, T_search, Brent()) * u"°C"
end

"""
    median_lethal_temperature(m, duration) → Unitful temperature

Temperature at which 50% of individuals die after `duration`.
Alias for `lethal_temperature`.
"""
median_lethal_temperature(m::AbstractTDTModel, duration) = lethal_temperature(m, duration)

"""
    z_value(m) → Unitful K

Temperature increment for 10-fold change in survival time.
Direct field for `LogLinearTDTModel`/`ToleranceLandscape`; derived from
Arrhenius temperature for others.
"""
z_value(m::LogLinearTDTModel)  = m.z_value
z_value(m::ToleranceLandscape) = m.z_value
z_value(m::ArrheniusModel)     = (_K(m.T_ref)^2 / _K(m.T_A) * log(10)) * u"K"   # E = T_ref²/T_A; z = E*log(10)

function z_value(m::AbstractTDTModel)
    # Numeric: estimate slope of log10(t) ~ T
    T1, T2 = 35.0, 40.0
    t1 = _min_or_bare(survival_time(m, T1 * u"°C"))
    t2 = _min_or_bare(survival_time(m, T2 * u"°C"))
    (-1.0 / ((log10(t2) - log10(t1)) / (T2 - T1))) * u"K"
end

"""
    thermal_death_slope(m) → Float64

Slope b of the TDT curve in log10-time/°C units (= 1/z).
"""
thermal_death_slope(m::AbstractTDTModel) = log(10) / z_value(m)

# ── Constant Temperature Equivalent (CTE) ─────────────────────────────────────

"""
    constant_temperature_equivalent(m, T_series) → Unitful temperature

The single constant temperature that produces the same mean thermal response
as the time series `T_series`. Used in DEBtool parameter estimation when body
temperatures fluctuate (Kearney NicheMapR vignette, `getCTE`).

Algorithm:
1. Compute thermal response at each T → `responses = f.(m, T_series)`
2. Average → `mean_response`
3. Root-find `T_eq` such that `f(m, T_eq) = mean_response`

Analytic solution for `ArrheniusModel`; numerical for all others.

References:
- Kearney NicheMapR vignette, `getCTE` function
- Jensen's inequality: for convex f, CTE > arithmetic mean temperature
"""
function constant_temperature_equivalent(m::ArrheniusModel, T_series)
    mean_tc = mean(temperature_correction.(Ref(m), T_series))
    # Analytic: exp(T_A/T_ref - T_A/T_eq) = mean_tc → T_eq = T_A/(T_A/T_ref - log(mean_tc))
    m.T_A / (m.T_A / m.T_ref - log(mean_tc))
end

function constant_temperature_equivalent(m::AbstractArrheniusModel, T_series;
        T_bounds=nothing)
    Tk_series = _strict_K.(T_series)
    bounds = something(T_bounds, (minimum(Tk_series) - 5.0, maximum(Tk_series) + 5.0))
    mean_tc = mean(temperature_correction.(Ref(m), Tk_series .* u"K"))
    T_eq_K  = find_zero(T -> temperature_correction(m, T * u"K") - mean_tc, bounds, Brent())
    T_eq_K * u"K"
end

function constant_temperature_equivalent(m::AbstractPhenomenologicalModel, T_series;
        T_bounds=nothing)
    Tk_series = _strict_K.(T_series)
    bounds = something(T_bounds, (minimum(Tk_series) - 5.0, maximum(Tk_series) + 5.0))
    mean_perf = mean(thermal_performance.(Ref(m), Tk_series .* u"K"))
    T_eq_K    = find_zero(T -> thermal_performance(m, T * u"K") - mean_perf, bounds, Brent())
    T_eq_K * u"K"
end

function constant_temperature_equivalent(m::AbstractTDTModel, T_series;
        T_bounds=nothing)
    Tc_series = _strict_C.(T_series)
    bounds = something(T_bounds, (minimum(Tc_series) - 5.0, maximum(Tc_series) + 5.0))
    mean_damage_rate = mean(1.0 ./ _min_or_bare.(survival_time.(Ref(m), Tc_series .* u"°C")))
    T_eq_C = find_zero(Tc -> 1.0 / _min_or_bare(survival_time(m, Tc * u"°C")) - mean_damage_rate, bounds, Brent())
    T_eq_C * u"°C"
end

"""
    mean_correction_factor(m, T_series) → Float64

Mean Arrhenius temperature-correction factor over `T_series`.
"""
mean_correction_factor(m::AbstractArrheniusModel, T_series) =
    mean(temperature_correction.(Ref(m), T_series))

"""
    mean_thermal_performance(m, T_series) → Float64

Mean thermal performance over `T_series`.
"""
mean_thermal_performance(m::AbstractTPCModel, T_series) =
    mean(thermal_performance.(Ref(m), T_series))
