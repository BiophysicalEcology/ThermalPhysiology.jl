using LsqFit: curve_fit
using Statistics: mean, median

# ── Data container types ───────────────────────────────────────────────────────

"""
    IndividualKnockdownData(; temperatures, knockdown_times)

Individual-level knockdown data: one entry per organism.
Required for `fit_tolerance_landscape`.

- `temperatures`: assay temperature per individual (°C)
- `knockdown_times`: time to knockdown per individual (min)
"""
struct IndividualKnockdownData
    temperatures::Vector{Float64}
    knockdown_times::Vector{Float64}
end

function IndividualKnockdownData(; temperatures, knockdown_times)
    IndividualKnockdownData(Float64.(_C.(temperatures)), Float64.(knockdown_times))
end

"""
    StaticKnockdownData(; temperatures, knockdown_times)

Group-summary TDT data: one mean/median knockdown time per assay temperature.
Sufficient for `fit_thermal_death_time_curve`.

- `temperatures`: assay temperatures (°C)
- `knockdown_times`: mean/median time to failure per temperature (min)
"""
struct StaticKnockdownData
    temperatures::Vector{Float64}
    knockdown_times::Vector{Float64}
end

function StaticKnockdownData(; temperatures, knockdown_times)
    StaticKnockdownData(Float64.(_C.(temperatures)), Float64.(knockdown_times))
end

"""
    DynamicKnockdownData(; ramp_rates, dynamic_ctmax_values, start_temperature)

Dynamic CTmax measurements at various ramp rates.

- `ramp_rates`: heating rates (°C/min)
- `dynamic_ctmax_values`: knockdown temperature at each ramp rate (°C)
- `start_temperature`: temperature at which ramp begins (°C)
"""
struct DynamicKnockdownData
    ramp_rates::Vector{Float64}
    dynamic_ctmax_values::Vector{Float64}
    start_temperature::Float64
end

function DynamicKnockdownData(; ramp_rates, dynamic_ctmax_values, start_temperature)
    DynamicKnockdownData(Float64.(ramp_rates), Float64.(_C.(dynamic_ctmax_values)),
                         _C(start_temperature))
end

# ── TPC fitting ────────────────────────────────────────────────────────────────

"""
    fit_thermal_performance_curve(ModelType, temperatures, rates; initial_parameters=nothing)

Fit a TPC model to data using nonlinear least squares (LsqFit.jl).

Returns a fitted instance of `ModelType`.
`temperatures` and `rates` are plain vectors (°C assumed for temperatures).
`initial_parameters` is an ordered vector matching the struct field order; estimated
from data heuristics if `nothing`.
"""
function fit_thermal_performance_curve(::Type{M}, temperatures, rates;
        initial_parameters=nothing,
        weights=nothing) where {M <: AbstractTPCModel}
    Tc = Float64.(_C.(temperatures))
    y  = Float64.(rates)
    p0 = isnothing(initial_parameters) ? _tpc_initial_parameters(M, Tc, y) : Float64.(initial_parameters)
    model_fn(T_vec, p) = [_tpc_predict(M, t, p) for t in T_vec]
    fit = if isnothing(weights)
        curve_fit(model_fn, Tc, y, p0)
    else
        curve_fit(model_fn, Tc, y, Float64.(weights), p0)
    end
    _tpc_from_params(M, fit.param)
end

# Dispatch helpers — add for each model type
function _tpc_predict(::Type{UniversalTPCModel}, T, p)
    T_opt, E, P_max = p
    x = (T + 273.15 - T_opt) / E
    P_max * exp(x) * (1 - x)
end

function _tpc_initial_parameters(::Type{UniversalTPCModel}, Tc, y)
    idx = argmax(y)
    T_opt = Tc[idx] + 273.15
    E = (maximum(Tc) - minimum(Tc)) / 4.0
    P_max = maximum(y)
    [T_opt, E, P_max]
end

function _tpc_from_params(::Type{UniversalTPCModel}, p)
    UniversalTPCModel(T_opt=p[1], E=p[2], maximum_performance=p[3])
end

function _tpc_predict(::Type{GaussianModel}, T, p)
    P_max, T_opt, width = p
    P_max * exp(-0.5 * ((T - T_opt) / width)^2)
end

function _tpc_initial_parameters(::Type{GaussianModel}, Tc, y)
    idx = argmax(y)
    [maximum(y), Tc[idx], (maximum(Tc) - minimum(Tc)) / 4.0]
end

function _tpc_from_params(::Type{GaussianModel}, p)
    GaussianModel(maximum_rate=p[1], optimal_temperature=p[2], width_parameter=p[3])
end

function _tpc_predict(::Type{DeutschModel}, T, p)
    P_max, T_opt, CTmax, width = p
    if T < T_opt
        P_max * exp(-((T - T_opt) / (2 * width))^2)
    else
        P_max * (1 - ((T - T_opt) / (T_opt - CTmax))^2)
    end
end

function _tpc_initial_parameters(::Type{DeutschModel}, Tc, y)
    idx = argmax(y)
    P_max = maximum(y)
    T_opt = Tc[idx]
    CTmax = maximum(Tc) + 2.0
    width = (T_opt - minimum(Tc)) / 2.0
    [P_max, T_opt, CTmax, width]
end

function _tpc_from_params(::Type{DeutschModel}, p)
    DeutschModel(maximum_rate=p[1], optimal_temperature=p[2],
                 critical_thermal_maximum=p[3], width_parameter=p[4])
end

# ── OLS helper (used for Schoolfield initial estimates and TDT fitting) ─────────

function _ols_intercept_slope(x::AbstractVector, y::AbstractVector)
    xbar = mean(x);  ybar = mean(y)
    b = sum((x .- xbar) .* (y .- ybar)) / sum((x .- xbar).^2)
    (a = ybar - b * xbar, b = b)
end

# ── Schoolfield graphical initial parameter estimation ──────────────────────────
#
# Algorithm (Schoolfield, Sharpe & Magnuson 1981, Fig. 2 / R vignette):
#  Split sorted data into cold / middle / hot fractions.
#  Fit OLS on ln(rate) ~ 1/T for each region.
#  T_A   from middle slope:   b_mid  = -T_A
#  T_AL  from cold slope:     b_cold = -(T_A + T_AL)  → T_AL = -b_cold - T_A
#  T_AH  from hot slope:      b_hot  = T_AH - T_A     → T_AH = b_hot + T_A
#  T_L / T_H: intersection of cold / hot line with the half-Arrhenius line
#    (middle line shifted down by ln 2, i.e. the temperature where the
#     inactivation term = 1 and the rate is halved).

function _schoolfield_initial_params(temps_C::AbstractVector, rates::AbstractVector;
        cold_frac::Real=0.30, hot_frac::Real=0.30, T_ref_K::Real=298.15)
    n   = length(temps_C)
    idx = sortperm(temps_C)
    Ts  = temps_C[idx]
    rs  = rates[idx]

    n_cold = max(2, round(Int, n * cold_frac))
    n_hot  = max(2, round(Int, n * hot_frac))
    n_mid  = n - n_cold - n_hot
    if n_mid < 2
        n_cold = max(2, div(n, 3))
        n_hot  = max(2, div(n, 3))
        n_mid  = n - n_cold - n_hot
    end

    inv_T_K(Tc) = 1.0 ./ (Tc .+ 273.15)

    cold = 1:n_cold
    mid  = (n_cold+1):(n_cold+n_mid)
    hot  = (n_cold+n_mid+1):n

    (a_cold, b_cold) = let r = _ols_intercept_slope(inv_T_K(Ts[cold]), log.(rs[cold]));  (r.a, r.b) end
    (a_mid,  b_mid)  = let r = _ols_intercept_slope(inv_T_K(Ts[mid]),  log.(rs[mid]));   (r.a, r.b) end
    (a_hot,  b_hot)  = let r = _ols_intercept_slope(inv_T_K(Ts[hot]),  log.(rs[hot]));   (r.a, r.b) end

    T_A  = max(100.0, -b_mid)
    T_AL = max(1000.0, -b_cold - T_A)    # = -b_cold + b_mid
    T_AH = max(1000.0,  b_hot  + T_A)    # = b_hot  - b_mid

    # Intersection of cold OLS line with half-Arrhenius line: 1/T_L = (a_mid - ln2 - a_cold)/(b_cold - b_mid)
    dL = b_cold - b_mid
    T_L = if abs(dL) > 1.0
        clamp((b_cold - b_mid) / (a_mid - log(2) - a_cold), 240.0, T_ref_K - 1.0)
    else
        Ts[n_cold] + 273.15
    end

    # Intersection of hot OLS line with half-Arrhenius line: 1/T_H = (a_mid - ln2 - a_hot)/(b_hot - b_mid)
    dH = b_hot - b_mid
    T_H = if abs(dH) > 1.0
        clamp((b_hot - b_mid) / (a_mid - log(2) - a_hot), T_ref_K + 1.0, 360.0)
    else
        Ts[n_cold + n_mid] + 273.15
    end

    rate_ref = exp(a_mid + b_mid / T_ref_K)

    (T_A=T_A, T_L=T_L, T_AL=T_AL, T_H=T_H, T_AH=T_AH, rate_at_reference=rate_ref)
end

# ── Schoolfield prediction helpers (used in specialized fit overload) ───────────

function _ss_predict(::Type{SharpSchoolFullModel}, Tk, p, T_ref_K)
    T_A, T_L, T_AL, T_H, T_AH, rate_ref = p
    rate_ref * (Tk / T_ref_K) * exp(T_A / T_ref_K - T_A / Tk) /
        (1 + exp(T_AL / Tk - T_AL / T_L) + exp(T_AH / T_H - T_AH / Tk))
end

function _ss_predict(::Type{SharpSchoolDEBModel}, Tk, p, T_ref_K)
    T_A, T_L, T_AL, T_H, T_AH, rate_ref = p
    boltz    = exp(T_A / T_ref_K - T_A / Tk)
    low_ref  = exp(T_AL / T_ref_K - T_AL / T_L)
    high_ref = exp(T_AH / T_H - T_AH / T_ref_K)
    low_T    = exp(T_AL / Tk   - T_AL / T_L)
    high_T   = exp(T_AH / T_H - T_AH / Tk)
    rate_ref * boltz * (1 + low_ref + high_ref) / (1 + low_T + high_T)
end

"""
    fit_thermal_performance_curve(SharpSchoolFullModel, temperatures, rates; ...)
    fit_thermal_performance_curve(SharpSchoolDEBModel, temperatures, rates; ...)

Fit a Sharpe-Schoolfield model using NLS (LsqFit.jl) with graphical initial parameter
estimates from the Arrhenius plot (Schoolfield, Sharpe & Magnuson 1981, Figure 2).

Keyword arguments:
- `T_ref`: reference temperature (Unitful quantity or bare Kelvin); default `298.15u"K"` (25 °C).
- `initial_parameters`: manual override `[T_A, T_L, T_AL, T_H, T_AH, rate_at_reference]`;
  if `nothing`, estimated from the graphical method.
- `log_transform`: fit in log-space (minimise sum of squared log-residuals); default `true`.
  This matches the Schoolfield et al. (1981) Marquardt fit and is appropriate when rates
  span orders of magnitude. Set `false` for absolute residuals (use with `weights`).
- `weights`: weight vector for LsqFit (only used when `log_transform=false`).
"""
function fit_thermal_performance_curve(::Type{M}, temperatures, rates;
        initial_parameters = nothing,
        T_ref              = 298.15u"K",
        log_transform      = true,
        weights            = nothing) where {M <: Union{SharpSchoolFullModel, SharpSchoolDEBModel}}
    Tc      = Float64.(_C.(temperatures))
    y       = Float64.(rates)
    T_ref_K = _K(T_ref)

    p0 = if isnothing(initial_parameters)
        init = _schoolfield_initial_params(Tc, y; T_ref_K=T_ref_K)
        [init.T_A, init.T_L, init.T_AL, init.T_H, init.T_AH, init.rate_at_reference]
    else
        Float64.(initial_parameters)
    end

    if log_transform
        log_model_fn(T_vec, p) = [log(_ss_predict(M, t + 273.15, p, T_ref_K)) for t in T_vec]
        fit = curve_fit(log_model_fn, Tc, log.(y), p0)
    else
        model_fn(T_vec, p) = [_ss_predict(M, t + 273.15, p, T_ref_K) for t in T_vec]
        fit = if isnothing(weights)
            curve_fit(model_fn, Tc, y, p0)
        else
            curve_fit(model_fn, Tc, y, Float64.(weights), p0)
        end
    end

    T_A, T_L, T_AL, T_H, T_AH, rate_ref = fit.param
    M(T_A=T_A, T_ref=T_ref_K, T_L=T_L, T_AL=T_AL, T_H=T_H, T_AH=T_AH,
      rate_at_reference=rate_ref)
end

# ── TDT fitting — static data ──────────────────────────────────────────────────

"""
    fit_thermal_death_time_curve(data::StaticKnockdownData; reference_duration=60.0,
                                  incipient_temperature=nothing)

Fit a `LogLinearTDTModel` to group-summary knockdown data via linear regression
on `log10(t) ~ T` (TDT_from_Static.R, lines 68–71).

Returns a `LogLinearTDTModel`. If `incipient_temperature` is `nothing`, defaults
to `minimum(temperatures) - 5.0`.
"""
function fit_thermal_death_time_curve(data::StaticKnockdownData;
        reference_duration::Real = 60.0,
        incipient_temperature    = nothing)
    Tc    = data.temperatures
    log_t = log10.(data.knockdown_times)
    (a, b) = let r = _ols_intercept_slope(Tc, log_t); (r.a, r.b) end
    z_val = -1.0 / b
    # reference_ctmax: T at which log10(t) = log10(reference_duration) → T = (log10(ref) - a)/b
    ref_ctmax = (log10(reference_duration) - a) / b
    T_inc = isnothing(incipient_temperature) ? minimum(Tc) - 5.0 : _C(incipient_temperature)
    LogLinearTDTModel(z_value=z_val, reference_ctmax=ref_ctmax,
                      reference_duration=Float64(reference_duration),
                      incipient_temperature=T_inc)
end

# ── TDT fitting — dynamic data ─────────────────────────────────────────────────

"""
    fit_thermal_death_time_curve(data::DynamicKnockdownData; reference_duration=60.0,
                                  incipient_temperature=nothing, initial_z=2.5)

Fit a `LogLinearTDTModel` from dynamic CTmax measurements (Jørgensen 2021 Eq. 7a).

- ≥3 ramp rates: NLS via LsqFit.jl
- 2 ramp rates: Roots.jl (multiroot)
- 1 ramp rate: algebraic inverse of Eq. 7a (requires reference_ctmax guess = mean dCTmax - 1°C)
"""
function fit_thermal_death_time_curve(data::DynamicKnockdownData;
        reference_duration::Real = 60.0,
        incipient_temperature    = nothing,
        initial_z::Real          = 2.5)
    rr   = data.ramp_rates
    dctm = data.dynamic_ctmax_values
    T0   = data.start_temperature
    T_inc = isnothing(incipient_temperature) ? T0 : _C(incipient_temperature)

    # Predict dCTmax given parameters (Jørgensen 2021 Eq. 7a)
    function predict_dctmax(z_val, ref_ctmax)
        m_tmp = LogLinearTDTModel(z_value=z_val, reference_ctmax=ref_ctmax,
                                  reference_duration=Float64(reference_duration),
                                  incipient_temperature=Float64(T_inc))
        [dynamic_ctmax(m_tmp, r; start_temperature=T0) for r in rr]
    end

    if length(rr) >= 3
        p0 = [initial_z, mean(dctm) - 1.0]
        fit = curve_fit((_, p) -> predict_dctmax(p[1], p[2]), ones(length(rr)), dctm, p0)
        z_fit, ref_fit = fit.param
    elseif length(rr) == 2
        # Two-equation system: scan z, solve for ref_ctmax at each ramp rate
        function two_rate_residual(p)
            z_val, ref_ctmax = p[1], p[2]
            pred = predict_dctmax(z_val, ref_ctmax)
            [pred[1] - dctm[1], pred[2] - dctm[2]]
        end
        # Numerical solve: scan for z, solve for ref_ctmax at each
        z_fit = initial_z
        ref_fit = mean(dctm) - 1.0
        for z_try in range(0.5, 10.0, length=50)
            try
                ref_try = find_zero(
                    rc -> predict_dctmax(z_try, rc)[1] - dctm[1],
                    (minimum(dctm) - 5.0, maximum(dctm) + 5.0), Brent())
                if abs(predict_dctmax(z_try, ref_try)[2] - dctm[2]) <
                   abs(predict_dctmax(z_fit, ref_fit)[2] - dctm[2])
                    z_fit, ref_fit = z_try, ref_try
                end
            catch; end
        end
    else
        # Single ramp rate: algebraic (fix z to initial_z, solve for ref_ctmax)
        z_fit = initial_z
        ref_fit = find_zero(
            rc -> predict_dctmax(z_fit, rc)[1] - dctm[1],
            (minimum(dctm) - 10.0, maximum(dctm) + 10.0), Brent())
    end

    LogLinearTDTModel(z_value=z_fit, reference_ctmax=ref_fit,
                      reference_duration=Float64(reference_duration),
                      incipient_temperature=Float64(T_inc))
end

# ── ToleranceLandscape fitting ─────────────────────────────────────────────────

"""
    fit_tolerance_landscape(data::IndividualKnockdownData; n_bins=1000)

Build a `ToleranceLandscape` from individual-level knockdown data (Rezende et al. 2014,
Funct Ecol; corresponds to `tolerance.landscape()` in the R reference code).

Algorithm:
1. Linear regression: `log10(time) ~ temperature` to estimate `z_value` (slope = -1/z)
   and `temperature_maximum` (T at log10(t) = 0, i.e. t = 1 min).
2. Compute mean assay temperature `T_mean = mean(temperatures)`.
3. Z-shift each individual knockdown time to `T_mean`:
   `shifted_time = t * 10^((T - T_mean) / z_value)`
4. Sort shifted times; compute empirical survival fractions S(τ).
5. Store as `n_bins`-row `survival_curve` matrix via linear interpolation for consistency.

Returns `ToleranceLandscape`.
"""
function fit_tolerance_landscape(data::IndividualKnockdownData; n_bins::Int=1000)
    Tc    = data.temperatures
    t     = data.knockdown_times
    log_t = log10.(t)

    # Step 1: OLS log10(t) ~ T → z_value, T_max
    (a, b) = let r = _ols_intercept_slope(Tc, log_t); (r.a, r.b) end
    z_val = -1.0 / b
    T_max = -a / b         # T at log10(t)=0, i.e. t=1 min

    # Step 2: mean assay temperature
    T_mean = mean(Tc)

    # Step 3: z-shift all knockdown times to T_mean
    shifted = t .* 10.0 .^ ((Tc .- T_mean) ./ z_val)

    # Step 4: empirical survival function
    sorted_t = sort(shifted)
    n_ind    = length(sorted_t)
    surv_raw = [(n_ind - i) / n_ind for i in 0:(n_ind-1)]  # S(t) at each sorted time

    # Step 5: interpolate to n_bins rows for a regular grid
    t_lo, t_hi = sorted_t[1], sorted_t[end]
    t_grid  = range(t_lo, t_hi, length=n_bins)
    s_grid  = [_interp_survival(sorted_t, surv_raw, ti) for ti in t_grid]
    curve   = hcat(collect(t_grid), s_grid)

    ToleranceLandscape(z_val, T_max, T_mean, curve)
end
