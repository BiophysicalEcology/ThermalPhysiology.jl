using Statistics: mean, median

# ── Thermal death time models ──────────────────────────────────────────────────
#
# Two tiers:
#   LogLinearTDTModel  — parametric, median knockdown times, simple injury accumulation
#                        (Jørgensen et al. 2021, Sci Reports)
#   ToleranceLandscape — semi-parametric, full S(τ) distribution, iterative curve-shifting
#                        (Rezende et al. 2020, Science; Rezende et al. 2014, Funct Ecol)

"""
    survival_time(model, T) → Float64

Return the median time to knockdown (minutes) at constant temperature `T`
for an [`AbstractTDTModel`](@ref).

`T` may be a Unitful temperature quantity or a bare `Float64` (°C assumed).
"""
function survival_time end

# ── LogLinearTDTModel ─────────────────────────────────────────────────────────

"""
    LogLinearTDTModel(; z_value, reference_ctmax, reference_duration, incipient_temperature)

Log-linear thermal death time model (Jørgensen et al. 2021, Sci Reports).

Survival time at constant temperature T:
    t(T) = reference_duration * 10^((reference_ctmax - T) / z_value)

Parameters:
- `z_value`: °C for 10× change in knockdown time (= -1/slope of log10(t) ~ T)
- `reference_ctmax`: sCTmax at `reference_duration` (°C)
- `reference_duration`: exposure duration defining `reference_ctmax` (min)
- `incipient_temperature`: Tc* below which thermal injury is negligible (°C)

Rezende's T_max (mean τ = 1 min) relates by: T_max = reference_ctmax + z * log10(reference_duration).
"""
@kwdef struct LogLinearTDTModel <: AbstractTDTModel
    z_value::Float64              = 4.0
    reference_ctmax::Float64      = 39.0   # sCTmax at reference_duration (°C)
    reference_duration::Float64   = 60.0   # minutes
    incipient_temperature::Float64 = 30.0  # °C
end

survival_time(m::LogLinearTDTModel, T) =
    m.reference_duration * 10^((m.reference_ctmax - _C(T)) / m.z_value)
(m::LogLinearTDTModel)(T) = survival_time(m, T)

"""
    temperature_maximum(m::LogLinearTDTModel) → Float64

Temperature (°C) at which the mean knockdown time equals 1 minute (Rezende parameterisation).

    T_max = reference_ctmax + z_value × log10(reference_duration)
"""
temperature_maximum(m::LogLinearTDTModel) =
    m.reference_ctmax + m.z_value * log10(m.reference_duration)

function log_linear_tdt(;
    z_value              = 4.0,
    reference_ctmax      = 39.0,
    reference_duration   = 60.0,
    incipient_temperature = 30.0,
)
    LogLinearTDTModel(
        z_value=Float64(z_value),
        reference_ctmax=_C(reference_ctmax),
        reference_duration=Float64(reference_duration),
        incipient_temperature=_C(incipient_temperature),
    )
end

# ── Dual-use Arrhenius: survival_time via temperature_correction ───────────────

survival_time(m::ArrheniusModel, T; critical_damage=1.0) =
    critical_damage / temperature_correction(m, T)

# ── Repair models ───────────────────────────────────────────────────────────────
# Damage accrual (LogLinearTDTModel's own incipient_temperature/survival_time)
# and recovery are independent, composable processes -- same idea as an
# organism's diapause clock breaking while development proceeds in parallel.
# `step_injury` is the one shared per-step update every injury function below
# is built from; `repair_rate`/`resets_injury` are what a repair model plugs
# into it.

"""
    repair_rate(model, temperature, dt_minutes) → Float64

Continuous reduction in injury over `dt_minutes` at `temperature` (same units
as injury itself, i.e. a fraction of the lethal dose). `0.0` for models with
no continuous repair (an instantaneous [`resets_injury`](@ref) reset, or none
at all).
"""
function repair_rate end

"""
    resets_injury(model, temperature) → Bool

Whether injury should be wiped to 0 outright this step (rather than reduced
gradually via [`repair_rate`](@ref)) at `temperature`.
"""
function resets_injury end

"""
    NoRepair()

No recovery from thermal injury -- accumulation is strictly monotone
(Jørgensen et al. 2021's own convention). Used by [`accumulated_injury`](@ref).
"""
struct NoRepair <: AbstractRepairModel end
repair_rate(::NoRepair, temperature, dt_minutes) = 0.0
resets_injury(::NoRepair, temperature) = false

"""
    FullRepairBelowThreshold(threshold)
    FullRepairBelowThreshold(; threshold)

Injury resets fully to 0 whenever temperature drops below `threshold` (°C),
unless the lethal dose (injury = 1.0) has already been reached. The simplest
possible repair model -- full, instant recovery below some safe temperature
(e.g. overnight cooling), rather than [`NoRepair`](@ref)'s never-recovers
accumulation or a gradual [`repair_rate`](@ref). Used by
[`resettable_injury`](@ref) (with `threshold = m.incipient_temperature`).
"""
struct FullRepairBelowThreshold <: AbstractRepairModel
    threshold::Float64
end
resets_injury(m::FullRepairBelowThreshold, temperature) = _C(temperature) < m.threshold
repair_rate(::FullRepairBelowThreshold, temperature, dt_minutes) = 0.0

"""
    full_repair_below_threshold(threshold)

Named constructor for [`FullRepairBelowThreshold`](@ref); `threshold` may be
Unitful or a bare Float64 (assumed °C).
"""
full_repair_below_threshold(threshold) = FullRepairBelowThreshold(_C(threshold))

# ── Injury accumulation (LogLinearTDTModel) ────────────────────────────────────

"""
    step_injury(m::LogLinearTDTModel, repair::AbstractRepairModel, injury, temperature, dt_minutes) → Float64

One time-step update of accumulated thermal injury: adds `dt_minutes /
survival_time(m, temperature)` when `temperature >= m.incipient_temperature`
(else 0), clamped to a lethal dose of 1.0, then applies `repair` (an
instantaneous reset via [`resets_injury`](@ref), or a continuous reduction
via [`repair_rate`](@ref)) -- unless the lethal dose was just reached, which
no repair model can undo.

This is the one shared step [`accumulated_injury`](@ref) and
[`resettable_injury`](@ref) are built from, and is also what a driving
simulation (e.g. an hourly development-model loop) should call directly when
injury needs to persist as part of that simulation's own state, rather than
being recomputed from a whole trajectory after the fact.
"""
function step_injury(m::LogLinearTDTModel, repair::AbstractRepairModel, injury, temperature, dt_minutes)
    increment = _C(temperature) >= m.incipient_temperature ? dt_minutes / survival_time(m, temperature) : 0.0
    new_injury = min(1.0, injury + increment)
    new_injury >= 1.0 && return new_injury
    resets_injury(repair, temperature) ? 0.0 :
        clamp(new_injury - repair_rate(repair, temperature, dt_minutes), 0.0, 1.0)
end

function _run_injury(m::LogLinearTDTModel, repair::AbstractRepairModel, T_series, dt_minutes)
    injury = Vector{Float64}(undef, length(T_series))
    cumulative = 0.0
    for (i, T) in enumerate(T_series)
        cumulative = step_injury(m, repair, cumulative, T, dt_minutes)
        injury[i] = cumulative
    end
    injury
end

"""
    accumulated_injury(m, T_series, dt_minutes)

Cumulative thermal injury under a fluctuating temperature series (Jørgensen 2021 Eq. 3),
via [`step_injury`](@ref) with [`NoRepair`](@ref) (monotone, never recovers).
Injury per step = `dt / survival_time(T)` when T ≥ incipient_temperature, else 0.
Injury = 1.0 corresponds to lethal dose.
Returns a cumulative vector of length `length(T_series)`.
"""
accumulated_injury(m::LogLinearTDTModel, T_series, dt_minutes) = _run_injury(m, NoRepair(), T_series, dt_minutes)

"""
    time_to_failure(m, T_series, dt_minutes; critical_injury=1.0, resettable=false)

Time (minutes) at which accumulated injury reaches `critical_injury`.
Returns `Inf` if the organism survives the entire series.

`resettable=false` (default) uses [`accumulated_injury`](@ref) (monotone,
never recovers). `resettable=true` uses [`resettable_injury`](@ref) instead
(injury resets to 0 whenever `T < incipient_temperature`, unless already at
`critical_injury`).
"""
function time_to_failure(m::LogLinearTDTModel, T_series, dt_minutes;
                          critical_injury=1.0, resettable::Bool=false)
    injuries = resettable ? resettable_injury(m, T_series, dt_minutes) :
                            accumulated_injury(m, T_series, dt_minutes)
    idx = findfirst(>=(critical_injury), injuries)
    isnothing(idx) ? Inf : idx * dt_minutes
end

"""
    resettable_injury(m, T_series, dt_minutes)

Cumulative thermal injury under a fluctuating temperature series, like
[`accumulated_injury`](@ref), but injury resets fully to 0 whenever
`T < m.incipient_temperature` and the organism hasn't yet reached the lethal
dose (injury < 1.0) -- i.e. a brief drop back below the incipient temperature
is treated as complete repair, rather than [`accumulated_injury`](@ref)'s
strictly monotone (never-recovers) accumulation.

This is an individual-level, "hard reset" alternative to
[`dynamic_survival`](@ref)'s gradual, `recovery_model`-driven recovery --
appropriate when sub-incipient exposure is assumed to fully clear
accumulated heat injury (e.g. overnight cooling), rather than partially.

Injury = 1.0 corresponds to lethal dose. Returns a vector of length
`length(T_series)`.
"""
resettable_injury(m::LogLinearTDTModel, T_series, dt_minutes) =
    _run_injury(m, FullRepairBelowThreshold(m.incipient_temperature), T_series, dt_minutes)

# ── Static ↔ dynamic CTmax conversions (Jørgensen 2021 Eqs. 7a/7b) ───────────

"""
    dynamic_ctmax(m, ramp_rate_per_min; start_temperature=m.incipient_temperature)

Predict the dynamic CTmax (knockdown temperature in a ramping assay) from TDT parameters.
Reference: Jørgensen et al. 2021 Eq. 7a; `TDT_from_Static.R`.
"""
function dynamic_ctmax(m::LogLinearTDTModel, ramp_rate_per_min;
                       start_temperature=m.incipient_temperature)
    k  = log(10) / m.z_value
    T0 = _C(start_temperature)
    Tc = _C(m.incipient_temperature)
    T0 + (1/k) * log(k * ramp_rate_per_min * m.reference_duration *
                     exp(k * (m.reference_ctmax - T0)) + exp(k * (Tc - T0)))
end

"""
    static_ctmax_from_dynamic(m, dctmax, ramp_rate_per_min; start_temperature=m.incipient_temperature)

Recover static sCTmax from a dynamic CTmax measurement and ramp rate.
Reference: Jørgensen et al. 2021 Eq. 7b.
"""
function static_ctmax_from_dynamic(m::LogLinearTDTModel, dctmax, ramp_rate_per_min;
                                   start_temperature=m.incipient_temperature)
    k  = log(10) / m.z_value
    T0 = _C(start_temperature)
    Tc = _C(m.incipient_temperature)
    T0 + (1/k) * log((1 / (k * ramp_rate_per_min * m.reference_duration)) *
                     (exp(k * (_C(dctmax) - T0)) - exp(k * (Tc - T0))))
end

"""
    ctmax_at_duration(m, duration_minutes)

Static sCTmax corresponding to a given exposure duration.
"""
ctmax_at_duration(m::LogLinearTDTModel, duration_minutes) =
    m.reference_ctmax + m.z_value * log10(m.reference_duration / duration_minutes)

# ── ToleranceLandscape ────────────────────────────────────────────────────────

"""
    ToleranceLandscape(; z_value, temperature_maximum, mean_assay_temperature, survival_curve)

Semi-parametric thermal tolerance landscape (Rezende et al. 2014, Funct Ecol;
Rezende et al. 2020, Science 369, 1242–1245).

Built from **individual-level** knockdown data via [`fit_tolerance_landscape`](@ref).
Stores the full empirical survival probability curve S(τ) collapsed to the mean
assay temperature using z-shifting.

**Individual vs. population framing:**
`dynamic_survival` tracks the survival probability of *one individual* moving through
S(τ) as heat exposure accumulates. The state variable can move in both directions:
down under heat stress, up during recovery (via `recovery_model`). This contrasts
with `accumulated_injury` (Jørgensen), which is monotone and population-level.

Fields:
- `z_value`: °C for 10× change in knockdown time
- `temperature_maximum`: T_max, temperature at which mean τ = 1 min (°C)
- `mean_assay_temperature`: T_mean, mean assay temperature for S(τ) curve (°C)
- `survival_curve`: n×2 matrix [time_minutes, survival_fraction (0–1)]
"""
struct ToleranceLandscape <: AbstractTDTModel
    z_value::Float64
    temperature_maximum::Float64        # T_max (°C); Rezende's CTmax definition
    mean_assay_temperature::Float64     # T_mean (°C)
    survival_curve::Matrix{Float64}     # [time_minutes, survival_fraction]
end

function tolerance_landscape(;
    z_value,
    temperature_maximum,
    mean_assay_temperature,
    survival_curve,
)
    ToleranceLandscape(
        Float64(z_value),
        _C(temperature_maximum),
        _C(mean_assay_temperature),
        survival_curve,
    )
end

"""
    survival_time(tl::ToleranceLandscape, T)

Median knockdown time at temperature T, obtained by z-shifting the stored S(τ) curve.
"""
function survival_time(tl::ToleranceLandscape, T)
    Tc    = _C(T)
    shift = 10^((Tc - tl.mean_assay_temperature) / tl.z_value)
    times = tl.survival_curve[:, 1]
    surv  = tl.survival_curve[:, 2]
    # Median = time at survival = 0.5 in the shifted curve
    shifted_times = times .* shift
    # Interpolate: find time where survival crosses 0.5
    idx = searchsortedfirst(surv, 0.5, rev=true)
    idx = clamp(idx, 1, length(surv) - 1)
    # Linear interpolation between adjacent points
    s1, s2 = surv[idx], surv[idx+1]
    t1, t2 = shifted_times[idx], shifted_times[idx+1]
    s1 ≈ s2 ? t1 : t1 + (0.5 - s1) * (t2 - t1) / (s2 - s1)
end

(tl::ToleranceLandscape)(T) = survival_time(tl, T)

# ── Dynamic survival prediction (Rezende 2020 iterative curve-shifting) ────────

"""
    dynamic_survival(tl, T_series; dt_minutes=1.0, recovery_model=nothing)

Predict individual survival probability under a fluctuating temperature series
using the Rezende (2020) iterative curve-shifting algorithm.

At each time step:
1. Compute `shift = 10^((T_mean - T_i) / z)` — scales time axis for current temperature
2. Find current time-equivalent position in the shifted S(τ) and advance by `dt_minutes`
3. Read new survival fraction via linear interpolation
4. If `recovery_model` is provided, add `thermal_performance(recovery_model, T_i) * dt_minutes`
   and clamp to 1.0 (matches `ArrFunc5` recovery in `dynamic.landscape1.R`)

Returns a vector of survival fractions (0–1) at each time step.

Reference: `dynamic.landscape()` in `Thermal landscape functions.R` (lines 76–89);
`dynamic.landscape1()` in `dynamic.lansdcape1.R`.
"""
function dynamic_survival(tl::ToleranceLandscape, T_series;
                          dt_minutes::Real = 1.0,
                          recovery_model = nothing)
    surv_ref  = tl.survival_curve[:, 2]   # survival fractions at T_mean (decreasing)
    time_ref  = tl.survival_curve[:, 1]   # times at T_mean

    alive_vec = Vector{Float64}(undef, length(T_series))
    time_rel  = 0.0   # current effective time position in reference S(τ)
    alive     = 1.0   # current survival fraction

    for (i, T) in enumerate(T_series)
        alive <= 0.0 && (alive_vec[i:end] .= 0.0; break)

        Tc    = _C(T)
        shift = 10^((tl.mean_assay_temperature - Tc) / tl.z_value)
        shifted_times = time_ref .* shift

        # New effective time after dt_minutes at current temperature
        new_time_rel = time_rel + dt_minutes

        # Interpolate: find survival at new_time_rel in shifted curve
        new_alive = _interp_survival(shifted_times, surv_ref, new_time_rel)

        # Optional recovery (per-step addition, clamped to 1.0)
        if recovery_model !== nothing
            rec = if recovery_model isa AbstractArrheniusModel
                temperature_correction(recovery_model, T) * dt_minutes
            else
                thermal_performance(recovery_model, T) * dt_minutes
            end
            new_alive = min(1.0, new_alive + rec)
        end

        alive_vec[i] = new_alive
        alive = new_alive

        # Update time_rel: find effective position for next step
        if i < length(T_series)
            next_T = T_series[i+1]
            next_shift = 10^((tl.mean_assay_temperature - _C(next_T)) / tl.z_value)
            next_shifted = time_ref .* next_shift
            time_rel = _interp_time(surv_ref, next_shifted, alive)
        end
    end

    alive_vec
end

# Linear interpolation helpers
function _interp_survival(times, surv, t)
    # times ascending, surv decreasing; find survival at time t
    t <= times[1]   && return surv[1]
    t >= times[end] && return 0.0
    idx = searchsortedfirst(times, t) - 1
    idx = clamp(idx, 1, length(times) - 1)
    t1, t2 = times[idx], times[idx+1]
    s1, s2 = surv[idx], surv[idx+1]
    t1 ≈ t2 ? s1 : s1 + (t - t1) * (s2 - s1) / (t2 - t1)
end

function _interp_time(surv, times, s)
    # surv decreasing, times ascending; find time at survival s
    s >= surv[1]   && return times[1]
    s <= surv[end] && return times[end]
    idx = searchsortedfirst(surv, s, rev=true)
    idx = clamp(idx, 1, length(surv) - 1)
    s1, s2 = surv[idx], surv[idx+1]
    t1, t2 = times[idx], times[idx+1]
    s1 ≈ s2 ? t1 : t1 + (s - s1) * (t2 - t1) / (s2 - s1)
end

# ── Multi-day mortality with overnight recovery ────────────────────────────────

"""
    daily_mortality(tl, T_series_24h; dt_minutes=1.0, recovery_model=nothing)

Fraction of individuals that die during one day's temperature exposure.
"""
function daily_mortality(tl::ToleranceLandscape, T_series_24h;
                         dt_minutes::Real = 1.0, recovery_model=nothing)
    surv = dynamic_survival(tl, T_series_24h; dt_minutes, recovery_model)
    1.0 - surv[end]
end

"""
    cumulative_survival(tl, T_series_per_day; dt_minutes=1.0, recovery_model=nothing)

Cumulative survival probability over multiple days (Rezende 2020 Eq. S9).
`T_series_per_day` is a vector of 24-h temperature vectors (one per day).
Full overnight recovery is assumed between days; pass `recovery_model` for
within-day partial recovery via `dynamic_survival`.
"""
function cumulative_survival(tl::ToleranceLandscape, T_series_per_day::Vector;
                             dt_minutes::Real = 1.0, recovery_model=nothing)
    daily = [last(dynamic_survival(tl, day; dt_minutes, recovery_model))
             for day in T_series_per_day]
    cumprod(daily)
end

# ── TPC ↔ TDT conversion ──────────────────────────────────────────────────────

"""
    tdt_from_tpc(m::UniversalTPCModel; reference_duration=60.0, incipient_temperature)

Convert UTPC parameters to a `LogLinearTDTModel`.
UTPC thermal breadth E and TDT z-value are linked: `z = E × log(10)`.
"""
function tdt_from_tpc(m::UniversalTPCModel;
                      reference_duration::Real    = 60.0,
                      incipient_temperature::Real = (m.T_opt - 273.15) - 20.0)
    z              = m.E * log(10)
    T_opt_C        = m.T_opt - 273.15   # m.T_opt stored as bare K
    reference_ctmax = T_opt_C + m.E
    LogLinearTDTModel(z_value=z, reference_ctmax=reference_ctmax,
                     reference_duration=Float64(reference_duration),
                     incipient_temperature=Float64(incipient_temperature))
end

"""
    thermal_breadth_from_tdt(m::LogLinearTDTModel)

UTPC thermal breadth E from TDT z-value: `E = z / log(10)`.
"""
thermal_breadth_from_tdt(m::LogLinearTDTModel) = m.z_value / log(10)
thermal_breadth_from_tdt(m::ToleranceLandscape) = m.z_value / log(10)
