# ── ThermalPhysiology.jl — worked examples ────────────────────────────────────
#
# Run from the package root with:
#   julia --project=. examples/tpc_tdt_overview.jl

using ThermalPhysiology
using Unitful

# ─────────────────────────────────────────────────────────────────────────────
# 1. Thermal performance curves (TPC)
# ─────────────────────────────────────────────────────────────────────────────

println("=== 1. Thermal Performance Curves ===\n")

# Universal TPC (Arnoldi et al. 2025): two-parameter, scale-invariant
# T_opt in °C; E (thermal breadth) in K
m_utpc = utpc(optimal_temperature=30.0u"°C", thermal_breadth=10.0u"K",
              maximum_performance=2.5)
println("Universal TPC at T_opt (30°C): ", thermal_performance(m_utpc, 30.0))
println("Universal TPC at 40°C:         ", round(thermal_performance(m_utpc, 40.0), digits=4))
println("Universal TPC at 20°C:         ", round(thermal_performance(m_utpc, 20.0), digits=4))

# Callable as a function
temps_C = 0.0:1.0:50.0
perf_utpc = m_utpc.(temps_C)
println("Peak performance (via maximum_rate): ", round(maximum_rate(m_utpc), digits=4))
println("Optimal temperature:                 ", optimal_temperature(m_utpc))
println("CTmax (performance → 0):             ", round(ustrip(u"°C", uconvert(u"°C", critical_thermal_maximum(m_utpc))), digits=2), " °C")
println("CTmin:                               ", round(ustrip(u"°C", uconvert(u"°C", critical_thermal_minimum(m_utpc))), digits=2), " °C")
println("Thermal breadth (CTmax - CTmin):     ", round(ustrip(u"K", thermal_breadth(m_utpc)), digits=2), " K")
println()

# Arrhenius temperature-correction model (DEBtool convention)
# Returns dimensionless factor = 1 at T_ref
m_arr = arrhenius(activation_energy=0.65u"eV", reference_temperature=20.0u"°C")
println("Arrhenius correction at T_ref (20°C): ", temperature_correction(m_arr, 20.0))
println("Arrhenius correction at 30°C:          ", round(temperature_correction(m_arr, 30.0), digits=4))
println("Arrhenius correction at 10°C:          ", round(temperature_correction(m_arr, 10.0), digits=4))
println()

# Sharpe-Schoolfield (full model: low + high deactivation)
m_ss = sharpe_schoolfield(
    activation_energy        = 0.65u"eV",
    reference_temperature    = 20.0u"°C",
    low_temperature          = 0.0u"°C",
    low_deactivation_energy  = 1.5u"eV",
    high_temperature         = 42.0u"°C",
    high_deactivation_energy = 5.0u"eV",
    rate_at_reference        = 1.0,
)
println("Sharpe-Schoolfield at 20°C: ", round(temperature_correction(m_ss, 20.0), digits=4))
println("Sharpe-Schoolfield at 40°C: ", round(temperature_correction(m_ss, 40.0), digits=4))
println("Sharpe-Schoolfield at 45°C: ", round(temperature_correction(m_ss, 45.0), digits=4), " (collapsed)")
println()

# Deutsch 2008 — widely used in macrophysiology
m_deutsch = deutsch(maximum_rate=1.0, optimal_temperature=25.0,
                    critical_thermal_maximum=40.0, width_parameter=6.0)
println("Deutsch at CTmax (40°C):  ", thermal_performance(m_deutsch, 40.0))
println("Deutsch at T_opt (25°C):  ", thermal_performance(m_deutsch, 25.0))
println()

# Pawar 2018 — activation/deactivation energies in eV
m_pawar = pawar(activation_energy=0.65, deactivation_energy=1.15,
                peak_temperature=32.0, reference_temperature=20.0)
println("Pawar at peak:  ", round(thermal_performance(m_pawar, 32.0), digits=4))
println("Pawar at 20°C:  ", round(thermal_performance(m_pawar, 20.0), digits=4))
println()

# ─────────────────────────────────────────────────────────────────────────────
# 2. Q10 and Constant Temperature Equivalent
# ─────────────────────────────────────────────────────────────────────────────

println("=== 2. Q10 and Constant Temperature Equivalent ===\n")

println("Q10 for Arrhenius at 10°C:  ", round(q10(m_utpc, 15.0), digits=3))
println("  (ratio of performance at 25°C vs 15°C)")

# Constant temperature equivalent: Jensen's inequality means CTE > arithmetic mean
T_series = collect(range(15.0, 35.0, length=1000))  # sine-like range
cte = constant_temperature_equivalent(m_arr, T_series)
T_arith = sum(T_series) / length(T_series)           # arithmetic mean = 25°C
println()
println("Temperature series: 15–35°C (mean = $(T_arith)°C)")
println("Constant temperature equivalent (Arrhenius):  ",
        round(ustrip(u"°C", uconvert(u"°C", cte)), digits=3), " °C")
println("  (CTE > 25°C because Arrhenius is convex — Jensen's inequality)")
println("Mean correction factor: ", round(mean_correction_factor(m_arr, T_series), digits=4))
println()

# ─────────────────────────────────────────────────────────────────────────────
# 3. Thermal death time (TDT) — LogLinearTDTModel
# ─────────────────────────────────────────────────────────────────────────────

println("=== 3. Thermal Death Time — Log-linear TDT ===\n")

# Drosophila melanogaster parameters (Jørgensen et al. 2021 Fig. 3A)
m_tdt = log_linear_tdt(
    z_value              = 4.26,   # °C per decade change in knockdown time
    reference_ctmax      = 39.1,   # sCTmax at 60-min exposure (°C)
    reference_duration   = 60.0,   # minutes
    incipient_temperature = 30.0,  # °C below which no injury accumulates
)

println("Survival time at 39.1°C (reference CTmax): ",
        round(survival_time(m_tdt, 39.1), digits=1), " min")
println("Survival time at 35°C:  ", round(survival_time(m_tdt, 35.0), digits=0), " min")
println("Survival time at 43°C:  ", round(survival_time(m_tdt, 43.0), digits=1), " min")
println("z-value:                ", z_value(m_tdt), " °C/decade")
println("T_max (mean τ = 1 min): ", round(temperature_maximum(m_tdt), digits=2), " °C")
println()

# CTmax at various exposure durations (TDT curve)
println("sCTmax at 5 min:   ", round(ctmax_at_duration(m_tdt, 5.0), digits=2), " °C")
println("sCTmax at 60 min:  ", round(ctmax_at_duration(m_tdt, 60.0), digits=2), " °C  (reference)")
println("sCTmax at 600 min: ", round(ctmax_at_duration(m_tdt, 600.0), digits=2), " °C")
println()

# Dynamic CTmax prediction (ramping assay, Jørgensen 2021 Eq. 7a)
println("Predicted dCTmax at 0.5°C/min ramp: ",
        round(dynamic_ctmax(m_tdt, 0.5/60.0), digits=2), " °C")  # 0.5°C/min in °C/min
println("Predicted dCTmax at 0.1°C/min ramp: ",
        round(dynamic_ctmax(m_tdt, 0.1/60.0), digits=2), " °C")

# Recover sCTmax from a dynamic assay
dctmax_obs = dynamic_ctmax(m_tdt, 0.5/60.0)
recovered_sctmax = static_ctmax_from_dynamic(m_tdt, dctmax_obs, 0.5/60.0)
println("Round-trip static→dynamic→static: recovered sCTmax = ",
        round(recovered_sctmax, digits=3), " °C (true = 39.1°C)")
println()

# Lethal temperature for a given exposure duration
LT50_10min = lethal_temperature(m_tdt, 10.0)
println("LT50 at 10 min exposure: ", round(LT50_10min, digits=2), " °C")
println()

# ─────────────────────────────────────────────────────────────────────────────
# 4. Injury accumulation under fluctuating temperatures
# ─────────────────────────────────────────────────────────────────────────────

println("=== 4. Injury Accumulation (Jørgensen 2021) ===\n")

# Simulate a 3-hour temperature exposure: 2h cool, then 1h hot
T_cool = fill(32.0, 120)   # 120 min at 32°C — below sCTmax, some injury
T_hot  = fill(41.0, 60)    # 60 min at 41°C — above sCTmax
T_exposure = vcat(T_cool, T_hot)

injury = accumulated_injury(m_tdt, T_exposure, 1.0)   # dt = 1 min
ttf    = time_to_failure(m_tdt, T_exposure, 1.0)

println("3-h exposure: 120 min @ 32°C, then 60 min @ 41°C")
println("  Injury after 120 min (32°C phase): ", round(injury[120], digits=4))
println("  Injury after 180 min (41°C phase): ", round(injury[180], digits=4))
println("  Time to failure: ", ttf == Inf ? "survives" : "$(round(ttf, digits=0)) min")
println()

# Pure cool exposure: survives
T_safe = fill(28.0, 240)
ttf_safe = time_to_failure(m_tdt, T_safe, 1.0)
println("240 min @ 28°C: ", ttf_safe == Inf ? "survives (no injury below incipient temp)" : "fails")
println()

# ─────────────────────────────────────────────────────────────────────────────
# 5. ToleranceLandscape — distributional TDT (Rezende 2020)
# ─────────────────────────────────────────────────────────────────────────────

println("=== 5. Tolerance Landscape (Rezende 2020) ===\n")

# Build a synthetic tolerance landscape from individual-level knockdown data
# Simulate 200 individuals at 4 assay temperatures (50 each)
using Random
Random.seed!(42)

# true_z=4, T_max=42°C (temperature where median knockdown time = 1 min)
# At reference_ctmax=35°C, survival time = 10^((42-35)/4) = 10^1.75 ≈ 56 min
true_z    = 4.0
true_tmax = 42.0    # T_max: where mean τ = 1 min (Rezende's CTmax definition)
assay_temps = repeat([32.0, 34.0, 36.0, 38.0], inner=50)

# True median survival time at each temperature: t = 1 * 10^((T_max - T)/z)
t_medians = [1.0 * 10^((true_tmax - T) / true_z) for T in assay_temps]
# Individual variation: log-normal spread
knockdown_times = t_medians .* exp.(0.3 .* randn(200))

data = IndividualKnockdownData(temperatures=assay_temps, knockdown_times=knockdown_times)
tl   = fit_tolerance_landscape(data; n_bins=500)

println("Fitted z-value:              ", round(tl.z_value, digits=2), " °C  (true = 4.0)")
println("Fitted T_max:                ", round(tl.temperature_maximum, digits=2), " °C  (true = 42.0)")
println("Mean assay temperature:      ", round(tl.mean_assay_temperature, digits=2), " °C")
println("Survival curve rows:         ", size(tl.survival_curve, 1))
println()

# Median knockdown time at assay temperatures
for T in [32.0, 34.0, 36.0, 38.0]
    true_t50 = 1.0 * 10^((true_tmax - T) / true_z)
    fitted_t50 = survival_time(tl, T)
    println("  T = $(T)°C: fitted τ₅₀ = $(round(fitted_t50, digits=1)) min  ",
            "(true = $(round(true_t50, digits=1)) min)")
end
println()

# Dynamic survival prediction under a fluctuating temperature series
# Simulate a warm afternoon (temperature ramps up and back)
n_hours = 8
t_minutes = 0:60:(n_hours*60 - 60)
T_daily = 30.0 .+ 8.0 .* sin.(π .* (0:length(t_minutes)-1) ./ (length(t_minutes)-1))
surv_vec = dynamic_survival(tl, T_daily; dt_minutes=60.0)

println("Dynamic survival over $(n_hours)h warm afternoon (peak $(round(maximum(T_daily),digits=1))°C):")
for (i, (t, T, s)) in enumerate(zip(t_minutes, T_daily, surv_vec))
    println("  t=$(lpad(t,3))min  T=$(round(T,digits=1))°C  survival=$(round(s, digits=3))")
end
println()

# ─────────────────────────────────────────────────────────────────────────────
# 6. TPC ↔ TDT bridge
# ─────────────────────────────────────────────────────────────────────────────

println("=== 6. TPC ↔ TDT Bridge ===\n")

# Convert a Universal TPC to a TDT model
m_tpc_bridge = utpc(optimal_temperature=30.0u"°C", thermal_breadth=10.0u"K")
m_tdt_bridge = tdt_from_tpc(m_tpc_bridge; reference_duration=60.0)

println("UTPC thermal breadth E: ", m_tpc_bridge.E, " K")
println("TDT z-value:            ", round(m_tdt_bridge.z_value, digits=3), " °C")
println("  (z = E × log(10) = 10 × 2.303 = $(round(10.0 * log(10), digits=3)))")
println("TDT reference CTmax:    ", round(m_tdt_bridge.reference_ctmax, digits=2),
        " °C  (T_opt + E = 30 + 10 = 40°C)")
println()
println("Recover E from TDT:     ", round(thermal_breadth_from_tdt(m_tdt_bridge), digits=3), " K")
println()

# ─────────────────────────────────────────────────────────────────────────────
# 7. Fitting TPC and TDT models to data
# ─────────────────────────────────────────────────────────────────────────────

println("=== 7. Curve Fitting ===\n")

# --- Fit a UTPC to synthetic TPC data ---
m_true_tpc = utpc(optimal_temperature=28.0u"°C", thermal_breadth=12.0u"K",
                  maximum_performance=3.0)
tpc_temps = collect(5.0:5.0:50.0)
tpc_rates = m_true_tpc.(tpc_temps) .+ 0.02 .* randn(length(tpc_temps))
tpc_rates = max.(tpc_rates, 0.0)

m_fit_tpc = fit_thermal_performance_curve(UniversalTPCModel, tpc_temps, tpc_rates)
println("UTPC fit to synthetic data (true: T_opt=28°C, E=12K, P_max=3.0):")
println("  Fitted T_opt: ", round(m_fit_tpc.T_opt - 273.15, digits=2), " °C")
println("  Fitted E:     ", round(m_fit_tpc.E, digits=2), " K")
println("  Fitted P_max: ", round(m_fit_tpc.maximum_performance, digits=2))
println()

# --- Fit a LogLinearTDTModel from static knockdown data ---
m_true_tdt = log_linear_tdt(z_value=4.0, reference_ctmax=39.0, reference_duration=60.0)
static_temps = [35.0, 37.0, 39.0, 41.0, 43.0]
static_times = [survival_time(m_true_tdt, T) for T in static_temps]

data_static = StaticKnockdownData(temperatures=static_temps, knockdown_times=static_times)
m_fit_tdt = fit_thermal_death_time_curve(data_static; reference_duration=60.0)

println("LogLinearTDT fit to static data (true: z=4.0, CTmax=39.0):")
println("  Fitted z:      ", round(m_fit_tdt.z_value, digits=3))
println("  Fitted CTmax:  ", round(m_fit_tdt.reference_ctmax, digits=3))
println()

# --- Fit Sharpe-Schoolfield to E. coli data (Schoolfield et al. 1981, Table 1) ---
# O'Donovan et al. (1965) E. coli growth rates; original paper used weights = 1/rate
# Expected: T_A≈5015K, T_L≈291.2K, T_AL≈25924K, T_H≈316.4K, T_AH≈107700K, ρ(25°C)≈0.273

ecoli_temps = [44.56, 42.32, 39.12, 36.93, 29.64, 25.42, 21.67, 19.05, 13.76, 10.38]
ecoli_rates = [0.2397, 0.5726, 0.5779, 0.5934, 0.3580, 0.2516, 0.2115, 0.1231, 0.0416, 0.0151]

m_ecoli = fit_thermal_performance_curve(
    SharpSchoolFullModel, ecoli_temps, ecoli_rates;
    T_ref = 298.15u"K",
)   # log_transform=true (default): minimises sum of squared log-residuals

# Note: Schoolfield et al. (1981) Table 1 reports T_A=5015K, T_L=291.2K, T_AL=25924K,
# T_H=316.4K, T_AH=107700K from graphical initial estimates + Marquardt.  Our log-space
# NLS finds a lower-SSR solution at slightly different parameter values — both are valid
# fits; the difference reflects a shallower optimum landscape for T_AL and T_AH.
println("Sharpe-Schoolfield fit to E. coli data (O'Donovan 1965; cf. Schoolfield 1981 Table 1):")
println("  T_A  = $(round(m_ecoli.T_A,  digits=0)) K    (Table 1 ref: 5015 K)")
println("  T_L  = $(round(m_ecoli.T_L,  digits=1)) K    (Table 1 ref: 291.2 K)")
println("  T_AL = $(round(m_ecoli.T_AL, digits=0)) K    (Table 1 ref: 25924 K)")
println("  T_H  = $(round(m_ecoli.T_H,  digits=1)) K    (Table 1 ref: 316.4 K)")
println("  T_AH = $(round(m_ecoli.T_AH, digits=0)) K    (Table 1 ref: 107700 K)")
println("  rate at 25°C = $(round(temperature_correction(m_ecoli, 25.0), digits=4))")
println()

# ─────────────────────────────────────────────────────────────────────────────
# 8. Model registry
# ─────────────────────────────────────────────────────────────────────────────

println("=== 8. Model Registry ===\n")

println("All registered models:")
for k in model_names()
    e = THERMAL_REGISTRY[k]
    println("  :$(lpad(string(k), 24))  [$(e.family)]  $(e.name)")
end
println()
println("TDT models only: ", model_names(family=:tdt))
println("Arrhenius family: ", model_names(family=:arrhenius))
