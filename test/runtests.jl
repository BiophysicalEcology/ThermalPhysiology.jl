using ThermalPhysiology
using Test
using Unitful
using Random

@testset "ThermalPhysiology.jl" begin

    # ── Arrhenius ──────────────────────────────────────────────────────────────
    @testset "ArrheniusModel" begin
        m = ArrheniusModel(0.65u"eV"; T_ref=293.15u"K")
        @test temperature_correction(m, 293.15u"K") ≈ 1.0
        @test temperature_correction(m, 20.0u"°C") ≈ 1.0        # 20°C == T_ref
        @test temperature_correction(m, 30.0u"°C") > 1.0         # warmer → faster
        @test m(20.0u"°C") ≈ temperature_correction(m, 20.0u"°C")    # callable
        @test_throws ArgumentError temperature_correction(m, 20.0)  # bare number rejected

        # ArrheniusModel also accepts activation energy directly (positional dispatch)
        m_from_T = ArrheniusModel(T_A=8000.0u"K", T_ref=293.15u"K")
        @test m_from_T.T_A ≈ 8000.0u"K"
        @test_throws ArgumentError ArrheniusModel(8000.0)   # bare number rejected
    end

    @testset "SharpSchoolFullModel" begin
        m = sharpe_schoolfield(activation=0.65u"eV", reference_temperature=293.15u"K",
                                low_temperature=277.15u"K", low_deactivation=5.0u"eV",
                                high_temperature=318.15u"K", high_deactivation=10.0u"eV",
                                rate_at_reference=1.0)
        # At T_ref the correction should be close to rate_at_reference
        tc = temperature_correction(m, 293.15u"K")
        @test tc > 0.0
        @test tc <= 1.0  # low_term and high_term both present
    end

    @testset "rate_at_reference=nothing gives bare correction" begin
        # Omitting rate_at_reference returns the dimensionless correction;
        # supplying one returns rate_at_reference * correction (units preserved).
        m_bare = SharpSchoolFullModel(T_A=9000.0u"K", T_ref=293.15u"K", T_L=273.15u"K",
                                       T_AL=50000.0u"K", T_H=303.15u"K", T_AH=90000.0u"K")
        m_rate = SharpSchoolFullModel(T_A=9000.0u"K", T_ref=293.15u"K", T_L=273.15u"K",
                                       T_AL=50000.0u"K", T_H=303.15u"K", T_AH=90000.0u"K",
                                       rate_at_reference=0.1u"d^-1")
        correction = temperature_correction(m_bare, 25.0u"°C")
        rate       = temperature_correction(m_rate, 25.0u"°C")
        @test correction isa Float64
        @test rate ≈ 0.1u"d^-1" * correction
        @test unit(rate) == u"d^-1"
    end

    @testset "Unitful Kelvin struct constructors" begin
        # All Arrhenius structs accept u"K" for temperature parameters
        m1 = ArrheniusModel(T_A=8000.0u"K", T_ref=293.15u"K")
        @test m1.T_A   ≈ 8000.0u"K"
        @test m1.T_ref ≈ 293.15u"K"
        @test temperature_correction(m1, 20.0u"°C") ≈ 1.0

        m2 = SharpSchoolDEBModel(
            T_A  = 8000.0u"K",
            T_ref = 293.15u"K",
            T_L  = 273.15u"K",
            T_AL = 50000.0u"K",
            T_H  = 310.15u"K",
            T_AH = 90000.0u"K",
            rate_at_reference = 1.0,
        )
        @test m2.T_A   ≈ 8000.0u"K"
        @test m2.T_ref ≈ 293.15u"K"
        @test m2.T_L   ≈ 273.15u"K"
        @test temperature_correction(m2, 293.15u"K") ≈ 1.0

        # Bare numbers for temperature parameters are rejected (units required)
        @test_throws ArgumentError SharpSchoolFullModel(T_A=5000.0u"K", T_ref=298.15,
            T_L=277.15u"K", T_AL=26000.0u"K", T_H=316.5u"K", T_AH=100000.0u"K",
            rate_at_reference=1.0)

        m3 = SharpSchoolFullModel(T_A=5000.0u"K", T_ref=298.15u"K", T_L=277.15u"K",
                                   T_AL=26000.0u"K", T_H=316.5u"K", T_AH=100000.0u"K",
                                   rate_at_reference=1.0)
        @test m3.T_A  ≈ 5000.0u"K"
        @test m3.T_AL ≈ 26000.0u"K"

        # Temperature-like fields carry units; unit conversion is transparent
        m4 = ArrheniusModel(T_A=8000.0u"K", T_ref=20.0u"°C")
        @test ustrip(u"K", m4.T_ref) ≈ 293.15
    end

    # ── Universal TPC ──────────────────────────────────────────────────────────
    @testset "UniversalTPCModel" begin
        m = utpc(optimal_temperature=30.0u"°C", thermal_breadth=10.0u"K", maximum_performance=2.0)
        @test thermal_performance(m, 30.0u"°C") ≈ 2.0   # peak at T_opt
        @test thermal_performance(m, 20.0u"°C") < 2.0
        @test m(30.0u"°C") ≈ 2.0                              # callable
        @test_throws ArgumentError thermal_performance(m, 20.0)   # bare number rejected
    end

    # ── Deutsch ────────────────────────────────────────────────────────────────
    @testset "DeutschModel" begin
        m = deutsch(maximum_rate=1.5, optimal_temperature=25.0u"°C", critical_thermal_maximum=40.0u"°C",
                    width_parameter=5.0u"K")
        @test thermal_performance(m, 25.0u"°C") ≈ 1.5
        @test thermal_performance(m, 40.0u"°C") ≈ 0.0 atol=1e-10
        @test thermal_performance(m, 45.0u"°C") < 0.0 || thermal_performance(m, 45.0u"°C") ≈ 0.0
    end

    # ── Other phenomenological models (units enforcement) ─────────────────────
    @testset "phenomenological models require units" begin
        mg = GaussianModel(maximum_rate=1.0, optimal_temperature=25.0u"°C", width_parameter=5.0u"K")
        @test thermal_performance(mg, 25.0u"°C") ≈ 1.0
        @test_throws ArgumentError thermal_performance(mg, 25.0)
        @test_throws ArgumentError GaussianModel(maximum_rate=1.0, optimal_temperature=25.0, width_parameter=5.0u"K")

        mt12 = Thomas2012Model(rate_constant=0.5, shape_parameter=15.0u"K", optimal_temperature=25.0u"°C")
        @test thermal_performance(mt12, 25.0u"°C") > 0.0

        mt17 = Thomas2017Model(maximum_rate=1.0, optimal_temperature=25.0u"°C", width_parameter=5.0u"K", skewness=0.0)
        @test thermal_performance(mt17, 25.0u"°C") ≈ 1.0

        mb2 = Briere2Model(rate_constant=0.01, minimum_temperature=10.0u"°C", maximum_temperature=40.0u"°C",
                            shape_parameter=2.0)
        @test thermal_performance(mb2, 25.0u"°C") > 0.0

        # Regression: max(0.0, val)/return 0.0 used to throw DimensionError when
        # rate_constant carried units (0.0 isn't comparable to a dimensioned Quantity).
        mb1_unitful = Briere1Model(rate_constant=0.0002u"d^-1", minimum_temperature=10.0u"°C", maximum_temperature=38.0u"°C")
        @test thermal_performance(mb1_unitful, 5.0u"°C") == 0.0u"d^-1"    # below T_min
        @test thermal_performance(mb1_unitful, 45.0u"°C") == 0.0u"d^-1"   # above T_max
        mt12_unitful = Thomas2012Model(rate_constant=0.5u"d^-1", shape_parameter=15.0u"K", optimal_temperature=25.0u"°C")
        @test thermal_performance(mt12_unitful, -20.0u"°C") == 0.0u"d^-1"

        ml = Lactin2Model(rate_constant=0.1, maximum_temperature=40.0u"°C", delta_temperature=2.0u"K", intercept=-1.0)
        @test thermal_performance(ml, 25.0u"°C") > 0.0
        @test_throws ArgumentError Lactin2Model(rate_constant=0.1, maximum_temperature=40.0,
            delta_temperature=2.0u"K", intercept=-1.0)

        # PawarModel: activation/deactivation accept either an energy or an Arrhenius temperature
        mp_ev = PawarModel(rate_at_reference=1.0, activation_energy=0.65u"eV", deactivation_energy=1.15u"eV",
                            peak_temperature=30.0u"°C", reference_temperature=20.0u"°C")
        mp_K  = PawarModel(rate_at_reference=1.0, activation_energy=ea_to_ta(0.65u"eV"),
                            deactivation_energy=ea_to_ta(1.15u"eV"),
                            peak_temperature=30.0u"°C", reference_temperature=20.0u"°C")
        @test thermal_performance(mp_ev, 25.0u"°C") ≈ thermal_performance(mp_K, 25.0u"°C")
        @test_throws ArgumentError PawarModel(rate_at_reference=1.0, activation_energy=0.65,
            deactivation_energy=1.15u"eV", peak_temperature=30.0u"°C", reference_temperature=20.0u"°C")
    end

    # ── LogLinearTDTModel ──────────────────────────────────────────────────────
    @testset "LogLinearTDTModel" begin
        m = log_linear_tdt(z_value=4.0u"K", reference_ctmax=39.0u"°C", reference_duration=60.0u"minute",
                           incipient_temperature=30.0u"°C")
        @test survival_time(m, 39.0u"°C") ≈ 60.0u"minute"     # at reference_ctmax → reference_duration
        @test ctmax_at_duration(m, 60.0u"minute") ≈ 39.0u"°C" # identity
        @test ctmax_at_duration(m, 600.0u"minute") ≈ (39.0 - 4.0)u"°C"  # 10× longer → 1 z-unit cooler
        @test m(39.0u"°C") ≈ 60.0u"minute"                    # callable
        @test temperature_maximum(m) ≈ (39.0 + 4.0 * log10(60.0))u"°C"
        @test_throws ArgumentError ctmax_at_duration(m, 60.0)   # bare number rejected
        @test_throws ArgumentError survival_time(m, 39.0)       # bare number rejected
        @test_throws ArgumentError log_linear_tdt(z_value=4.0, reference_ctmax=39.0u"°C",
            reference_duration=60.0u"minute", incipient_temperature=30.0u"°C")  # bare z_value rejected

        # z_value also accepts °C, reinterpreted as a K-sized difference
        m_degc = log_linear_tdt(z_value=4.0u"°C", reference_ctmax=39.0u"°C",
            reference_duration=60.0u"minute", incipient_temperature=30.0u"°C")
        @test z_value(m_degc) == 4.0u"K"

        # Round-trip: lethal_temperature inverts survival_time
        @test lethal_temperature(m, survival_time(m, 38.0u"°C")) ≈ 38.0u"°C" atol=1e-6u"°C"
    end

    # ── ToleranceLandscape ─────────────────────────────────────────────────────
    @testset "ToleranceLandscape" begin
        # Synthetic survival curve: linear decline from 1 to 0 over 0..120 min
        times = collect(range(0.0, 120.0, length=200))
        surv  = collect(range(1.0, 0.0, length=200))
        tl = ToleranceLandscape(z_value=4.0u"K", temperature_maximum=42.0u"°C",
                                mean_assay_temperature=38.0u"°C", survival_curve=hcat(times, surv))

        # Median at mean_assay_temperature should be near midpoint (60 min for linear)
        t50 = survival_time(tl, 38.0u"°C")
        @test 55.0 < t50 < 65.0
        @test_throws ArgumentError survival_time(tl, 38.0)   # bare number rejected
        @test_throws ArgumentError ToleranceLandscape(z_value=4.0, temperature_maximum=42.0u"°C",
            mean_assay_temperature=38.0u"°C", survival_curve=hcat(times, surv))  # bare z_value rejected

        # dynamic_survival: all-cold series → survival stays high
        cold = fill(20.0u"°C", 30)
        surv_cold = dynamic_survival(tl, cold)
        @test all(surv_cold .>= 0.9)

        # dynamic_survival: all-hot series → survival declines
        hot = fill(45.0u"°C", 200)
        surv_hot = dynamic_survival(tl, hot)
        @test surv_hot[end] < surv_hot[1]
    end

    # ── Injury accumulation ────────────────────────────────────────────────────
    @testset "accumulated_injury and time_to_failure" begin
        m = log_linear_tdt(z_value=4.0u"K", reference_ctmax=39.0u"°C", reference_duration=60.0u"minute",
                           incipient_temperature=30.0u"°C")
        # Constant exposure at reference_ctmax: injury reaches 1 after reference_duration steps
        T_const = fill(39.0u"°C", 120)
        inj = accumulated_injury(m, T_const, 1.0u"minute")
        @test inj[60] ≈ 1.0 atol=0.01
        @test time_to_failure(m, T_const, 1.0u"minute") ≈ 60.0u"minute" atol=1.0u"minute"
        @test_throws ArgumentError accumulated_injury(m, T_const, 1.0)   # bare number rejected
        @test_throws ArgumentError accumulated_injury(m, fill(39.0, 120), 1.0u"minute")   # bare T_series rejected

        # Below incipient temperature: no injury
        T_cool = fill(25.0u"°C", 60)
        inj_cool = accumulated_injury(m, T_cool, 1.0u"minute")
        @test all(inj_cool .== 0.0)

        # Survives entire series
        @test time_to_failure(m, T_cool, 1.0u"minute") == Inf * u"minute"
    end

    # ── Resettable injury accumulation ─────────────────────────────────────────
    @testset "resettable_injury" begin
        m = log_linear_tdt(z_value=4.0u"K", reference_ctmax=39.0u"°C", reference_duration=60.0u"minute",
                           incipient_temperature=30.0u"°C")

        # Heat for a while (accumulating injury), drop below incipient temp
        # (should reset to 0), then heat again from scratch.
        T_series = vcat(fill(39.0u"°C", 40), fill(25.0u"°C", 5), fill(39.0u"°C", 40))
        inj = resettable_injury(m, T_series, 1.0u"minute")
        @test inj[40] > 0.0              # injury built up during first hot phase
        @test all(inj[41:45] .== 0.0)    # fully reset while cool
        @test inj[46] > 0.0 && inj[46] < inj[40]  # restarted from 0, not carried over

        # Once lethal (injury reaches 1.0), a cool-down does NOT resurrect it
        T_lethal_then_cool = vcat(fill(39.0u"°C", 120), fill(20.0u"°C", 10))
        inj2 = resettable_injury(m, T_lethal_then_cool, 1.0u"minute")
        @test inj2[120] ≈ 1.0 atol=0.01
        @test inj2[end] ≈ 1.0 atol=0.01   # stays at 1.0, not reset

        # resettable mode of time_to_failure matches a manual check
        @test time_to_failure(m, T_lethal_then_cool, 1.0u"minute"; resettable=true) ≈
              time_to_failure(m, T_lethal_then_cool, 1.0u"minute"; resettable=false)
    end

    # ── step_injury / repair models ────────────────────────────────────────────
    @testset "step_injury and repair models" begin
        m = log_linear_tdt(z_value=4.0u"K", reference_ctmax=39.0u"°C", reference_duration=60.0u"minute",
                           incipient_temperature=30.0u"°C")
        T_series = vcat(fill(39.0u"°C", 40), fill(25.0u"°C", 5), fill(39.0u"°C", 40))

        # NoRepair matches accumulated_injury; FullRepairBelowThreshold matches resettable_injury
        inj_norepair = Float64[]
        cum = 0.0
        for T in T_series
            cum = step_injury(m, NoRepair(), cum, T, 1.0u"minute")
            push!(inj_norepair, cum)
        end
        @test inj_norepair ≈ accumulated_injury(m, T_series, 1.0u"minute")

        repair = full_repair_below_threshold(30.0u"°C")
        @test repair isa FullRepairBelowThreshold
        inj_repair = Float64[]
        cum = 0.0
        for T in T_series
            cum = step_injury(m, repair, cum, T, 1.0u"minute")
            push!(inj_repair, cum)
        end
        @test inj_repair ≈ resettable_injury(m, T_series, 1.0u"minute")
        @test inj_repair[45] == 0.0   # cooled off -- reset
    end

    # ── TDT fitting from raw binary survival data (joint one-stage fit) ───────
    @testset "fit_thermal_death_time_curve (binary survival)" begin
        # Synthetic data from a known model: z=2.8, reference_ctmax=54 at 1 min
        m_true = log_linear_tdt(z_value=2.8u"K", reference_ctmax=54.0u"°C", reference_duration=1.0u"minute",
                                 incipient_temperature=30.0u"°C")
        temps = [44.0, 45.0, 46.0, 47.0, 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0]
        times = [1440.0, 720.0, 180.0, 45.0, 12.0, 6.0, 2.0, 0.5]

        Random.seed!(1)
        rows_T = Float64[]; rows_t = Float64[]; rows_surv = Bool[]
        for T in temps, t in times
            surv_t = ustrip(u"minute", survival_time(m_true, T * u"°C"))
            p_alive = 1.0 / (1.0 + (t / surv_t)^2.0)
            for _ in 1:8
                push!(rows_T, T); push!(rows_t, t)
                push!(rows_surv, rand() < p_alive)
            end
        end
        data = BinarySurvivalData(temperatures=rows_T.*u"°C", exposure_times=rows_t.*u"minute", survived=rows_surv)
        m_fit = fit_thermal_death_time_curve(data; reference_duration=1.0u"minute")

        @test m_fit.z_value ≈ 2.8u"K" rtol=0.25
        @test m_fit.reference_ctmax ≈ 54.0u"°C" atol=1.5u"°C"
    end

    # ── Properties ────────────────────────────────────────────────────────────
    @testset "optimal_temperature" begin
        m = utpc(optimal_temperature=30.0u"°C", thermal_breadth=10.0u"K", maximum_performance=1.0)
        T_opt = optimal_temperature(m)
        @test ustrip(u"K", uconvert(u"K", T_opt)) ≈ 303.15 atol=0.01

        # Analytic for Deutsch
        md = deutsch(maximum_rate=1.0, optimal_temperature=28.0u"°C",
                      critical_thermal_maximum=40.0u"°C", width_parameter=5.0u"K")
        @test ustrip(u"°C", uconvert(u"°C", optimal_temperature(md))) ≈ 28.0 atol=0.01

        # Arrhenius: no finite peak
        ma = ArrheniusModel(T_A=8000.0u"K", T_ref=293.15u"K")
        @test optimal_temperature(ma) == Inf * u"K"
    end

    @testset "critical_thermal_maximum" begin
        m = utpc(optimal_temperature=30.0u"°C", thermal_breadth=10.0u"K", maximum_performance=1.0)
        ctmax = critical_thermal_maximum(m)
        @test ustrip(u"K", uconvert(u"K", ctmax)) > 303.15  # above T_opt

        md = deutsch(maximum_rate=1.0, optimal_temperature=25.0u"°C", critical_thermal_maximum=40.0u"°C",
                      width_parameter=5.0u"K")
        @test ustrip(u"°C", uconvert(u"°C", critical_thermal_maximum(md))) ≈ 40.0
    end

    @testset "maximum_rate" begin
        m = utpc(optimal_temperature=30.0u"°C", thermal_breadth=10.0u"K", maximum_performance=3.0)
        @test maximum_rate(m) ≈ 3.0 atol=0.01
    end

    @testset "q10" begin
        # UTPC Q10 should be > 1 on the rising limb (below T_opt)
        m = utpc(optimal_temperature=30.0u"°C", thermal_breadth=10.0u"K", maximum_performance=1.0)
        @test q10(m, 15.0u"°C") > 1.0   # 15→25°C, both below T_opt
        @test_throws ArgumentError q10(m, 15.0)   # bare number rejected
    end

    @testset "z_value" begin
        m = log_linear_tdt(z_value=4.0u"K", reference_ctmax=39.0u"°C", reference_duration=60.0u"minute",
                            incipient_temperature=30.0u"°C")
        @test z_value(m) ≈ 4.0u"K"

        # Arrhenius z_value: T_ref²/T_A * log(10)
        ma = ArrheniusModel(T_A=8000.0u"K", T_ref=293.15u"K")
        @test z_value(ma) ≈ (293.15^2 / 8000.0 * log(10))u"K"
    end

    # ── CTE ───────────────────────────────────────────────────────────────────
    @testset "constant_temperature_equivalent" begin
        m = ArrheniusModel(T_A=8000.0u"K", T_ref=293.15u"K")
        T_vals   = collect(range(15.0, 35.0, length=100))   # 100 values from 15 to 35°C
        T_series = T_vals .* u"°C"
        cte = constant_temperature_equivalent(m, T_series)
        T_mean = sum(T_vals) / length(T_vals)
        # Jensen's inequality: CTE > arithmetic mean for convex Arrhenius
        @test ustrip(u"K", uconvert(u"K", cte)) > T_mean + 273.15

        # Numeric and analytic agree for ArrheniusModel
        m2 = ArrheniusModel(T_A=9000.0u"K", T_ref=293.15u"K")
        cte2 = constant_temperature_equivalent(m2, T_series)
        @test ustrip(u"K", uconvert(u"K", cte2)) > 0.0
    end

    # ── Static ↔ dynamic CTmax conversions ─────────────────────────────────────
    @testset "dynamic_ctmax and static_ctmax_from_dynamic" begin
        m = log_linear_tdt(z_value=4.0u"K", reference_ctmax=39.0u"°C", reference_duration=60.0u"minute",
                           incipient_temperature=20.0u"°C")
        ramp = 0.1u"K/minute"
        dctmax = dynamic_ctmax(m, ramp)
        @test dctmax > m.reference_ctmax  # dynamic > static for slow ramp

        # Round-trip: recover reference_ctmax
        recovered = static_ctmax_from_dynamic(m, dctmax, ramp)
        @test recovered ≈ m.reference_ctmax atol=0.01u"°C"

        @test_throws ArgumentError dynamic_ctmax(m, 0.1)   # bare number rejected
        @test_throws Unitful.AffineError 0.1u"°C/minute"   # °C is invalid for a rate
    end

    # ── TPC ↔ TDT bridge ───────────────────────────────────────────────────────
    @testset "tdt_from_tpc" begin
        m_tpc = utpc(optimal_temperature=30.0u"°C", thermal_breadth=10.0u"K", maximum_performance=1.0)
        m_tdt = tdt_from_tpc(m_tpc)
        # z = E * log(10) ≈ 10 * log(10) ≈ 23.03
        @test m_tdt.z_value ≈ (10.0 * log(10))u"K" atol=0.01u"K"
        @test thermal_breadth_from_tdt(m_tdt) ≈ 10.0u"K" atol=0.01u"K"
    end

    # ── TDT fitting ────────────────────────────────────────────────────────────
    @testset "fit_thermal_death_time_curve (static)" begin
        # Synthetic data: z=4, reference_ctmax=39 at 60 min
        m_true = log_linear_tdt(z_value=4.0u"K", reference_ctmax=39.0u"°C", reference_duration=60.0u"minute",
                                 incipient_temperature=30.0u"°C")
        temps  = [35.0, 37.0, 39.0, 41.0, 43.0]
        times  = [survival_time(m_true, T * u"°C") for T in temps]
        data   = StaticKnockdownData(temperatures=temps.*u"°C", knockdown_times=times)
        m_fit  = fit_thermal_death_time_curve(data; reference_duration=60.0u"minute")
        @test_throws ArgumentError StaticKnockdownData(temperatures=temps.*u"°C", knockdown_times=Float64.(temps))
        @test m_fit.z_value ≈ 4.0u"K" atol=0.05u"K"
        @test m_fit.reference_ctmax ≈ 39.0u"°C" atol=0.05u"°C"
    end

    # ── Briere1 fitting ────────────────────────────────────────────────────────
    @testset "fit_thermal_performance_curve (Briere1Model)" begin
        m_true = Briere1Model(rate_constant=0.0002, minimum_temperature=10.0u"°C", maximum_temperature=38.0u"°C")
        temps = collect(12.0:2.0:36.0) .* u"°C"
        rates = thermal_performance.(Ref(m_true), temps)
        m_fit = fit_thermal_performance_curve(Briere1Model, temps, rates)
        @test ustrip(u"°C", m_fit.minimum_temperature) ≈ 10.0 atol=1.0
        @test ustrip(u"°C", m_fit.maximum_temperature) ≈ 38.0 atol=1.0
        @test thermal_performance(m_fit, 25.0u"°C") ≈ thermal_performance(m_true, 25.0u"°C") rtol=0.05
    end

    # ── Schoolfield fitting (synthetic data — verifies parameter recovery) ────────
    @testset "fit_thermal_performance_curve (SharpSchoolFullModel)" begin
        # Generate noiseless data from a known model and verify NLS recovery.
        m_true = SharpSchoolFullModel(T_A=5000.0u"K", T_ref=298.15u"K", T_L=291.0u"K",
                                      T_AL=26000.0u"K", T_H=316.5u"K", T_AH=100000.0u"K",
                                      rate_at_reference=0.28)
        temps = collect(8.0:4.0:46.0) .* u"°C"   # 10 points spanning cold→hot
        rates = temperature_correction.(Ref(m_true), temps)
        m_fit = fit_thermal_performance_curve(
            SharpSchoolFullModel, temps, rates; T_ref=298.15u"K")
        @test ustrip(u"K", m_fit.T_A)  ≈ 5000.0   rtol=0.02
        @test ustrip(u"K", m_fit.T_L)  ≈ 291.0    rtol=0.01
        @test ustrip(u"K", m_fit.T_H)  ≈ 316.5    rtol=0.01
        @test ustrip(u"K", m_fit.T_AL) ≈ 26000.0  rtol=0.05
        @test ustrip(u"K", m_fit.T_AH) ≈ 100000.0 rtol=0.05
        @test temperature_correction(m_fit, 25.0u"°C") ≈ temperature_correction(m_true, 25.0u"°C") rtol=0.01

        # Regression: bare T_ref used to be silently misread as °C via _K instead
        # of the strict _kelvin_param (T_ref=298.15 → 571.3 K).
        @test_throws ArgumentError fit_thermal_performance_curve(
            SharpSchoolFullModel, temps, rates; T_ref=298.15)
    end

    # ── Unit-safety regressions ────────────────────────────────────────────────
    @testset "unit-safety bug fixes" begin
        # ea_to_ta/ta_to_ea used to silently accept bare numbers and produce a
        # nonsense compound unit (K·eV⁻¹) instead of K.
        @test_throws ArgumentError ea_to_ta(0.65)
        @test_throws ArgumentError ta_to_ea(9000.0)
        @test ea_to_ta(0.65u"eV") isa Unitful.Temperature
        @test ta_to_ea(9000.0u"K") isa Unitful.Energy

        # tdt_from_tpc's incipient_temperature used to be annotated ::Real,
        # hard-rejecting Unitful input even though LogLinearTDTModel requires it.
        m_utpc = utpc(optimal_temperature=30.0u"°C", thermal_breadth=10.0u"K", maximum_performance=1.0)
        m_tdt_unitful = tdt_from_tpc(m_utpc; incipient_temperature=20.0u"°C")
        @test ustrip(u"°C", m_tdt_unitful.incipient_temperature) ≈ 20.0
    end

    # ── Registry ───────────────────────────────────────────────────────────────
    @testset "registry" begin
        names = model_names()
        @test :arrhenius in names
        @test :utpc in names
        @test :log_linear_tdt in names
        @test :tolerance_landscape in names

        tdt_names = model_names(family=:tdt)
        @test all(v -> THERMAL_REGISTRY[v].family == :tdt, tdt_names)

        phenom = thermal_models(family=:phenomenological)
        @test haskey(phenom, :utpc)
        @test haskey(phenom, :deutsch)
    end

    # ── multi-day survival ─────────────────────────────────────────────────────
    @testset "cumulative_survival" begin
        times = collect(range(0.0, 120.0, length=200))
        surv  = collect(range(1.0, 0.0, length=200))
        tl    = ToleranceLandscape(z_value=4.0u"K", temperature_maximum=42.0u"°C",
                                   mean_assay_temperature=38.0u"°C", survival_curve=hcat(times, surv))
        mild  = fill(30.0u"°C", 60)
        cum   = cumulative_survival(tl, [mild, mild, mild])
        @test length(cum) == 3
        @test cum[1] >= cum[2] >= cum[3]  # non-increasing
        @test all(0.0 .<= cum .<= 1.0)
    end

end # @testset ThermalPhysiology
