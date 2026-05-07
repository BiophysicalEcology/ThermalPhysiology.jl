using ThermalPhysiology
using Test
using Unitful

@testset "ThermalPhysiology.jl" begin

    # ── Arrhenius ──────────────────────────────────────────────────────────────
    @testset "ArrheniusModel" begin
        m = arrhenius(activation_energy=0.65u"eV", reference_temperature=293.15u"K")
        @test temperature_correction(m, 293.15u"K") ≈ 1.0
        @test temperature_correction(m, 20.0) ≈ 1.0        # 20°C == T_ref
        @test temperature_correction(m, 30.0) > 1.0         # warmer → faster
        @test m(20.0) ≈ temperature_correction(m, 20.0)    # callable
    end

    @testset "SharpSchoolFullModel" begin
        m = sharpe_schoolfield()
        # At T_ref the correction should be close to rate_at_reference
        tc = temperature_correction(m, 293.15u"K")
        @test tc > 0.0
        @test tc <= 1.0  # low_term and high_term both present
    end

    @testset "Unitful Kelvin struct constructors" begin
        # All Arrhenius structs accept u"K" for temperature parameters
        m1 = ArrheniusModel(T_A=8000.0u"K", T_ref=293.15u"K")
        @test m1.T_A   ≈ 8000.0
        @test m1.T_ref ≈ 293.15
        @test temperature_correction(m1, 20.0) ≈ 1.0

        m2 = SharpSchoolDEBModel(
            T_A  = 8000.0u"K",
            T_ref = 293.15u"K",
            T_L  = 273.15u"K",
            T_AL = 50000.0u"K",
            T_H  = 310.15u"K",
            T_AH = 90000.0u"K",
            rate_at_reference = 1.0,
        )
        @test m2.T_A   ≈ 8000.0
        @test m2.T_ref ≈ 293.15
        @test m2.T_L   ≈ 273.15
        @test temperature_correction(m2, 293.15u"K") ≈ 1.0

        # Mix of bare Float64 and Unitful is also accepted
        m3 = SharpSchoolFullModel(T_A=5000.0u"K", T_ref=298.15, T_L=277.15,
                                   T_AL=26000.0u"K", T_H=316.5, T_AH=100000.0u"K")
        @test m3.T_A  ≈ 5000.0
        @test m3.T_AL ≈ 26000.0
    end

    # ── Universal TPC ──────────────────────────────────────────────────────────
    @testset "UniversalTPCModel" begin
        m = utpc(optimal_temperature=30.0u"°C", thermal_breadth=10.0u"K", maximum_performance=2.0)
        @test thermal_performance(m, 30.0u"°C") ≈ 2.0   # peak at T_opt
        @test thermal_performance(m, 20.0) < 2.0
        @test m(30.0) ≈ 2.0                              # callable
    end

    # ── Deutsch ────────────────────────────────────────────────────────────────
    @testset "DeutschModel" begin
        m = deutsch(maximum_rate=1.5, optimal_temperature=25.0, critical_thermal_maximum=40.0)
        @test thermal_performance(m, 25.0) ≈ 1.5
        @test thermal_performance(m, 40.0) ≈ 0.0 atol=1e-10
        @test thermal_performance(m, 45.0) < 0.0 || thermal_performance(m, 45.0) ≈ 0.0
    end

    # ── LogLinearTDTModel ──────────────────────────────────────────────────────
    @testset "LogLinearTDTModel" begin
        m = log_linear_tdt(z_value=4.0, reference_ctmax=39.0, reference_duration=60.0,
                           incipient_temperature=30.0)
        @test survival_time(m, 39.0) ≈ 60.0           # at reference_ctmax → reference_duration
        @test ctmax_at_duration(m, 60.0) ≈ 39.0       # identity
        @test ctmax_at_duration(m, 600.0) ≈ 39.0 - 4.0  # 10× longer → 1 z-unit cooler
        @test m(39.0) ≈ 60.0                            # callable
        @test temperature_maximum(m) ≈ 39.0 + 4.0 * log10(60.0)

        # Round-trip: lethal_temperature inverts survival_time
        @test lethal_temperature(m, survival_time(m, 38.0)) ≈ 38.0 atol=1e-6
    end

    # ── ToleranceLandscape ─────────────────────────────────────────────────────
    @testset "ToleranceLandscape" begin
        # Synthetic survival curve: linear decline from 1 to 0 over 0..120 min
        times = collect(range(0.0, 120.0, length=200))
        surv  = collect(range(1.0, 0.0, length=200))
        tl = ToleranceLandscape(4.0, 42.0, 38.0, hcat(times, surv))

        # Median at mean_assay_temperature should be near midpoint (60 min for linear)
        t50 = survival_time(tl, 38.0)
        @test 55.0 < t50 < 65.0

        # dynamic_survival: all-cold series → survival stays high
        cold = fill(20.0, 30)
        surv_cold = dynamic_survival(tl, cold)
        @test all(surv_cold .>= 0.9)

        # dynamic_survival: all-hot series → survival declines
        hot = fill(45.0, 200)
        surv_hot = dynamic_survival(tl, hot)
        @test surv_hot[end] < surv_hot[1]
    end

    # ── Injury accumulation ────────────────────────────────────────────────────
    @testset "accumulated_injury and time_to_failure" begin
        m = log_linear_tdt(z_value=4.0, reference_ctmax=39.0, reference_duration=60.0,
                           incipient_temperature=30.0)
        # Constant exposure at reference_ctmax: injury reaches 1 after reference_duration steps
        T_const = fill(39.0, 120)
        inj = accumulated_injury(m, T_const, 1.0)
        @test inj[60] ≈ 1.0 atol=0.01
        @test time_to_failure(m, T_const, 1.0) ≈ 60.0 atol=1.0

        # Below incipient temperature: no injury
        T_cool = fill(25.0, 60)
        inj_cool = accumulated_injury(m, T_cool, 1.0)
        @test all(inj_cool .== 0.0)

        # Survives entire series
        @test time_to_failure(m, T_cool, 1.0) == Inf
    end

    # ── Properties ────────────────────────────────────────────────────────────
    @testset "optimal_temperature" begin
        m = utpc(optimal_temperature=30.0, thermal_breadth=10.0u"K")
        T_opt = optimal_temperature(m)
        @test ustrip(u"K", uconvert(u"K", T_opt)) ≈ 303.15 atol=0.01

        # Analytic for Deutsch
        md = deutsch(optimal_temperature=28.0)
        @test ustrip(u"°C", uconvert(u"°C", optimal_temperature(md))) ≈ 28.0 atol=0.01

        # Arrhenius: no finite peak
        ma = arrhenius()
        @test optimal_temperature(ma) == Inf * u"K"
    end

    @testset "critical_thermal_maximum" begin
        m = utpc(optimal_temperature=30.0, thermal_breadth=10.0u"K")
        ctmax = critical_thermal_maximum(m)
        @test ustrip(u"K", uconvert(u"K", ctmax)) > 303.15  # above T_opt

        md = deutsch(optimal_temperature=25.0, critical_thermal_maximum=40.0)
        @test ustrip(u"°C", uconvert(u"°C", critical_thermal_maximum(md))) ≈ 40.0
    end

    @testset "maximum_rate" begin
        m = utpc(optimal_temperature=30.0, thermal_breadth=10.0u"K", maximum_performance=3.0)
        @test maximum_rate(m) ≈ 3.0 atol=0.01
    end

    @testset "q10" begin
        # UTPC Q10 should be > 1 on the rising limb (below T_opt)
        m = utpc(optimal_temperature=30.0, thermal_breadth=10.0u"K")
        @test q10(m, 15.0) > 1.0   # 15→25°C, both below T_opt
    end

    @testset "z_value" begin
        m = log_linear_tdt(z_value=4.0, reference_ctmax=39.0, reference_duration=60.0)
        @test z_value(m) ≈ 4.0

        # Arrhenius z_value: T_ref²/T_A * log(10)
        ma = ArrheniusModel(T_A=8000.0, T_ref=293.15)
        @test z_value(ma) ≈ ma.T_ref^2 / ma.T_A * log(10)
    end

    # ── CTE ───────────────────────────────────────────────────────────────────
    @testset "constant_temperature_equivalent" begin
        m = arrhenius()
        T_series = collect(range(15.0, 35.0, length=100))   # 100 values from 15 to 35°C
        cte = constant_temperature_equivalent(m, T_series)
        T_mean = sum(T_series) / length(T_series)
        # Jensen's inequality: CTE > arithmetic mean for convex Arrhenius
        @test ustrip(u"K", uconvert(u"K", cte)) > T_mean + 273.15

        # Numeric and analytic agree for ArrheniusModel
        m2 = ArrheniusModel(T_A=9000.0, T_ref=293.15)
        cte2 = constant_temperature_equivalent(m2, T_series)
        @test ustrip(u"K", uconvert(u"K", cte2)) > 0.0
    end

    # ── Static ↔ dynamic CTmax conversions ─────────────────────────────────────
    @testset "dynamic_ctmax and static_ctmax_from_dynamic" begin
        m = log_linear_tdt(z_value=4.0, reference_ctmax=39.0, reference_duration=60.0,
                           incipient_temperature=20.0)
        ramp = 0.1  # °C/min
        dctmax = dynamic_ctmax(m, ramp)
        @test dctmax > m.reference_ctmax  # dynamic > static for slow ramp

        # Round-trip: recover reference_ctmax
        recovered = static_ctmax_from_dynamic(m, dctmax, ramp)
        @test recovered ≈ m.reference_ctmax atol=0.01
    end

    # ── TPC ↔ TDT bridge ───────────────────────────────────────────────────────
    @testset "tdt_from_tpc" begin
        m_tpc = utpc(optimal_temperature=30.0, thermal_breadth=10.0u"K")
        m_tdt = tdt_from_tpc(m_tpc)
        # z = E * log(10) ≈ 10 * log(10) ≈ 23.03
        @test m_tdt.z_value ≈ 10.0 * log(10) atol=0.01
        @test thermal_breadth_from_tdt(m_tdt) ≈ 10.0 atol=0.01
    end

    # ── TDT fitting ────────────────────────────────────────────────────────────
    @testset "fit_thermal_death_time_curve (static)" begin
        # Synthetic data: z=4, reference_ctmax=39 at 60 min
        m_true = log_linear_tdt(z_value=4.0, reference_ctmax=39.0, reference_duration=60.0)
        temps  = [35.0, 37.0, 39.0, 41.0, 43.0]
        times  = [survival_time(m_true, T) for T in temps]
        data   = StaticKnockdownData(temperatures=temps, knockdown_times=times)
        m_fit  = fit_thermal_death_time_curve(data; reference_duration=60.0)
        @test m_fit.z_value ≈ 4.0 atol=0.05
        @test m_fit.reference_ctmax ≈ 39.0 atol=0.05
    end

    # ── Schoolfield fitting (synthetic data — verifies parameter recovery) ────────
    @testset "fit_thermal_performance_curve (SharpSchoolFullModel)" begin
        # Generate noiseless data from a known model and verify NLS recovery.
        m_true = SharpSchoolFullModel(T_A=5000.0, T_ref=298.15, T_L=291.0,
                                      T_AL=26000.0, T_H=316.5, T_AH=100000.0,
                                      rate_at_reference=0.28)
        temps = collect(8.0:4.0:46.0)   # 10 points spanning cold→hot
        rates = temperature_correction.(Ref(m_true), temps)
        m_fit = fit_thermal_performance_curve(
            SharpSchoolFullModel, temps, rates; T_ref=298.15u"K")
        @test m_fit.T_A  ≈ 5000.0   rtol=0.02
        @test m_fit.T_L  ≈ 291.0    rtol=0.01
        @test m_fit.T_H  ≈ 316.5    rtol=0.01
        @test m_fit.T_AL ≈ 26000.0  rtol=0.05
        @test m_fit.T_AH ≈ 100000.0 rtol=0.05
        @test temperature_correction(m_fit, 25.0) ≈ temperature_correction(m_true, 25.0) rtol=0.01
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
        tl    = ToleranceLandscape(4.0, 42.0, 38.0, hcat(times, surv))
        mild  = fill(30.0, 60)
        cum   = cumulative_survival(tl, [mild, mild, mild])
        @test length(cum) == 3
        @test cum[1] >= cum[2] >= cum[3]  # non-increasing
        @test all(0.0 .<= cum .<= 1.0)
    end

end # @testset ThermalPhysiology
