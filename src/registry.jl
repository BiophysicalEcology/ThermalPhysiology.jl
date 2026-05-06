# ── Model Registry ─────────────────────────────────────────────────────────────

struct ThermalModelEntry
    name::String
    model_type::Type
    family::Symbol          # :arrhenius, :phenomenological, :tdt
    description::String
    parameters::Vector{Symbol}
    reference::String
end

const THERMAL_REGISTRY = Dict{Symbol, ThermalModelEntry}(
    :arrhenius => ThermalModelEntry(
        "Arrhenius",
        ArrheniusModel,
        :arrhenius,
        "One-parameter temperature correction: exp(T_A/T_ref - T_A/T)",
        [:T_A, :T_ref],
        "Kooijman (2000) DEB Theory"
    ),
    :sharpe_schoolfield_high => ThermalModelEntry(
        "Sharpe-Schoolfield (high deactivation)",
        SharpSchoolHighModel,
        :arrhenius,
        "Arrhenius with high-temperature enzyme deactivation",
        [:T_A, :T_ref, :T_H, :T_AH, :rate_at_reference],
        "Schoolfield, Sharpe & Magnuson (1981) J Theor Biol 88:719"
    ),
    :sharpe_schoolfield_low => ThermalModelEntry(
        "Sharpe-Schoolfield (low suppression)",
        SharpSchoolLowModel,
        :arrhenius,
        "Arrhenius with low-temperature enzyme suppression",
        [:T_A, :T_ref, :T_L, :T_AL, :rate_at_reference],
        "Schoolfield, Sharpe & Magnuson (1981) J Theor Biol 88:719"
    ),
    :sharpe_schoolfield => ThermalModelEntry(
        "Sharpe-Schoolfield (full)",
        SharpSchoolFullModel,
        :arrhenius,
        "Schoolfield (1981): Arrhenius ÷ (1 + low_term + high_term); rate_at_reference assumes T_ref in Arrhenius zone",
        [:T_A, :T_ref, :T_L, :T_AL, :T_H, :T_AH, :rate_at_reference],
        "Schoolfield, Sharpe & Magnuson (1981) J Theor Biol 88:719"
    ),
    :sharpe_schoolfield_deb => ThermalModelEntry(
        "Sharpe-Schoolfield (DEBtool)",
        SharpSchoolDEBModel,
        :arrhenius,
        "DEBtool-normalised: numerator correction ensures k(T_ref) = rate_at_reference exactly; used in tempcorr",
        [:T_A, :T_ref, :T_L, :T_AL, :T_H, :T_AH, :rate_at_reference],
        "Kooijman (2010) DEB Theory §2.6; DEBtool_J tempcorr.m"
    ),
    :johnson_lewin => ThermalModelEntry(
        "Johnson-Lewin",
        JohnsonLewinModel,
        :arrhenius,
        "Original enzyme kinetics model with high-T deactivation",
        [:T_A, :T_ref, :T_H, :T_AH, :rate_at_reference],
        "Johnson & Lewin (1946) J Cell Comp Physiol 28:47"
    ),
    :utpc => ThermalModelEntry(
        "Universal TPC",
        UniversalTPCModel,
        :phenomenological,
        "y(x)=exp(x)*(1-x), x=(T-T_opt)/E; two-parameter, scale-invariant",
        [:T_opt, :E, :maximum_performance],
        "Arnoldi, Jackson, Peralta-Maraver & Payne (2025) PNAS"
    ),
    :deutsch => ThermalModelEntry(
        "Deutsch",
        DeutschModel,
        :phenomenological,
        "Gaussian rise below T_opt, quadratic decline above T_opt to CTmax",
        [:maximum_rate, :optimal_temperature, :critical_thermal_maximum, :width_parameter],
        "Deutsch et al. (2008) PNAS 105:6668"
    ),
    :briere1 => ThermalModelEntry(
        "Briere 1",
        Briere1Model,
        :phenomenological,
        "rate = c*T*(T-T_min)*sqrt(T_max-T)",
        [:rate_constant, :minimum_temperature, :maximum_temperature],
        "Briere et al. (1999) Environ Entomol 28:22"
    ),
    :briere2 => ThermalModelEntry(
        "Briere 2",
        Briere2Model,
        :phenomenological,
        "Generalised Briere with shape exponent: rate = c*T*(T-T_min)*(T_max-T)^(1/d)",
        [:rate_constant, :minimum_temperature, :maximum_temperature, :shape_parameter],
        "Briere et al. (1999) Environ Entomol 28:22"
    ),
    :gaussian => ThermalModelEntry(
        "Gaussian",
        GaussianModel,
        :phenomenological,
        "Symmetric Gaussian: rate = P_max*exp(-0.5*((T-T_opt)/width)^2)",
        [:maximum_rate, :optimal_temperature, :width_parameter],
        "classic"
    ),
    :thomas2012 => ThermalModelEntry(
        "Thomas 2012",
        Thomas2012Model,
        :phenomenological,
        "Quadratic TPC: rate = c*(T-T_opt+b)*(b-(T-T_opt))",
        [:rate_constant, :shape_parameter, :optimal_temperature],
        "Thomas et al. (2012) Ecology 93:2483"
    ),
    :thomas2017 => ThermalModelEntry(
        "Thomas 2017",
        Thomas2017Model,
        :phenomenological,
        "Asymmetric skewed-Gaussian TPC",
        [:maximum_rate, :optimal_temperature, :width_parameter, :skewness],
        "Thomas et al. (2017) Global Change Biol 23:4566"
    ),
    :pawar => ThermalModelEntry(
        "Pawar",
        PawarModel,
        :phenomenological,
        "Metabolic TPC: Arrhenius rise with high-T deactivation, biologically parameterised",
        [:rate_at_reference, :activation_energy, :deactivation_energy, :peak_temperature, :reference_temperature],
        "Pawar et al. (2018) Funct Ecol 32:1874"
    ),
    :lactin2 => ThermalModelEntry(
        "Lactin 2",
        Lactin2Model,
        :phenomenological,
        "exp(ρ*T) − exp(ρ*T_max − (T_max−T)/δ) + λ",
        [:rate_constant, :maximum_temperature, :delta_temperature, :intercept],
        "Lactin et al. (1995) Environ Entomol 24:68"
    ),
    :log_linear_tdt => ThermalModelEntry(
        "Log-linear TDT",
        LogLinearTDTModel,
        :tdt,
        "Parametric TDT: t(T) = t_ref * 10^((CTmax - T)/z)",
        [:z_value, :reference_ctmax, :reference_duration, :incipient_temperature],
        "Jørgensen et al. (2021) Sci Reports 11:12962"
    ),
    :tolerance_landscape => ThermalModelEntry(
        "Tolerance Landscape",
        ToleranceLandscape,
        :tdt,
        "Semi-parametric TDT: full empirical S(τ) distribution with z-shifting",
        [:z_value, :temperature_maximum, :mean_assay_temperature, :survival_curve],
        "Rezende et al. (2014) Funct Ecol 28:799; Rezende et al. (2020) Science 369:1242"
    ),
)

"""
    model_names(; family=nothing) → Vector{Symbol}

Return all registered model keys, optionally filtered by `family`
(`:arrhenius`, `:phenomenological`, `:tdt`, or `nothing` for all).
"""
function model_names(; family=nothing)
    if isnothing(family)
        sort!(collect(keys(THERMAL_REGISTRY)))
    else
        sort!([k for (k, v) in THERMAL_REGISTRY if v.family == family])
    end
end

"""
    thermal_models(; family=nothing) → Dict{Symbol, ThermalModelEntry}

Return registry entries, optionally filtered by `family`.
"""
function thermal_models(; family=nothing)
    isnothing(family) ? copy(THERMAL_REGISTRY) :
        Dict(k => v for (k, v) in THERMAL_REGISTRY if v.family == family)
end
