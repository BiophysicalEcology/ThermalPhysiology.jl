module ThermalPhysiology

# ── Includes (order matters: types → utils → models → properties → fitting → registry) ──

include("types.jl")
include("utils.jl")
include("arrhenius.jl")
include("phenomenological.jl")
include("thermal_death.jl")
include("properties.jl")
include("fitting.jl")
include("registry.jl")

# ── Exports ────────────────────────────────────────────────────────────────────

# Abstract types
export AbstractTPCModel, AbstractArrheniusModel, AbstractPhenomenologicalModel, AbstractTDTModel,
       AbstractRepairModel

# Arrhenius-family structs
export ArrheniusModel, SharpSchoolHighModel, SharpSchoolLowModel, SharpSchoolFullModel,
       SharpSchoolDEBModel, JohnsonLewinModel

# Phenomenological TPC structs
export UniversalTPCModel, DeutschModel, Briere1Model, Briere2Model,
       GaussianModel, Thomas2012Model, Thomas2017Model, PawarModel, Lactin2Model

# TDT structs + data input types
export LogLinearTDTModel, ToleranceLandscape
export IndividualKnockdownData, StaticKnockdownData, DynamicKnockdownData, BinarySurvivalData

# Primary dispatch functions
export thermal_performance, temperature_correction, survival_time

# Named constructors
export arrhenius, sharpe_schoolfield, sharpe_schoolfield_deb,
       sharpe_schoolfield_high, sharpe_schoolfield_low,
       johnson_lewin, utpc, deutsch, briere, gaussian, thomas, pawar, lactin2,
       log_linear_tdt, tolerance_landscape

# TPC properties
export optimal_temperature, critical_thermal_maximum, critical_thermal_minimum,
       thermal_breadth, maximum_rate, q10

# TDT properties
export lethal_temperature, median_lethal_temperature, z_value, thermal_death_slope,
       ctmax_at_duration, temperature_maximum

# TDT fluctuating-temperature functions (LogLinearTDTModel)
export accumulated_injury, resettable_injury, time_to_failure, dynamic_ctmax, static_ctmax_from_dynamic

# Repair models (damage recovery, evaluated in parallel with accrual)
export NoRepair, FullRepairBelowThreshold, full_repair_below_threshold,
       repair_rate, resets_injury, step_injury

# ToleranceLandscape functions (Rezende 2020)
export dynamic_survival, daily_mortality, cumulative_survival

# CTE and diagnostics
export constant_temperature_equivalent, mean_correction_factor, mean_thermal_performance

# Cross-family conversions
export tdt_from_tpc, thermal_breadth_from_tdt

# Fitting
export fit_thermal_performance_curve, fit_thermal_death_time_curve, fit_tolerance_landscape

# Registry
export THERMAL_REGISTRY, model_names, thermal_models

# Unit utilities
export ea_to_ta, ta_to_ea

end # module ThermalPhysiology
