# API

## Model functions

```@docs
temperature_correction
thermal_performance
survival_time
```

## Abstract types

```@docs
AbstractTPCModel
AbstractArrheniusModel
AbstractPhenomenologicalModel
AbstractTDTModel
AbstractRepairModel
```

## Arrhenius models

```@docs
ArrheniusModel
SharpSchoolHighModel
sharpe_schoolfield_high
SharpSchoolLowModel
sharpe_schoolfield_low
SharpSchoolFullModel
sharpe_schoolfield
SharpSchoolDEBModel
sharpe_schoolfield_deb
JohnsonLewinModel
johnson_lewin
```

## Thermal performance curves

```@docs
UniversalTPCModel
utpc
DeutschModel
deutsch
GaussianModel
gaussian
Briere1Model
briere
Briere2Model
Thomas2012Model
thomas
Thomas2017Model
PawarModel
pawar
Lactin2Model
lactin2
```

### Properties

```@docs
optimal_temperature
maximum_rate
critical_thermal_maximum
critical_thermal_minimum
thermal_breadth
q10
```

## Thermal death time

```@docs
LogLinearTDTModel
log_linear_tdt
ToleranceLandscape
tolerance_landscape
```

### Properties

```@docs
lethal_temperature
median_lethal_temperature
ctmax_at_duration
temperature_maximum
z_value
thermal_death_slope
dynamic_ctmax
static_ctmax_from_dynamic
```

### Injury at varying temperatures

```@docs
accumulated_injury
resettable_injury
time_to_failure
step_injury
NoRepair
FullRepairBelowThreshold
full_repair_below_threshold
repair_rate
resets_injury
```

### Tolerance landscapes

```@docs
dynamic_survival
daily_mortality
cumulative_survival
```

## Fluctuating temperatures

```@docs
constant_temperature_equivalent
mean_correction_factor
mean_thermal_performance
```

## Linking performance and death

```@docs
tdt_from_tpc
thermal_breadth_from_tdt
```

## Fitting

```@docs
fit_thermal_performance_curve
fit_thermal_death_time_curve
fit_tolerance_landscape
StaticKnockdownData
DynamicKnockdownData
BinarySurvivalData
IndividualKnockdownData
```

## Registry

```@docs
THERMAL_REGISTRY
model_names
thermal_models
```

## Units

```@docs
ea_to_ta
ta_to_ea
```
