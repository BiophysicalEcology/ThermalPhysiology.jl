# API Reference

## Abstract types

```@docs
AbstractTPCModel
AbstractArrheniusModel
AbstractPhenomenologicalModel
AbstractTDTModel
```

## Primary dispatch functions

```@docs
thermal_performance
temperature_correction
survival_time
```

## Arrhenius-family structs

```@docs
ArrheniusModel
SharpSchoolHighModel
SharpSchoolLowModel
SharpSchoolFullModel
SharpSchoolDEBModel
JohnsonLewinModel
```

## Phenomenological TPC structs

```@docs
UniversalTPCModel
DeutschModel
Briere1Model
Briere2Model
GaussianModel
Thomas2012Model
Thomas2017Model
PawarModel
Lactin2Model
```

## TDT structs

```@docs
LogLinearTDTModel
ToleranceLandscape
```

## TPC properties

```@docs
optimal_temperature
critical_thermal_maximum
critical_thermal_minimum
thermal_breadth
maximum_rate
q10
```

## TDT properties

```@docs
lethal_temperature
median_lethal_temperature
z_value
thermal_death_slope
ctmax_at_duration
temperature_maximum
```

## TDT fluctuating-temperature functions

```@docs
accumulated_injury
time_to_failure
dynamic_ctmax
static_ctmax_from_dynamic
```

## ToleranceLandscape functions

```@docs
dynamic_survival
daily_mortality
cumulative_survival
```

## Constant temperature equivalent

```@docs
constant_temperature_equivalent
mean_correction_factor
mean_thermal_performance
```

## Cross-family conversions

```@docs
tdt_from_tpc
thermal_breadth_from_tdt
```

## Curve fitting

```@docs
fit_thermal_performance_curve
fit_thermal_death_time_curve
fit_tolerance_landscape
IndividualKnockdownData
StaticKnockdownData
DynamicKnockdownData
```

## Registry

```@docs
THERMAL_REGISTRY
model_names
thermal_models
```

## Unit utilities

```@docs
ea_to_ta
ta_to_ea
```
