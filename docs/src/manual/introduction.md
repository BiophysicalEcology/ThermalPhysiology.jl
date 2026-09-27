# Introduction

ThermalPhysiology.jl provides models of the effect of temperature on the rates of biological processes (thermal
performance curves and temperature corrections), and on survival at high temperatures (thermal death time). It is
built for the [BiophysicalEcology](https://github.com/BiophysicalEcology) ecosystem for mechanistic niche modelling,
for temperature corrections in Dynamic Energy Budget (DEB) theory, and for general use in studies of thermal responses.

```@setup intro
using ThermalPhysiology, Unitful
```

## A standalone package

ThermalPhysiology.jl is a complete package in its own right, with no dependence on the rest of the ecosystem. It can
be used wherever temperature affects rates or survival: to fit performance curves and thermal death time curves to
experimental data, to compare models, to temperature-correct physiological rates, or to estimate performance, injury
and survival from measured or modelled body temperatures, such as those from data loggers or other models. The
[Get started](../get_started.md) page and the rest of this manual use it on its own.

## ThermalPhysiology.jl in NicheMapper.jl

ThermalPhysiology.jl is also one of the packages of the [BiophysicalEcology](https://github.com/BiophysicalEcology)
organisation for mechanistic niche modelling, the Julia successor to [NicheMapR](https://github.com/mrke/NicheMapR).
These will be brought together under the umbrella package
[NicheMapper.jl](https://github.com/BiophysicalEcology/NicheMapper.jl), which, like AnimalMapper.jl, is in development.
Each package can be used on its own. A mechanistic niche model connects habitat features and weather to the metabolism of 
the organism, via the environment that the organism constructs through its morphology, physiology and behavior:

| Package | Role |
|:--------|:-----|
| [Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl) | Hourly microclimates above and below ground (radiation, wind, air and soil temperature and moisture, dew, frost, snow and pooling water) from site, weather and soil inputs |
| [MicroclimateMapper.jl](https://github.com/BiophysicalEcology/MicroclimateMapper.jl) | Runs Microclimate.jl over maps, with gridded data from RasterDataSources.jl and Rasters.jl |
| [HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl) | The heat budget of an organism (radiation, convection, conduction and evaporation) and its associated body temperature in a given microclimate |
| [BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl) | The choice of microhabitat, trait values and activity, and the resulting body temperatures, water loss rates, metabolic processes and activity times |
| [BiologicalScaling.jl](https://github.com/BiophysicalEcology/BiologicalScaling.jl) | Allometric scaling of metabolic rate, morphology, locomotion and life history with body size, used for the rates that ThermalPhysiology.jl corrects for temperature |
| ThermalPhysiology.jl | The consequences of body temperature: rates of performance, temperature corrections of physiological rates, and heat injury and survival |
| [AnimalMapper.jl](https://github.com/BiophysicalEcology/AnimalMapper.jl) (in development) | Animal niche models, running BiophysicalBehaviour.jl over time and space with the energy and water budgets of the animal, including Dynamic Energy Budget (DEB) models |

Dynamic Energy Budget models correct all their rates for body temperature.
[DEBtool_J.jl](https://github.com/add-my-pet/DEBtool_J.jl), the Julia version of the DEBtool software of the
[Add-my-Pet](https://www.bio.vu.nl/thb/deb/deblab/add_my_pet/) collection, will use ThermalPhysiology.jl for its
temperature corrections. [`ArrheniusModel`](@ref) and
[`SharpSchoolDEBModel`](@ref) follow the `tempcorr` convention of DEBtool, so parameters estimated for a species in
Add-my-Pet can be used directly. For DEB parameter estimation from data collected at fluctuating temperatures, the
[constant temperature equivalent](fluctuating.md) converts a body temperature series to the constant temperature DEBtool
expects for 'zero-variate' data.

The body temperatures from HeatExchange.jl and BiophysicalBehaviour.jl, and air/soil temperatures of Microclimate.jl 
and MicroclimateMapper.jl, can be used as input to the models of this package.
From a body temperature time series, for example hourly over a day, ThermalPhysiology.jl gives:

- the rate of any process, from a [thermal performance curve](performance_curves.md), and the time available for
  activity at a given level of performance,
- the [temperature correction](arrhenius.md) of physiological rates, such as those of a DEB model, and its mean
  over a day, or the [constant temperature equivalent](fluctuating.md) of the fluctuating body temperature,
- the [heat injury](thermal_death.md) accumulated when the body temperature exceeds the thermal limits, and the
  survival of heatwaves.

For example, with an hourly body temperature from a thermoregulating ectotherm (here made up as a sine wave in place of
the output of BiophysicalBehaviour.jl):

```@example intro
body_temperature = [24.0 + 14.0 * max(0.0, sinpi((h - 6) / 12)) for h in 0:23] .* u"°C"

performance = utpc(optimal_temperature = 33.0u"°C", thermal_breadth = 6.0u"K", maximum_performance = 1.0)
correction = ArrheniusModel(0.65u"eV"; T_ref = 20.0u"°C")
tdt = log_linear_tdt(z_value = 3.5u"K", reference_ctmax = 40.0u"°C",
                     reference_duration = 60.0u"minute", incipient_temperature = 36.0u"°C")

(
    mean_performance = mean_thermal_performance(performance, body_temperature),
    mean_correction = mean_correction_factor(correction, body_temperature),
    injury = last(accumulated_injury(tdt, body_temperature, 60.0u"minute")),
)
```

In AnimalMapper.jl, these are evaluated for every hour and location of a simulation. The
[tutorials](../tutorials/metabolic_rate.md) combine ThermalPhysiology.jl with BiologicalScaling.jl for metabolic rates,
and with Microclimate.jl for natural thermal conditions.

## Model families

There are three families of models, each with its own function:

| Family | Function | Output | Parameters in |
|:-------|:---------|:-------|:--------------|
| [Arrhenius](arrhenius.md) | [`temperature_correction`](@ref) | Factor equal to 1 at ``T_{ref}``, or a corrected rate | K |
| [Phenomenological](performance_curves.md) | [`thermal_performance`](@ref) | Absolute or relative rate | °C or K |
| [Thermal death time](thermal_death.md) | [`survival_time`](@ref) | Median time to knockdown | °C or K, and time |

The families share a mathematical core: the Arrhenius rate law underlies both temperature corrections and the rate of
heat damage, so the breadth of a performance curve and the slope of a thermal death time curve are related. See
[Linking performance and death](tpc_tdt_bridge.md).

## Models as structs

Each model is a struct containing its parameters. The same function name dispatches to the formula of each model
type, and every model can also be called directly:

```@example intro
m = utpc(optimal_temperature = 30.0u"°C", thermal_breadth = 10.0u"K", maximum_performance = 1.0)
m(25.0u"°C") == thermal_performance(m, 25.0u"°C")
```

Models can be built from the struct's own keyword constructor, with the short field names of the literature, or from
a named constructor with descriptive keywords:

::: tabs

== Named constructor

```@example intro
sharpe_schoolfield_high(activation = 0.65u"eV", reference_temperature = 20.0u"°C",
                        high_temperature = 45.0u"°C", high_deactivation = 7.0u"eV")
```

== Struct constructor

```@example intro
SharpSchoolHighModel(T_A = 7543.0u"K", T_ref = 293.15u"K", T_H = 318.15u"K", T_AH = 81232.0u"K")
```

:::

No parameter has a default value: every parameter is biologically meaningful, and must be chosen deliberately. The
docstring of each struct in the [API](../api.md) gives an example.

## Units

Every temperature, time, energy and rate is a [Unitful.jl](https://github.com/PainterQubits/Unitful.jl)
quantity, both in the model parameters and in the temperatures a model is evaluated at, including vectors of
temperatures. Bare numbers are rejected with an `ArgumentError`, rather than silently assumed to be in some unit.
Any compatible unit can be used:

::: tabs

== °C

```@example intro
arrhenius = ArrheniusModel(T_A = 8000.0u"K", T_ref = 20.0u"°C")
arrhenius(30.0u"°C")
```

== K

```@example intro
arrhenius(303.15u"K")
```

== °F

```@example intro
arrhenius(86.0u"°F")
```

:::

Some conventions follow from the units:

- Temperature differences, such as a thermal breadth or a z-value, are given in K. A value in °C is accepted and
  taken as a difference of the same size.
- Heating rates for thermal death time curves must be in K per unit of time, e.g. `0.1u"K/minute"`, since °C is an 
affine unit (arbitrary zero point) and cannot be used in a rate.
- Activation energies can be given as an energy, e.g. `0.65u"eV"`, or as an Arrhenius temperature
  ``T_A = E_a / k_B``, e.g. `7500.0u"K"`. [`ea_to_ta`](@ref) and [`ta_to_ea`](@ref) convert between them:

```@example intro
ea_to_ta(0.65u"eV"), ta_to_ea(7500.0u"K")
```

- A few curve-shape coefficients, such as the `rate_constant` and `intercept` of [`Lactin2Model`](@ref), are
  defined against the numeric value of the temperature in °C, and are bare numbers.
- Rate parameters, such as `maximum_rate` or `rate_at_reference`, can be bare numbers (a relative rate) or Unitful
  rates, and the units are kept in the output.

## Rates and corrections

The Arrhenius-family models other than [`ArrheniusModel`](@ref) have a `rate_at_reference` parameter. It is not a
stand-in value but a switch: without it [`temperature_correction`](@ref) returns the dimensionless correction factor,
and with it the corrected rate, in the units of the rate given:

::: tabs

== Correction factor

```@example intro
correction = sharpe_schoolfield_high(activation = 0.65u"eV", reference_temperature = 20.0u"°C",
                                     high_temperature = 45.0u"°C", high_deactivation = 7.0u"eV")
correction(30.0u"°C")
```

== Corrected rate

```@example intro
rate = sharpe_schoolfield_high(activation = 0.65u"eV", reference_temperature = 20.0u"°C",
                               high_temperature = 45.0u"°C", high_deactivation = 7.0u"eV",
                               rate_at_reference = 0.2u"d^-1")
rate(30.0u"°C")
```

:::

The switch is a type parameter of the struct, so both forms are type-stable.

## Model registry

All models are listed in [`THERMAL_REGISTRY`](@ref), with a description, the parameters and a reference.
[`model_names`](@ref) lists the keys, optionally by family:

::: tabs

== All

```@example intro
model_names()
```

== Arrhenius

```@example intro
model_names(family = :arrhenius)
```

== Phenomenological

```@example intro
model_names(family = :phenomenological)
```

== Thermal death time

```@example intro
model_names(family = :tdt)
```

:::

```@example intro
entry = THERMAL_REGISTRY[:utpc]
entry.description, entry.reference
```
