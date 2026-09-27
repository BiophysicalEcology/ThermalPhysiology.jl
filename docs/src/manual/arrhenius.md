# Arrhenius models

The Arrhenius-family models give the factor by which a rate at a reference temperature ``T_{ref}`` is multiplied to
give the rate at temperature ``T``. They are used to correct rates for temperature in Dynamic Energy Budget (DEB)
models (Kooijman 2010), and describe the rise of rates with temperature from enzyme kinetics, with or without the
inactivation of enzymes at low and high temperatures (Sharpe and DeMichele 1977, Schoolfield et al. 1981).

```@setup arrhenius
using Main.FigureHelpers
using CairoMakie, ThermalPhysiology, Unitful
```

All temperatures in these models are absolute, in K, and the models are evaluated with
[`temperature_correction`](@ref), or by calling the model. The parameters are the Arrhenius temperature ``T_A``,
the reference temperature ``T_{ref}``, and for inactivation the temperatures ``T_L`` and ``T_H`` at which half the
enzyme is inactive at low and high temperatures, and the Arrhenius temperatures of the inactivation, ``T_{AL}`` and
``T_{AH}``. The Arrhenius temperature is an activation energy ``E_a`` divided by the Boltzmann constant
``k_B``, ``T_A = E_a / k_B``, and the named constructors accept either (see [Units](introduction.md#Units)).

Kooijman (2010) prefers the Arrhenius temperature to the activation energy, as calling the Arrhenius relation
"mechanistic" for whole organisms overstates the case.

## Arrhenius

The simplest model has a single parameter, ``T_A``:

```math
f(T) = \exp\left(\frac{T_A}{T_{ref}} - \frac{T_A}{T}\right)
```

The correction is 1 at ``T_{ref}``, and increases without limit with temperature. [`ArrheniusModel`](@ref) takes
either an Arrhenius temperature or an activation energy:

::: tabs

== Arrhenius temperature

```@example arrhenius
m = ArrheniusModel(T_A = 8000.0u"K", T_ref = 20.0u"°C")
m.([10.0, 20.0, 30.0] .* u"°C")
```

== Activation energy

```@example arrhenius
m = ArrheniusModel(0.65u"eV"; T_ref = 20.0u"°C")
m.([10.0, 20.0, 30.0] .* u"°C")
```

:::

The Arrhenius temperature sets how steeply rates rise, and so the ``Q_{10}``, the factor by which a rate increases
over 10 °C:

```@example arrhenius
temperatures = collect(0.0:0.5:40.0) .* u"°C"
fig, ax = figure_axis("Temperature (°C)", "Correction factor")
for T_A in (4000.0, 8000.0, 12000.0) .* u"K"
    m = ArrheniusModel(; T_A, T_ref = 20.0u"°C")
    lines!(ax, ustrip.(temperatures), m.(temperatures); linewidth = 2,
           label = "T_A = $(T_A), Q10 = $(round(m(30.0u"°C") / m(20.0u"°C"); digits = 2))")
end
axislegend(ax; position = :lt)
fig
```

## Inactivation at high and low temperatures

At high temperatures enzymes denature, and at low temperatures they may be inactivated too, so that rates fall away
from the Arrhenius relation at both ends. The Sharpe-Schoolfield models 
(Sharpe and DeMichele 1977, Schoolfield et al. 1981) describe these with extra terms in the denominator:

::: tabs

== High

[`SharpSchoolHighModel`](@ref), built with [`sharpe_schoolfield_high`](@ref), has inactivation at high temperatures:

```math
r(T) = r_{ref} \frac{\exp\left(T_A/T_{ref} - T_A/T\right)}{1 + \exp\left(T_{AH}/T_H - T_{AH}/T\right)}
```

```@example arrhenius
high = sharpe_schoolfield_high(activation = 0.65u"eV", reference_temperature = 20.0u"°C",
                               high_temperature = 35.0u"°C", high_deactivation = 7.0u"eV")
high.([20.0, 30.0, 40.0] .* u"°C")
```

== Low

[`SharpSchoolLowModel`](@ref), built with [`sharpe_schoolfield_low`](@ref), has inactivation at low temperatures:

```math
r(T) = r_{ref} \frac{\exp\left(T_A/T_{ref} - T_A/T\right)}{1 + \exp\left(T_{AL}/T - T_{AL}/T_L\right)}
```

```@example arrhenius
low = sharpe_schoolfield_low(activation = 0.65u"eV", reference_temperature = 20.0u"°C",
                             low_temperature = 5.0u"°C", low_deactivation = 3.0u"eV")
low.([0.0, 10.0, 20.0] .* u"°C")
```

== Full

[`SharpSchoolFullModel`](@ref), built with [`sharpe_schoolfield`](@ref), is the full model of Schoolfield et al. (1981,
eq. 4) with inactivation at both ends, and the ``T/T_{ref}`` factor from collision theory:

```math
r(T) = r_{ref} \frac{T}{T_{ref}} \frac{\exp\left(T_A/T_{ref} - T_A/T\right)}
       {1 + \exp\left(T_{AL}/T - T_{AL}/T_L\right) + \exp\left(T_{AH}/T_H - T_{AH}/T\right)}
```

Here ``r_{ref}`` is the rate at ``T_{ref}`` in the absence of inactivation, so ``T_{ref}`` should be in the range
where the rate follows the Arrhenius relation.

```@example arrhenius
full = sharpe_schoolfield(activation = 0.65u"eV", reference_temperature = 20.0u"°C",
                          low_temperature = 5.0u"°C", low_deactivation = 3.0u"eV",
                          high_temperature = 35.0u"°C", high_deactivation = 7.0u"eV")
full.([0.0, 20.0, 40.0] .* u"°C")
```

== DEB

[`SharpSchoolDEBModel`](@ref), built with [`sharpe_schoolfield_deb`](@ref), is normalised as in the `tempcorr`
function of DEBtool, without the ``T/T_{ref}`` factor, so that the rate at ``T_{ref}`` is exactly ``r_{ref}``:

```math
r(T) = r_{ref} \exp\left(\frac{T_A}{T_{ref}} - \frac{T_A}{T}\right)
       \frac{1 + \exp\left(T_{AL}/T_{ref} - T_{AL}/T_L\right) + \exp\left(T_{AH}/T_H - T_{AH}/T_{ref}\right)}
            {1 + \exp\left(T_{AL}/T - T_{AL}/T_L\right) + \exp\left(T_{AH}/T_H - T_{AH}/T\right)}
```

```@example arrhenius
deb = sharpe_schoolfield_deb(activation = 0.65u"eV", reference_temperature = 20.0u"°C",
                             low_temperature = 5.0u"°C", low_deactivation = 3.0u"eV",
                             high_temperature = 35.0u"°C", high_deactivation = 7.0u"eV")
deb.([0.0, 20.0, 40.0] .* u"°C")
```

== Johnson-Lewin

[`JohnsonLewinModel`](@ref), built with [`johnson_lewin`](@ref), is the original enzyme kinetics model of Johnson and
Lewin (1946), of the same form as the model with high-temperature inactivation:

```math
r(T) = r_{ref} \frac{\exp\left(T_A/T_{ref} - T_A/T\right)}{1 + \exp\left(T_{AH}/T_H - T_{AH}/T\right)}
```

```@example arrhenius
jl = johnson_lewin(activation = 0.65u"eV", reference_temperature = 20.0u"°C",
                   high_temperature = 35.0u"°C", high_deactivation = 7.0u"eV")
jl.([20.0, 30.0, 40.0] .* u"°C")
```

:::

With the same parameters, the models are:

```@example arrhenius
arrhenius = ArrheniusModel(0.65u"eV"; T_ref = 20.0u"°C")
temperatures = collect(-5.0:0.25:45.0) .* u"°C"
fig, ax = figure_axis("Temperature (°C)", "Correction factor"; limits = (nothing, (0, 3)))
for (model, label) in ((arrhenius, "Arrhenius"), (high, "High inactivation"), (low, "Low inactivation"),
                       (full, "Full"), (deb, "DEB"))
    lines!(ax, ustrip.(temperatures), model.(temperatures); linewidth = 2, label)
end
axislegend(ax; position = :lt)
fig
```

At ``T_{ref}`` = 20 °C the DEB model is exactly 1, while the full model is a little lower, as there is some
inactivation even at the reference temperature.

## Arrhenius plots

On an Arrhenius plot, the logarithm of the rate against the inverse of the absolute temperature, the Arrhenius model
is a straight line with slope ``-T_A``, and inactivation bends the line down at each end:

```@example arrhenius
temperatures = collect(0.0:0.25:42.0) .* u"°C"
inverse_T = 1000 ./ ustrip.(u"K", temperatures)
fig, ax = figure_axis("1000 / T (1/K)", "ln(correction factor)")
lines!(ax, inverse_T, log.(arrhenius.(temperatures)); linewidth = 2, label = "Arrhenius")
lines!(ax, inverse_T, log.(full.(temperatures)); linewidth = 2, label = "Full")
axislegend(ax; position = :lb)
fig
```

Schoolfield et al. (1981) used the slopes of the three regions of this plot to estimate the parameters, which is how
the fits of these models are started, see [Fitting to data](fitting.md#Sharpe-Schoolfield-models).

## Rates

Given a `rate_at_reference`, the models return the corrected rate rather than the correction factor, in the units of
the rate:

```@example arrhenius
ingestion = sharpe_schoolfield_deb(activation = 0.65u"eV", reference_temperature = 20.0u"°C",
                                   low_temperature = 5.0u"°C", low_deactivation = 3.0u"eV",
                                   high_temperature = 35.0u"°C", high_deactivation = 7.0u"eV",
                                   rate_at_reference = 12.0u"J/d")
ingestion(20.0u"°C"), ingestion(30.0u"°C")
```
