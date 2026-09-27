# Fluctuating temperatures

Body temperatures fluctuate, and because rates are nonlinear functions of temperature, the mean rate over a
fluctuating temperature is not the rate at the mean temperature. Where the rate curve is convex, as for the Arrhenius
model, fluctuations raise the mean rate (Jensen's inequality), and near and above the optimum of a performance curve,
where it is concave, they lower it (Ruel and Ayres 1999).

```@setup fluctuating
using Main.FigureHelpers
using CairoMakie, ThermalPhysiology, Unitful
```

## Mean rates

[`mean_correction_factor`](@ref) and [`mean_thermal_performance`](@ref) average a model over a series of temperatures.
For a daily cycle of body temperature of 20 ± 8 °C:

```@example fluctuating
body_temperature = [20.0 + 8.0 * sinpi(2h / 24) for h in 0:0.5:23.5] .* u"°C"
arrhenius = ArrheniusModel(0.65u"eV"; T_ref = 20.0u"°C")
mean_correction_factor(arrhenius, body_temperature), arrhenius(20.0u"°C")
```

The mean correction is higher than the correction at the mean temperature of 20 °C, which is 1.

## Constant temperature equivalent

The constant temperature equivalent (CTE) is the constant temperature that gives the same mean rate as the
fluctuating temperature series ``\{T_i\}``: the temperature ``T_{eq}`` for which

```math
f(T_{eq}) = \frac{1}{n} \sum_{i=1}^{n} f(T_i)
```

It is used to estimate DEB parameters from data collected at fluctuating temperatures, since DEBtool expects a constant
temperature (as in the `getCTE` function of NicheMapR). [`constant_temperature_equivalent`](@ref) has an
analytic solution for [`ArrheniusModel`](@ref), and solves the equation numerically for other models, including thermal
death time models, where the rate is that of injury, ``1/t(T)``:

::: tabs

== Arrhenius

```@example fluctuating
uconvert(u"°C", constant_temperature_equivalent(arrhenius, body_temperature))
```

== Sharpe-Schoolfield

```@example fluctuating
schoolfield = sharpe_schoolfield_deb(activation = 0.65u"eV", reference_temperature = 20.0u"°C",
                                     low_temperature = 5.0u"°C", low_deactivation = 3.0u"eV",
                                     high_temperature = 35.0u"°C", high_deactivation = 7.0u"eV")
uconvert(u"°C", constant_temperature_equivalent(schoolfield, body_temperature))
```

== Performance curve

```@example fluctuating
performance = gaussian(maximum_rate = 1.0, optimal_temperature = 32.0u"°C", width_parameter = 6.0u"K")
uconvert(u"°C", constant_temperature_equivalent(performance, body_temperature))
```

== Thermal death time

```@example fluctuating
tdt = log_linear_tdt(z_value = 4.0u"K", reference_ctmax = 39.0u"°C",
                     reference_duration = 60.0u"minute", incipient_temperature = 30.0u"°C")
hot_body_temperature = body_temperature .+ 12.0u"K"
constant_temperature_equivalent(tdt, hot_body_temperature)
```

:::

The CTE of the Arrhenius model is above the mean temperature, and increases with the amplitude of the fluctuation:

```@example fluctuating
amplitudes = 0.0:0.5:15.0
cte(model, amplitude) = ustrip(u"°C", constant_temperature_equivalent(model,
    [20.0 + amplitude * sinpi(2h / 24) for h in 0:0.5:23.5] .* u"°C"))
fig, ax = figure_axis("Amplitude of daily cycle (K)", "Constant temperature equivalent (°C)")
for T_A in (4000.0, 8000.0, 12000.0) .* u"K"
    model = ArrheniusModel(; T_A, T_ref = 20.0u"°C")
    lines!(ax, amplitudes, cte.(Ref(model), amplitudes); linewidth = 2, label = "T_A = $T_A")
end
hlines!(ax, [20.0]; color = :gray, linestyle = :dash)
axislegend(ax; position = :lt)
fig
```

The dashed line is the mean temperature. A higher Arrhenius temperature makes the curve more convex, and raises the
CTE further.

For a performance curve, the CTE is only defined where the mean performance is reached by the curve within the range
of temperatures; the search is bounded by the range of the series, widened by 5 K, and can be set with the `T_bounds`
keyword, in K.
