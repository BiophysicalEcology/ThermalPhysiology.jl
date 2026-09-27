# Get started

ThermalPhysiology.jl describes how the rates of biological processes, and the time organisms survive, depend on
temperature. Inputs and outputs are [Unitful.jl](https://github.com/PainterQubits/Unitful.jl) quantities.

```julia
using Pkg
Pkg.add("ThermalPhysiology")
```

## Three families of models

There are three families of models, each with its own function:

| Family | Function | Returns |
|:-------|:---------|:--------|
| [Arrhenius](manual/arrhenius.md) | [`temperature_correction`](@ref) | a correction factor, 1 at the reference temperature, or a corrected rate |
| [Thermal performance curves](manual/performance_curves.md) | [`thermal_performance`](@ref) | an absolute or relative rate |
| [Thermal death time](manual/thermal_death.md) | [`survival_time`](@ref) | the median time to knockdown at a constant temperature |

A model is a struct holding its parameters, built with a keyword constructor. Every model can also be called like a
function of temperature:

::: tabs

== Arrhenius

```@example get_started
using ThermalPhysiology, Unitful

arrhenius = ArrheniusModel(0.65u"eV"; T_ref = 20.0u"°C")
temperature_correction(arrhenius, 30.0u"°C")
```

== Performance curve

```@example get_started
tpc = utpc(optimal_temperature = 30.0u"°C", thermal_breadth = 10.0u"K", maximum_performance = 1.0)
thermal_performance(tpc, 25.0u"°C")
```

== Thermal death time

```@example get_started
tdt = log_linear_tdt(z_value = 4.0u"K", reference_ctmax = 39.0u"°C",
                     reference_duration = 60.0u"minute", incipient_temperature = 30.0u"°C")
survival_time(tdt, 41.0u"°C")
```

== Called as a function

```@example get_started
arrhenius(30.0u"°C"), tpc(25.0u"°C"), tdt(41.0u"°C")
```

:::

## Units

Temperatures can be given in any unit, and are converted as needed:

```@example get_started
tpc(298.15u"K") == tpc(25.0u"°C")
```

Bare numbers are rejected rather than assumed to be in some unit:

```@example get_started
try
    tpc(25.0)
catch err
    err
end
```

## Properties of a curve

The optimum, the critical thermal limits and the breadth of a performance curve are found analytically where
possible, and numerically otherwise:

```@example get_started
optimal_temperature(tpc), critical_thermal_maximum(tpc)
```

```@example get_started
q10(tpc, 20.0u"°C")
```

## Plotting

The models broadcast over vectors of temperatures, so they are easy to plot, here with
[CairoMakie](https://docs.makie.org):

```@example get_started
using CairoMakie

temperatures = collect(0.0:0.25:40.0) .* u"°C"
fig = Figure(size = (700, 400))
ax = Axis(fig[1, 1]; xlabel = "Temperature (°C)", ylabel = "Relative performance")
lines!(ax, ustrip.(temperatures), tpc.(temperatures); linewidth = 2)
fig
```

## Thermal death time

A thermal death time model gives the survival time at a constant temperature, the temperature that is lethal after
a given exposure, and the injury accumulated under a fluctuating temperature:

```@example get_started
lethal_temperature(tdt, 120.0u"minute")
```

```@example get_started
heatwave = [30.0 + 12.0 * sinpi(h / 24) for h in 0:0.1:24] .* u"°C"
time_to_failure(tdt, heatwave, 6.0u"minute")
```

## Fluctuating temperatures

Because rates are nonlinear in temperature, the mean rate under a fluctuating temperature is not the rate at the mean
temperature. The [`constant_temperature_equivalent`](@ref) (CTE) is the constant temperature with the same mean rate. For
example, here the CTE for a fluctuating temperature regime with a mean of 20 °C is: 

```@example get_started
body_temperatures = [20.0 + 8.0 * sinpi(2h / 24) for h in 0:23] .* u"°C"
uconvert(u"°C", constant_temperature_equivalent(arrhenius, body_temperatures))
```

## Fitting

Models can be fitted to data, for example a performance curve to measured rates:

```@example get_started
temperatures = [5, 10, 15, 20, 25, 30, 35, 40] .* u"°C"
rates = [0.05, 0.2, 0.4, 0.62, 0.85, 1.0, 0.72, 0.0]
fitted = fit_thermal_performance_curve(UniversalTPCModel, temperatures, rates)
```

See [Fitting to data](manual/fitting.md) for thermal death time data and tolerance landscapes.
