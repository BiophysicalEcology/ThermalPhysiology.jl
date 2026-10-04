# Thermal death time

At high temperatures, organisms survive for a limited time, and the time falls steeply as temperature rises. Thermal
death time (TDT) models describe this time, and the heat injury accumulated at varying temperatures. You can think
of it like a toxin - there's a threshold where it starts to act and a dosage that kills. These processes are
evaluated with [`survival_time`](@ref), the median time to knockdown (or death) at a constant temperature.

```@setup tdt
using Main.FigureHelpers
using CairoMakie, Random, ThermalPhysiology, Unitful
```

## Log-linear thermal death time

The logarithm of the knockdown time falls linearly with temperature (Rezende et al. 2014, Jørgensen et al. 2021).
[`LogLinearTDTModel`](@ref), built with [`log_linear_tdt`](@ref), gives the time at temperature ``T`` as

```math
t(T) = t_{ref}\, 10^{(CT_{max} - T) / z}
```

where ``CT_{max}`` is the static critical thermal maximum, the temperature at which knockdown takes the reference
duration ``t_{ref}``, and the z-value ``z`` is the rise in temperature that shortens the time ten-fold. Below the
incipient temperature, heat injury is negligible.

```@example tdt
tdt = log_linear_tdt(z_value = 4.0u"K", reference_ctmax = 39.0u"°C",
                     reference_duration = 60.0u"minute", incipient_temperature = 32.0u"°C")
survival_time.(Ref(tdt), [38.0, 40.0, 42.0] .* u"°C")
```

The z-value sets how quickly survival time falls with temperature. On a logarithmic axis, the curves are straight
lines crossing at ``CT_{max}`` and ``t_{ref}``:

```@example tdt
temperatures = collect(34.0:0.1:46.0) .* u"°C"
fig, ax = figure_axis("Temperature (°C)", "Survival time (minutes)"; yscale = log10)
for z in (2.0, 4.0, 6.0) .* u"K"
    m = log_linear_tdt(z_value = z, reference_ctmax = 39.0u"°C",
                       reference_duration = 60.0u"minute", incipient_temperature = 32.0u"°C")
    lines!(ax, ustrip.(u"°C", temperatures), ustrip.(u"minute", survival_time.(Ref(m), temperatures));
           linewidth = 2, label = "z = $z")
end
axislegend(ax; position = :rt)
fig
```

Times and temperatures can be in any units, which are kept:

```@example tdt
tdt_hours = log_linear_tdt(z_value = 4.0u"K", reference_ctmax = 39.0u"°C",
                           reference_duration = 1.0u"hr", incipient_temperature = 32.0u"°C")
survival_time(tdt_hours, 40.0u"°C"), uconvert(u"minute", survival_time(tdt_hours, 40.0u"°C"))
```

## Derived quantities

The critical thermal maximum depends on the duration of exposure. The model gives the static ``CT_{max}`` for any
duration, the lethal temperature for a duration, and the parameters in other forms:

::: tabs

== CTmax at a duration

[`ctmax_at_duration`](@ref) is the temperature at which knockdown takes the given time:

```@example tdt
ctmax_at_duration.(Ref(tdt), [1.0, 10.0, 60.0, 600.0] .* u"minute")
```

== Lethal temperature

[`lethal_temperature`](@ref), or [`median_lethal_temperature`](@ref), finds the same temperature numerically, from
[`survival_time`](@ref), so it works for any thermal death time model:

```@example tdt
lethal_temperature(tdt, 10.0u"minute")
```

== Maximum temperature

[`temperature_maximum`](@ref) is the ``T_{max}`` of Rezende et al. (2014), at which knockdown takes one minute:

```@example tdt
temperature_maximum(tdt)
```

== z-value and slope

[`z_value`](@ref) and [`thermal_death_slope`](@ref), the slope of the natural logarithm of time against temperature:

```@example tdt
z_value(tdt), thermal_death_slope(tdt)
```

:::

## Static and dynamic CTmax

The critical thermal maximum is often measured in a ramping assay, where temperature rises at a constant rate until
knockdown. The knockdown temperature of such an assay, the dynamic ``CT_{max}``, depends on the ramp rate, since heat
injury accumulates as the temperature rises. Jørgensen et al. (2021, eqs. 7a and 7b) give the dynamic ``CT_{max}`` from
the TDT parameters, and the static ``CT_{max}`` from a dynamic one:

::: tabs

== Dynamic from static

[`dynamic_ctmax`](@ref) predicts the knockdown temperature at a ramp rate:

```@example tdt
dynamic_ctmax.(Ref(tdt), [0.05, 0.1, 0.25, 1.0] .* u"K/minute")
```

== Static from dynamic

[`static_ctmax_from_dynamic`](@ref) recovers the static ``CT_{max}`` at the reference duration from a dynamic
``CT_{max}`` measured at a ramp rate:

```@example tdt
static_ctmax_from_dynamic(tdt, 42.9u"°C", 0.25u"K/minute")
```

:::

Faster ramps leave less time for injury to accumulate, so the dynamic ``CT_{max}`` rises with the ramp rate:

```@example tdt
ramp_rates = exp10.(range(-2, 1; length = 100)) .* u"K/minute"
fig, ax = figure_axis("Ramp rate (K/minute)", "Dynamic CTmax (°C)"; xscale = log10)
for z in (2.0, 4.0, 6.0) .* u"K"
    m = log_linear_tdt(z_value = z, reference_ctmax = 39.0u"°C",
                       reference_duration = 60.0u"minute", incipient_temperature = 32.0u"°C")
    lines!(ax, ustrip.(u"K/minute", ramp_rates), ustrip.(u"°C", dynamic_ctmax.(Ref(m), ramp_rates)); linewidth = 2,
           label = "z = $z")
end
axislegend(ax; position = :lt)
fig
```

The start temperature of the ramp is the incipient temperature unless the `start_temperature` keyword is given.

## Injury at varying temperatures

At a varying temperature, heat injury accumulates in each time step ``\Delta t`` as the fraction
``\Delta t / t(T)`` of the survival time at the temperature of that step, when the temperature is above the incipient
temperature (Jørgensen et al. 2021, eq. 3). An injury of 1 is the lethal dose. Whether injury is repaired at lower
temperatures is an open question, and there are two options:

::: tabs

== Accumulated

[`accumulated_injury`](@ref) never recovers, so injury only ever increases:

```@example tdt
day = [30.0 + 4.5 * max(0.0, sinpi((h - 6) / 12)) for h in 0:0.25:(24 * 3 - 0.25)] .* u"°C"
injury = accumulated_injury(tdt, day, 15.0u"minute")
maximum(injury)
```

== Resettable

[`resettable_injury`](@ref) is fully repaired whenever the temperature falls below the incipient temperature, unless
the lethal dose has been reached:

```@example tdt
reset_injury = resettable_injury(tdt, day, 15.0u"minute")
maximum(reset_injury)
```

:::

Over three warm days, with a daily maximum of 34.5 °C, injury accumulates each afternoon. Without repair the lethal dose
is reached on the third day, while with repair each day starts afresh and the lethal dose is never reached:

```@example tdt
hours = (0:length(day) - 1) ./ 4
fig = Figure(size = (700, 500))
ax1 = Axis(fig[1, 1]; ylabel = "Temperature (°C)")
lines!(ax1, hours, ustrip.(u"°C", day); linewidth = 2, color = :firebrick)
hlines!(ax1, [32.0]; color = :gray, linestyle = :dash)
ax2 = Axis(fig[2, 1]; xlabel = "Time (hours)", ylabel = "Injury")
lines!(ax2, hours, injury; linewidth = 2, label = "Accumulated")
lines!(ax2, hours, reset_injury; linewidth = 2, label = "Resettable")
axislegend(ax2; position = :lt)
linkxaxes!(ax1, ax2)
hidexdecorations!(ax1; grid = false)
fig
```

The dashed line is the incipient temperature. [`time_to_failure`](@ref) gives the time at which the lethal dose is
reached, or `Inf` if it is not:

```@example tdt
time_to_failure(tdt, day, 15.0u"minute"), time_to_failure(tdt, day, 15.0u"minute"; resettable = true)
```

### Repair models

Injury and repair are separate processes, combined in [`step_injury`](@ref), which updates the injury over one time
step. A repair model, a subtype of [`AbstractRepairModel`](@ref), defines the repair with [`repair_rate`](@ref), a
continuous reduction in each step, and [`resets_injury`](@ref), a complete reset. [`NoRepair`](@ref) is used by
[`accumulated_injury`](@ref), and [`FullRepairBelowThreshold`](@ref) by [`resettable_injury`](@ref), with the
incipient temperature as threshold.

[`step_injury`](@ref) is for simulations that keep injury as part of their own state, for example an hourly model of
development or behaviour, where the temperature of the next hour depends on the state of the organism:

```@example tdt
repair = full_repair_below_threshold(30.0u"°C")
injury_state = 0.0
for temperature in (35.0, 37.0, 38.0, 29.0, 37.0) .* u"°C"
    global injury_state = step_injury(tdt, repair, injury_state, temperature, 60.0u"minute")
    println(temperature, " => ", round(injury_state; digits = 3))
end
```

## Tolerance landscapes

The log-linear model describes the median knockdown time. Individuals vary, and the tolerance landscape of Rezende et
al. (2014, 2020), [`ToleranceLandscape`](@ref), keeps the whole distribution of knockdown times, the survival curve
``S(\tau)``, from individual knockdown data. The knockdown times at all assay temperatures are shifted with the z-value
to the mean assay temperature, and the survival curve at any other temperature is the same curve scaled in time.

A tolerance landscape is built from individual knockdown times with [`fit_tolerance_landscape`](@ref), here from
simulated data:

```@example tdt
Random.seed!(1)
assay_temperatures = repeat(36.0:1.0:42.0, inner = 20)
log_times = log10(60.0) .- (assay_temperatures .- 39.0) ./ 4.0 .+ 0.25 .* randn(length(assay_temperatures))
data = IndividualKnockdownData(temperatures = assay_temperatures .* u"°C",
                               knockdown_times = exp10.(log_times) .* u"minute")
landscape = fit_tolerance_landscape(data)
landscape.z_value, landscape.temperature_maximum, landscape.mean_assay_temperature
```

The median survival time at a temperature, in minutes, is given by [`survival_time`](@ref):

```julia
survival_time(landscape, 40.0u"°C")
```

### Dynamic survival

[`dynamic_survival`](@ref) follows the fraction of individuals surviving through a varying temperature, shifting the
survival curve in each time step to the current temperature (Rezende et al. 2020). In ramping assays at different
rates:

```@example tdt
fig, ax = figure_axis("Temperature (°C)", "Survival")
for ramp in (0.1, 0.25, 1.0)
    ramp_temperatures = collect(30.0:ramp:48.0) .* u"°C"
    survival = dynamic_survival(landscape, ramp_temperatures; dt_minutes = 1.0u"minute")
    lines!(ax, ustrip.(u"°C", ramp_temperatures), survival; linewidth = 2, label = "$ramp K/minute")
end
axislegend(ax; position = :lb)
fig
```

A `recovery_model` adds recovery in each time step, from the [`temperature_correction`](@ref) of an Arrhenius model or
the [`thermal_performance`](@ref) of a performance curve.

### Mortality over days

[`daily_mortality`](@ref) is the fraction dying over one day of temperatures, and [`cumulative_survival`](@ref) the
fraction surviving a series of days, with full recovery overnight (Rezende et al. 2020):

```julia
# one vector of temperatures, at one-minute steps, for each day
heatwave = [[30.0 + (7.0 + d) * max(0.0, sinpi((m / 60 - 6) / 12)) for m in 0:1439] .* u"°C" for d in 0:4]
daily_mortality.(Ref(landscape), heatwave)
cumulative_survival(landscape, heatwave)
```
