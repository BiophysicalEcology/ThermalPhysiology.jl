# Fitting to data

Performance curves, thermal death time curves and tolerance landscapes can be fitted to experimental data, by
nonlinear least squares with [LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl), or by ordinary least squares
where the model is linear. The examples here fit data simulated from known models, with added noise.

```@setup fitting
using Main.FigureHelpers
using CairoMakie, Random, ThermalPhysiology, Unitful
Random.seed!(1234)
```

## Performance curves

[`fit_thermal_performance_curve`](@ref) takes the model type, the temperatures, and the rates, and returns a fitted
model. Initial parameters are estimated from the data unless `initial_parameters` are given, in the order of the fields
of the struct, and `weights` gives weighted least squares. Data are simulated from a universal performance curve:

```@example fitting
truth = utpc(optimal_temperature = 30.0u"°C", thermal_breadth = 6.0u"K", maximum_performance = 1.0)
temperatures = collect(10.0:2.0:36.0) .* u"°C"
rates = max.(0.0, truth.(temperatures) .+ 0.04 .* randn(length(temperatures)))
nothing # hide
```

Four models can be fitted this way:

::: tabs

== Universal

```@example fitting
fit_universal = fit_thermal_performance_curve(UniversalTPCModel, temperatures, rates)
parameter_table("Fitted" => fit_universal, "True" => truth)
```

== Gaussian

```@example fitting
fit_gaussian = fit_thermal_performance_curve(GaussianModel, temperatures, rates)
parameter_table(fit_gaussian)
```

== Deutsch

```@example fitting
fit_deutsch = fit_thermal_performance_curve(DeutschModel, temperatures, rates)
parameter_table(fit_deutsch)
```

== Brière

```@example fitting
fit_briere = fit_thermal_performance_curve(Briere1Model, temperatures, rates)
parameter_table(fit_briere)
```

:::

```@example fitting
curve_temperatures = collect(8.0:0.1:38.0) .* u"°C"
fig, ax = figure_axis("Temperature (°C)", "Performance"; limits = (nothing, (-0.05, 1.15)))
scatter!(ax, ustrip.(temperatures), rates; color = :black, label = "Data")
for (model, label) in ((fit_universal, "Universal"), (fit_gaussian, "Gaussian"),
                       (fit_deutsch, "Deutsch"), (fit_briere, "Brière 1"))
    lines!(ax, ustrip.(curve_temperatures), model.(curve_temperatures); linewidth = 2, label)
end
axislegend(ax; position = :lt)
fig
```

The fitted models are ordinary models, so their properties can be found as for any other:

```@example fitting
optimal_temperature(fit_universal), critical_thermal_maximum(fit_universal)
```

## Sharpe-Schoolfield models

[`SharpSchoolFullModel`](@ref) and [`SharpSchoolDEBModel`](@ref) have their own method of
[`fit_thermal_performance_curve`](@ref), with the reference temperature as the `T_ref` keyword. The initial parameters
are estimated by the graphical method of Schoolfield et al. (1981, fig. 2): the data are split into cold, middle and hot
thirds, and a line is fitted by least squares to the logarithm of the rate against ``1/T`` in each. From the slopes,

```math
T_A \approx -b_{middle}, \qquad T_{AL} \approx -b_{cold} - T_A, \qquad T_{AH} \approx b_{hot} + T_A
```

and ``T_L`` and ``T_H`` are where the cold and hot lines cross the middle line lowered by ``\ln 2``, where half the
enzyme is inactive. Data are simulated with multiplicative noise:

```@example fitting
schoolfield = sharpe_schoolfield(activation = 0.65u"eV", reference_temperature = 25.0u"°C",
                                 low_temperature = 10.0u"°C", low_deactivation = 3.0u"eV",
                                 high_temperature = 36.0u"°C", high_deactivation = 8.0u"eV",
                                 rate_at_reference = 1.0)
ss_temperatures = collect(2.0:2.0:42.0) .* u"°C"
ss_rates = schoolfield.(ss_temperatures) .* exp.(0.08 .* randn(length(ss_temperatures)))
nothing # hide
```

By default the fit minimises the squared differences of the logarithms of the rates, as Schoolfield et al. (1981) did,
which suits rates that vary over orders of magnitude. With `log_transform = false` it minimises the differences of the
rates themselves, optionally with `weights`:

::: tabs

== Log-transformed

```@example fitting
fit_log = fit_thermal_performance_curve(SharpSchoolFullModel, ss_temperatures, ss_rates; T_ref = 25.0u"°C")
parameter_table("Fitted" => fit_log, "True" => schoolfield)
```

== Absolute

```@example fitting
fit_absolute = fit_thermal_performance_curve(SharpSchoolFullModel, ss_temperatures, ss_rates;
                                             T_ref = 25.0u"°C", log_transform = false)
parameter_table("Fitted" => fit_absolute, "True" => schoolfield)
```

== Weighted

```@example fitting
fit_weighted = fit_thermal_performance_curve(SharpSchoolFullModel, ss_temperatures, ss_rates;
                                             T_ref = 25.0u"°C", log_transform = false,
                                             weights = 1 ./ ss_rates)
parameter_table("Fitted" => fit_weighted, "True" => schoolfield)
```

:::

```@example fitting
curve_temperatures = collect(0.0:0.1:44.0) .* u"°C"
fig, ax = figure_axis("Temperature (°C)", "Rate")
scatter!(ax, ustrip.(ss_temperatures), ss_rates; color = :black, label = "Data")
lines!(ax, ustrip.(curve_temperatures), schoolfield.(curve_temperatures); color = :gray, linestyle = :dash,
       linewidth = 2, label = "True")
for (model, label) in ((fit_log, "Log-transformed"), (fit_absolute, "Absolute"), (fit_weighted, "Weighted"))
    lines!(ax, ustrip.(curve_temperatures), model.(curve_temperatures); linewidth = 2, label)
end
axislegend(ax; position = :lt)
fig
```

## Thermal death time curves

[`fit_thermal_death_time_curve`](@ref) fits a [`LogLinearTDTModel`](@ref) to data from static assays at constant
temperatures, from ramping assays, or from survival scored at fixed exposures, and
[`fit_tolerance_landscape`](@ref) builds a [`ToleranceLandscape`](@ref) from individual knockdown times. Each kind of
data has its own container, with Unitful temperatures and times. The data are simulated from a known model:

```@example fitting
tdt = log_linear_tdt(z_value = 3.5u"K", reference_ctmax = 40.0u"°C",
                     reference_duration = 60.0u"minute", incipient_temperature = 32.0u"°C")
nothing # hide
```

::: tabs

== Static

[`StaticKnockdownData`](@ref) holds the mean or median knockdown time at each assay temperature. A line is fitted
by least squares to the logarithm of the time against temperature:

```@example fitting
static_temperatures = collect(38.0:1.0:44.0) .* u"°C"
static_times = survival_time.(Ref(tdt), static_temperatures) .* exp.(0.1 .* randn(length(static_temperatures)))
static_data = StaticKnockdownData(temperatures = static_temperatures, knockdown_times = static_times)
static_fit = fit_thermal_death_time_curve(static_data; reference_duration = 60.0u"minute")
parameter_table("Fitted" => static_fit, "True" => tdt)
```

== Dynamic

[`DynamicKnockdownData`](@ref) holds the knockdown temperatures of ramping assays at several ramp rates, from a start
temperature. With three or more ramp rates, eq. 7a of Jørgensen et al. (2021) is fitted by nonlinear least squares:

```@example fitting
ramp_rates = [0.05, 0.1, 0.25, 0.5, 1.0] .* u"K/minute"
dynamic_ctmaxes = dynamic_ctmax.(Ref(tdt), ramp_rates) .+ 0.1u"K" .* randn(length(ramp_rates))
dynamic_data = DynamicKnockdownData(ramp_rates = ramp_rates, dynamic_ctmax_values = dynamic_ctmaxes,
                                    start_temperature = 32.0u"°C")
dynamic_fit = fit_thermal_death_time_curve(dynamic_data; reference_duration = 60.0u"minute")
parameter_table("Fitted" => dynamic_fit, "True" => tdt)
```

== Binary survival

[`BinarySurvivalData`](@ref) holds whether each individual survived an exposure of some duration at some temperature.
The probability of surviving is fitted as a log-logistic function of the exposure time relative to the survival time,
in one step, without first summarising the data by group:

```@example fitting
exposure_temperatures = repeat(collect(38.0:1.0:43.0) .* u"°C", inner = 5 * 10)
exposure_times = repeat(repeat([5.0, 15.0, 30.0, 60.0, 120.0] .* u"minute", inner = 10), 6)
p_alive = @. 1 / (1 + (exposure_times / survival_time($(Ref(tdt)), exposure_temperatures))^3)
binary_data = BinarySurvivalData(temperatures = exposure_temperatures, exposure_times = exposure_times,
                                 survived = rand(length(p_alive)) .< p_alive)
binary_fit = fit_thermal_death_time_curve(binary_data; reference_duration = 60.0u"minute")
parameter_table("Fitted" => binary_fit, "True" => tdt)
```

== Tolerance landscape

[`IndividualKnockdownData`](@ref) holds the knockdown time of each individual. [`fit_tolerance_landscape`](@ref) fits
a line to the logarithm of the knockdown times against temperature for the z-value and ``T_{max}``, shifts all times
to the mean assay temperature with the z-value, and stores their empirical survival curve on `n_bins` points
(Rezende et al. 2014):

```@example fitting
individual_temperatures = repeat(collect(38.0:1.0:43.0) .* u"°C", inner = 20)
individual_times = survival_time.(Ref(tdt), individual_temperatures) .* exp.(0.4 .* randn(length(individual_temperatures)))
individual_data = IndividualKnockdownData(temperatures = individual_temperatures,
                                          knockdown_times = individual_times)
landscape = fit_tolerance_landscape(individual_data; n_bins = 500)
parameter_table(landscape)
```

:::

The true model had a z-value of 3.5 K and a static ``CT_{max}`` of 40 °C at 60 minutes, and a ``T_{max}`` of:

```@example fitting
temperature_maximum(tdt)
```
