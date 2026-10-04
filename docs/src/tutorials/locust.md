# Fitting to data: the Australian plague locust

This tutorial fits a thermal performance curve and a thermal death time curve for the Australian plague locust,
*Chortoicetes terminifera*, and uses them to predict the fate of eggs, laid in the soil near the surface, in a heatwave.

```@setup locust
using Main.FigureHelpers
using CairoMakie
```

## Egg development

Gregg (1981) measured egg development time at constant temperatures:

```@example locust
using ThermalPhysiology, CSV, DataFrames, Statistics, Unitful

eggs = CSV.read(joinpath(pkgdir(ThermalPhysiology), "docs", "src", "data", "chortoicetes_egg_development.csv"), DataFrame)
```

The development rate is the inverse of the development time in days:

```@example locust
temperatures = eggs.temperature .* u"°C"
development_rate = 1 ./ eggs.development_time
nothing # hide
```

Three curves fitted to the rates, the Sharpe-Schoolfield fit starting from the graphical estimates of Schoolfield et al.
(1981) (see [Fitting to data](../manual/fitting.md#Sharpe-Schoolfield-models)):

::: tabs

== Sharpe-Schoolfield

```@example locust
egg_development = fit_thermal_performance_curve(SharpSchoolFullModel, temperatures, development_rate;
                                                T_ref = 28.5u"°C")
parameter_table(egg_development)
```

== Universal

```@example locust
egg_universal = fit_thermal_performance_curve(UniversalTPCModel, temperatures, development_rate)
parameter_table(egg_universal)
```

== Gaussian

```@example locust
egg_gaussian = fit_thermal_performance_curve(GaussianModel, temperatures, development_rate)
parameter_table(egg_gaussian)
```

:::

```@example locust
curve_temperatures = collect(10.0:0.1:44.0) .* u"°C"
fig, ax = figure_axis("Temperature (°C)", "Development rate (1/day)"; limits = (nothing, (0, 0.11)))
scatter!(ax, eggs.temperature, development_rate; color = :black, label = "Gregg (1981)")
for (model, label) in ((egg_development, "Sharpe-Schoolfield"), (egg_universal, "Universal"),
                       (egg_gaussian, "Gaussian"))
    lines!(ax, ustrip.(u"°C", curve_temperatures), max.(0.0, model.(curve_temperatures)); linewidth = 2, label)
end
axislegend(ax; position = :lt)
fig
```

Only the Sharpe-Schoolfield model follows both the slow development in the cold and the fall at 39 °C. Its optimum:

```@example locust
grid = collect(10.0:0.01:45.0) .* u"°C"
grid_rates = egg_development.(grid)
optimum = grid[argmax(grid_rates)]
markdown_table(["Optimal temperature", "Maximum development rate (1/day)"], [(optimum, maximum(grid_rates))])
```

## Heat tolerance

First-instar nymphs were exposed to 44–55 °C for 30 seconds to 24 hours, and scored as alive at the end of the exposure
or 24 hours later (Kearney and Aitkenhead, in prep.). We assume eggs are as tolerant as nymphs. Unscored individuals are
left out:

```@example locust
tolerance = CSV.read(joinpath(pkgdir(ThermalPhysiology), "docs", "src", "data", "chortoicetes_heat_tolerance.csv"), DataFrame)
tolerance = filter(row -> !(ismissing(row.end_score) && ismissing(row.recovery_score)), tolerance)
tolerance.survived = coalesce.(tolerance.end_score, 0) .== 1 .|| coalesce.(tolerance.recovery_score, 0) .== 1
groups = combine(groupby(tolerance, [:exposure_time, :temperature]), :survived => mean => :surviving, nrow => :n)
first(groups, 10)
```

Rather than estimating a median lethal time for each group and then regressing, [`BinarySurvivalData`](@ref) fits
the individual outcomes in one step:

```@example locust
survival_data = BinarySurvivalData(temperatures = tolerance.temperature .* u"°C",
                                   exposure_times = tolerance.exposure_time .* u"minute",
                                   survived = tolerance.survived)
egg_tdt = fit_thermal_death_time_curve(survival_data; reference_duration = 60.0u"minute")
parameter_table(egg_tdt)
```

Every 2.9 K kills ten times faster. Critical temperatures for different exposures:

```@example locust
exposures = [1.0, 10.0, 60.0, 720.0, 1440.0] .* u"minute"
markdown_table(["Exposure time", "Critical temperature"], zip(exposures, ctmax_at_duration.(Ref(egg_tdt), exposures)))
```

```@example locust
tdt_temperatures = collect(43.0:0.1:56.0) .* u"°C"
fig, ax = figure_axis("Temperature (°C)", "Exposure time (minutes)"; yscale = log10)
scatter!(ax, groups.temperature, groups.exposure_time; color = groups.surviving, colormap = :RdYlBu,
         colorrange = (0, 1), markersize = 14)
lines!(ax, ustrip.(u"°C", tdt_temperatures), ustrip.(u"minute", survival_time.(Ref(egg_tdt), tdt_temperatures));
       color = :black, linewidth = 2)
Colorbar(fig[1, 2]; colormap = :RdYlBu, limits = (0, 1), label = "Fraction surviving")
fig
```

## The incipient temperature

We take heat injury to begin where development has fallen to 90% of its maximum above the optimum:

```@example locust
incipient = grid[findfirst(i -> grid[i] > optimum && grid_rates[i] < 0.9 * maximum(grid_rates), eachindex(grid))]
egg_tdt = log_linear_tdt(z_value = egg_tdt.z_value, reference_ctmax = egg_tdt.reference_ctmax,
                         reference_duration = egg_tdt.reference_duration, incipient_temperature = incipient)
incipient
```

## Survival during development

Gregg (1981) also recorded egg survival over development, exposures of 10 to 120 days. Relative to survival at
benign temperatures (20–33 °C), and with development time as a fraction of the median survival time from the assays:

```@example locust
benign = filter(row -> 20 <= row.temperature <= 33, eggs)
control_survival = mean(benign.survival) / 100
relative_survival = min.(1.0, eggs.survival ./ 100 ./ control_survival)
exposure_fraction = [uconvert(NoUnits, row.development_time * u"d" / survival_time(egg_tdt, row.temperature * u"°C"))
                     for row in eachrow(eggs)]
markdown_table(["Temperature", "Development time (days)", "Relative survival", "Development time / median survival time"],
               zip(temperatures, eggs.development_time, relative_survival, exposure_fraction))
```

Mortality begins, as expected, just above the incipient temperature: none at 36 °C, a third at 39 °C. But at 39 °C
development takes only 14% of the predicted median survival time, so the model predicts almost no deaths.

Fitting the thermal death time model to these data (assuming 100 eggs per temperature, as numbers were not recorded)
gives a much shallower curve than the assays, though poorly determined, as only 39 °C shows heat mortality:

```@example locust
hot = filter(row -> row.temperature >= 25, eggs)
eggs_per_temperature = 100
gregg_temperatures = Float64[]
gregg_times = Float64[]
gregg_survived = Bool[]
for row in eachrow(hot)
    alive = round(Int, min(1.0, row.survival / 100 / control_survival) * eggs_per_temperature)
    append!(gregg_temperatures, fill(row.temperature, eggs_per_temperature))
    append!(gregg_times, fill(row.development_time, eggs_per_temperature))
    append!(gregg_survived, [fill(true, alive); fill(false, eggs_per_temperature - alive)])
end
gregg_data = BinarySurvivalData(temperatures = gregg_temperatures .* u"°C", exposure_times = gregg_times .* u"d",
                                survived = gregg_survived)
gregg_tdt = fit_thermal_death_time_curve(gregg_data; reference_duration = 60.0u"minute",
                                         initial_reference_ctmax = 45.0)
parameter_table("Development survival" => gregg_tdt, "Heat tolerance assays" => egg_tdt)
```

```@example locust
fig, ax = figure_axis("Temperature (°C)", "Exposure time (minutes)"; yscale = log10)
lines!(ax, ustrip.(u"°C", tdt_temperatures), ustrip.(u"minute", survival_time.(Ref(egg_tdt), tdt_temperatures));
       color = :black, linewidth = 2, label = "Assays, median survival time")
scatter!(ax, groups.temperature, groups.exposure_time; color = groups.surviving, colormap = :RdYlBu,
         colorrange = (0, 1), markersize = 12, label = "Assays")
scatter!(ax, eggs.temperature, eggs.development_time .* 1440; color = relative_survival, colormap = :RdYlBu,
         colorrange = (0, 1), marker = :diamond, markersize = 16, strokewidth = 1, label = "Development (Gregg 1981)")
xlims!(ax, 30, 56)
Colorbar(fig[1, 2]; colormap = :RdYlBu, limits = (0, 1), label = "Fraction surviving")
axislegend(ax; position = :rt)
fig
```

The assay model suits exposures of hours, such as hot afternoons, but underestimates mortality over days above the
incipient temperature.

## A heatwave

Microclimate.jl builds hourly temperatures from daily minima and maxima with a `DielCurve`: a sine rise from sunrise
to an hour after midday, and an exponential decay after sunset. A week near the soil surface, with sunrise 7 hours
before solar noon:

```@example locust
using Microclimate

minimum_temperature = [18.0, 21.0, 24.0, 26.0, 27.0, 25.0, 20.0] .* u"°C"
maximum_temperature = [36.0, 40.0, 44.0, 46.0, 47.0, 44.0, 38.0] .* u"°C"
days = length(minimum_temperature)
curve = DielCurve((Sine(Sunrise(), Midday(1)), Decay(Sunset(), Sunrise())); inputs = (Sunrise(), Midday(1)))
forcing = DielForcing(curve, (minimum_temperature, maximum_temperature))
solar = (hour_solar_noon = fill(12.0, days), hour_angle_sunrise = fill(7.0, days))
egg_temperature = uconvert.(u"°C", Microclimate.evaluate!(fill(0.0u"K", 24 * days), forcing, solar, true))
nothing # hide
```

Heat injury accumulates on the hottest afternoons while development continues:

```@example locust
dt = 60.0u"minute"
injury = accumulated_injury(egg_tdt, egg_temperature, dt)
reset_injury = resettable_injury(egg_tdt, egg_temperature, dt)
development = cumsum(egg_development.(egg_temperature) .* (1 / 24))
hours = (0:(24 * days - 1)) ./ 24
fig = Figure(size = (700, 750))
ax1 = Axis(fig[1, 1]; ylabel = "Temperature (°C)")
lines!(ax1, hours, ustrip.(u"°C", egg_temperature); color = :firebrick, linewidth = 2)
hlines!(ax1, [ustrip(u"°C", incipient)]; color = :gray, linestyle = :dash)
ax2 = Axis(fig[2, 1]; ylabel = "Heat injury")
lines!(ax2, hours, injury; linewidth = 2, label = "Accumulated")
lines!(ax2, hours, reset_injury; linewidth = 2, label = "Resettable")
axislegend(ax2; position = :lt)
ax3 = Axis(fig[3, 1]; xlabel = "Day", ylabel = "Development")
lines!(ax3, hours, development; color = :black, linewidth = 2)
linkxaxes!(ax1, ax2, ax3)
hidexdecorations!(ax1; grid = false)
hidexdecorations!(ax2; grid = false)
fig
```

The dashed line is the incipient temperature. Survival depends on overnight repair:

```@example locust
markdown_table(["Repair", "Time to lethal dose", "Development completed"],
               [("None", uconvert(u"d", time_to_failure(egg_tdt, egg_temperature, dt)), last(development)),
                ("Overnight", uconvert(u"d", time_to_failure(egg_tdt, egg_temperature, dt; resettable = true)),
                 last(development))])
```

In practice, egg temperatures would come from Microclimate.jl's soil temperature at the laying depth, driven by the
site's weather.
