# Metabolic rate in natural conditions

Organisms in nature rarely experience constant temperatures. This tutorial computes the metabolic rate of a lizard
in the microclimates of a real site, from [Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl), with
the metabolic model of the [previous tutorial](metabolic_rate.md).

```@setup natural
using Main.FigureHelpers
using CairoMakie
```

## The microclimate

Microclimate.jl computes hourly air temperatures above the ground and soil temperatures below it from weather, site
and soil data. Its example problem is for Madison, Wisconsin, USA, with one representative day in the middle of each
month:

```@example natural
using Microclimate, ThermalPhysiology, BiologicalScaling, Statistics, Unitful

problem = example_microclimate_problem()
micro = solve(problem)
problem.model.depths
```

The output has 24 hours for each of the 12 days. The air temperature is at the heights `problem.model.heights`,
here 1 cm and 2 m, and the soil temperature at each depth. The temperatures available to a lizard are those near the
ground surface, and those in a burrow:

```@example natural
depths = problem.model.depths
microhabitats = (
    "Air, 1 cm" => micro.profile.air_temperature[:, 1],
    "Soil, 10 cm" => micro.soil_temperature[:, findfirst(==(10.0u"cm"), depths)],
    "Soil, 30 cm" => micro.soil_temperature[:, findfirst(==(30.0u"cm"), depths)],
)
day(x, d) = x[(24 * (d - 1) + 1):(24 * d)]
months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
nothing # hide
```

## Metabolic rate over a day

The metabolic model is the Sharpe-Schoolfield model of the previous tutorial, for a 50 g lizard, with the standard
metabolic rate at 20 °C from BiologicalScaling.jl and the ``Q_{10}`` of Andrews and Pough (1985):

```@example natural
lizard_rate(T) = standard_metabolic_rate(Squamate(), 50.0u"g", T)
T_A = log(lizard_rate(30.0u"°C") / lizard_rate(20.0u"°C")) / (1 / 293.15u"K" - 1 / 303.15u"K")
metabolism = sharpe_schoolfield_high(activation = T_A, reference_temperature = 20.0u"°C",
                                     high_temperature = 40.0u"°C", high_deactivation = 7.0u"eV",
                                     rate_at_reference = lizard_rate(20.0u"°C"))
nothing # hide
```

On the July day, the temperature near the surface ranges widely, and the metabolic rate with it, while in the burrow
both are nearly constant:

```@example natural
hours = 0:23
fig = Figure(size = (700, 600))
ax1 = Axis(fig[1, 1]; ylabel = "Temperature (°C)", title = "July")
ax2 = Axis(fig[2, 1]; xlabel = "Hour", ylabel = "Metabolic rate (mW)")
for (label, temperature) in microhabitats
    july = day(temperature, 7)
    lines!(ax1, hours, ustrip.(u"°C", july); linewidth = 2, label)
    lines!(ax2, hours, ustrip.(u"mW", metabolism.(july)); linewidth = 2)
end
axislegend(ax1; position = :lt)
linkxaxes!(ax1, ax2)
hidexdecorations!(ax1; grid = false)
fig
```

## Energy use over the year

The energy used in a day is the sum of the hourly rates. A lizard that stays in one place all day uses:

```@example natural
daily_energy(temperature, d) = uconvert(u"kJ", sum(metabolism.(day(temperature, d))) * 1u"hr")
fig, ax = figure_axis("Month", "Energy used (kJ/day)"; xticks = (1:12, months))
for (label, temperature) in microhabitats
    scatterlines!(ax, 1:12, [ustrip(u"kJ", daily_energy(temperature, d)) for d in 1:12]; linewidth = 2, label)
end
axislegend(ax; position = :lt)
fig
```

In summer a lizard in a deep burrow uses less than half the energy of one near the surface. In winter the deep soil is the warmest place, and the lizard uses least near the frozen surface.

## Fluctuations and the mean

Because the metabolic rate is a convex function of temperature, the mean rate over a fluctuating day is higher than
the rate at the mean temperature. The [`constant_temperature_equivalent`](@ref) is the constant temperature with the
same mean rate:

```@example natural
surface = last(first(microhabitats))
markdown_table(["Month", "Mean temperature", "Constant temperature equivalent", "Mean rate / rate at mean temperature"],
    [(months[d],
      uconvert(u"°C", mean(day(surface, d))),
      uconvert(u"°C", constant_temperature_equivalent(metabolism, day(surface, d);
          T_bounds = ustrip.(u"K", extrema(day(surface, d))))),
      mean(metabolism.(day(surface, d))) / metabolism(mean(day(surface, d))))
     for d in (1, 4, 7, 10)])
```

Near the surface, the daily fluctuation raises the mean metabolic rate 15 to 25% above the rate at the mean temperature,
so estimates of energy use from mean temperatures are too low. The search for the constant
temperature equivalent is bounded here by the range of each day's temperatures, with `T_bounds` in K: the metabolic rate
peaks near 40 °C, and above the peak another temperature has the same rate.

## Thermoregulation

Lizards do not stay in one place: they move between microhabitats to keep their body temperature near a preferred
temperature. As a simple example, suppose the lizard can reach any height or depth of the microclimate, and each hour
chooses the one closest to a preferred body temperature of 33 °C:

```@example natural
all_temperatures = hcat(micro.profile.air_temperature, micro.soil_temperature)
preferred = uconvert(u"K", 33.0u"°C")
body_temperature = [all_temperatures[h, argmin(abs.(all_temperatures[h, :] .- preferred))]
                    for h in axes(all_temperatures, 1)]
thermoregulating_energy = [ustrip(u"kJ", uconvert(u"kJ", sum(metabolism.(day(body_temperature, d))) * 1u"hr"))
                           for d in 1:12]
round.(thermoregulating_energy; digits = 2)
```

```@example natural
fig, ax = figure_axis("Hour", "Body temperature (°C)"; title = "Thermoregulating lizard")
for d in (4, 7, 10)
    lines!(ax, hours, ustrip.(u"°C", day(body_temperature, d)); linewidth = 2, label = months[d])
end
hlines!(ax, [33.0]; color = :gray, linestyle = :dash)
axislegend(ax; position = :lt)
fig
```

In July the lizard holds about 33 °C from morning to evening, and uses more energy than it would in any one place, while in April and October no microhabitat is warm enough. Here
the body temperature is taken to be the temperature of the chosen microhabitat. A real body temperature depends on
radiation, wind and body size, which [HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl) computes,
and the choice of microhabitat is modelled in
[BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl).
