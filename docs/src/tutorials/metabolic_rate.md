# Metabolic rate and temperature

Metabolic rate depends on body size and on body temperature. Allometric equations describe the effect of size, and
are provided by [BiologicalScaling.jl](https://github.com/BiophysicalEcology/BiologicalScaling.jl), another package of
the BiophysicalEcology ecosystem. This tutorial combines them with the temperature models of ThermalPhysiology.jl.

```@setup metabolic
using Main.FigureHelpers
using CairoMakie
```

## Metabolic rate and body size

For squamate reptiles, BiologicalScaling.jl gives the standard metabolic rate (fasted and inactive) of Andrews and
Pough (1985) as a function of body mass and body temperature:

```@example metabolic
using BiologicalScaling, ThermalPhysiology, Unitful

standard_metabolic_rate(Squamate(), 50.0u"g", 20.0u"°C")
```

The rate scales with mass to the power 0.8, so on logarithmic axes it is a straight line:

```@example metabolic
masses = exp10.(range(0, 4; length = 50)) .* u"g"
fig, ax = figure_axis("Body mass (g)", "Standard metabolic rate (mW)"; xscale = log10, yscale = log10)
for temperature in [10.0, 20.0, 30.0] .* u"°C"
    rates = standard_metabolic_rate.(Ref(Squamate()), masses, temperature)
    lines!(ax, ustrip.(u"g", masses), ustrip.(u"mW", rates); linewidth = 2, label = "$temperature")
end
axislegend(ax; position = :lt)
fig
```

## The temperature dependence

In the equation of Andrews and Pough (1985) metabolic rate increases exponentially with temperature, by
``10^{0.038 T}``, so the ``Q_{10}`` is the same at all temperatures:

```@example metabolic
lizard_rate(T) = standard_metabolic_rate(Squamate(), 50.0u"g", T)
lizard_rate(30.0u"°C") / lizard_rate(20.0u"°C"), lizard_rate(40.0u"°C") / lizard_rate(30.0u"°C")
```

The Arrhenius model is the temperature dependence of metabolic rate in the metabolic theory of ecology (Gillooly et
al. 2001), and in DEB theory. It can be separated from the size dependence: the rate at a reference temperature is
taken from the allometry, and corrected for temperature with ThermalPhysiology.jl. The Arrhenius temperature that gives
the same ``Q_{10}`` of 2.4 between 20 and 30 °C is

```@example metabolic
Q10 = lizard_rate(30.0u"°C") / lizard_rate(20.0u"°C")
T_A = log(Q10) / (1 / 293.15u"K" - 1 / 303.15u"K")
```

or an activation energy of

```@example metabolic
ta_to_ea(T_A)
```

close to the 0.65 eV of the metabolic theory of ecology. With this, and the rate at 20 °C from the allometry:

```@example metabolic
reference_rate = lizard_rate(20.0u"°C")
arrhenius = ArrheniusModel(; T_A, T_ref = 20.0u"°C")
arrhenius_rate(T) = reference_rate * arrhenius(T)
arrhenius_rate(30.0u"°C"), lizard_rate(30.0u"°C")
```

The two agree at 20 and 30 °C, but the Arrhenius model has a ``Q_{10}`` that falls with temperature, while the
exponential model's does not.

## Inactivation at high temperatures

Neither model describes the fall in metabolic rate at high temperatures, as enzymes are inactivated. A
Sharpe-Schoolfield model adds this, here with half the enzyme inactive at 40 °C, and gives the rate in the units of the
reference rate, W:

```@example metabolic
schoolfield = sharpe_schoolfield_high(activation = T_A, reference_temperature = 20.0u"°C",
                                      high_temperature = 40.0u"°C", high_deactivation = 7.0u"eV",
                                      rate_at_reference = reference_rate)
schoolfield(30.0u"°C")
```

```@example metabolic
temperatures = collect(0.0:0.25:48.0) .* u"°C"
fig, ax = figure_axis("Body temperature (°C)", "Standard metabolic rate (mW)")
lines!(ax, ustrip.(u"°C", temperatures), ustrip.(u"mW", lizard_rate.(temperatures)); linewidth = 2,
       label = "Andrews and Pough (1985)")
lines!(ax, ustrip.(u"°C", temperatures), ustrip.(u"mW", arrhenius_rate.(temperatures)); linewidth = 2,
       label = "Arrhenius")
lines!(ax, ustrip.(u"°C", temperatures), ustrip.(u"mW", schoolfield.(temperatures)); linewidth = 2,
       label = "Sharpe-Schoolfield")
axislegend(ax; position = :lt)
fig
```

Below about 30 °C the three are nearly the same. They differ at the high body temperatures that an active lizard
reaches, where the choice of model matters most for energy budgets. The [next tutorial](microclimate.md) uses these
models in natural thermal conditions.

## Choosing a thermal response

The allometry gives the effect of size, and any thermal response in ThermalPhysiology.jl can be applied to it. The
reference rate is taken from the Andrews and Pough (1985) equation at a reference temperature, here for a 200 g lizard
at 25 °C:

```@example metabolic
reference_temperature = 25.0u"°C"
reference_rate = standard_metabolic_rate(Squamate(), 200.0u"g", reference_temperature)
```

and passed as the `rate_at_reference` of a model, which then returns the metabolic rate at any body temperature, in W:

::: tabs

== Arrhenius

The Arrhenius model with the activation energy of 0.65 eV of the metabolic theory of ecology. [`ArrheniusModel`](@ref)
gives the correction factor, which multiplies the reference rate:

```@example metabolic
mte = ArrheniusModel(0.65u"eV"; T_ref = reference_temperature)
mte_rate(T) = reference_rate * mte(T)
mte_rate(35.0u"°C")
```

== Sharpe-Schoolfield

The full Sharpe-Schoolfield model, with inactivation below 5 °C and above 40 °C:

```@example metabolic
full = sharpe_schoolfield(activation = 0.65u"eV", reference_temperature = reference_temperature,
                          low_temperature = 5.0u"°C", low_deactivation = 3.0u"eV",
                          high_temperature = 40.0u"°C", high_deactivation = 7.0u"eV",
                          rate_at_reference = reference_rate)
full(35.0u"°C")
```

== DEB

The Sharpe-Schoolfield model as normalised in DEBtool, which returns exactly the reference rate at the reference
temperature:

```@example metabolic
deb = sharpe_schoolfield_deb(activation = 0.65u"eV", reference_temperature = reference_temperature,
                             low_temperature = 5.0u"°C", low_deactivation = 3.0u"eV",
                             high_temperature = 40.0u"°C", high_deactivation = 7.0u"eV",
                             rate_at_reference = reference_rate)
deb(reference_temperature) == reference_rate, deb(35.0u"°C")
```

== Pawar

The metabolic curve of Pawar et al. (2016), which peaks below 40 °C:

```@example metabolic
pawar_model = pawar(rate_at_reference = reference_rate, activation_energy = 0.65u"eV",
                    deactivation_energy = 4.0u"eV", peak_temperature = 40.0u"°C",
                    reference_temperature = reference_temperature)
pawar_model(35.0u"°C")
```

:::

The rates from each at a range of body temperatures, in mW:

```@example metabolic
responses = ("Andrews and Pough" => T -> standard_metabolic_rate(Squamate(), 200.0u"g", T),
             "Arrhenius" => mte_rate, "Sharpe-Schoolfield" => full, "DEB" => deb, "Pawar" => pawar_model)
body_temperatures = [5.0, 15.0, 25.0, 35.0, 42.0] .* u"°C"
markdown_table(["Body temperature"; collect(first.(responses))],
               [[T; [uconvert(u"mW", f(T)) for f in last.(responses)]] for T in body_temperatures])
```

```@example metabolic
temperatures = collect(0.0:0.25:48.0) .* u"°C"
fig, ax = figure_axis("Body temperature (°C)", "Standard metabolic rate (mW)")
for (label, response) in responses
    lines!(ax, ustrip.(u"°C", temperatures), ustrip.(u"mW", response.(temperatures)); linewidth = 2, label)
end
axislegend(ax; position = :lt)
fig
```

All pass through, or near, the same rate at 25 °C, since they share the reference rate from the allometry, and differ
in how the rate changes away from it.

## Plants

BiologicalScaling.jl gives the dark respiration of C3 plants per unit biomass with the Arrhenius model and an
activation energy of 0.65 eV, normalised to 25 °C (Reich et al. 2006). This is the same as an
[`ArrheniusModel`](@ref) of ThermalPhysiology.jl:

```@example metabolic
biomass = 2.0u"kg"
respiration = ArrheniusModel(0.65u"eV"; T_ref = 25.0u"°C")
respiration_25 = allometric(BasalMetabolicRate(), C3Plant(), biomass, 25.0u"°C")
markdown_table(["Temperature", "BiologicalScaling.jl", "ThermalPhysiology.jl"],
               [(T, uconvert(u"mW", allometric(BasalMetabolicRate(), C3Plant(), biomass, T)),
                 uconvert(u"mW", respiration_25 * respiration(T))) for T in [5.0, 15.0, 25.0, 35.0] .* u"°C"])
```
