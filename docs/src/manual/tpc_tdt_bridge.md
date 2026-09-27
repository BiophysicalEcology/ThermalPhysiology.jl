# Linking performance and death

Thermal performance curves and thermal death time curves are usually measured and modelled separately, but they share
a mathematical core in the Arrhenius rate law, which links the breadth of a performance curve to the slope of a thermal
death time curve (Arnoldi et al. 2025).

```@setup bridge
using Main.FigureHelpers
using CairoMakie, ThermalPhysiology, Unitful
```

## From Arrhenius to thermal death time

If heat damage accumulates at a rate following the Arrhenius law,

```math
R(T) \propto \exp\left(\frac{T_A}{T_{ref}} - \frac{T_A}{T}\right)
```

then the time to a critical amount of damage ``D_{crit}`` is

```math
t = \frac{D_{crit}}{R(T)} \quad\Rightarrow\quad \ln t = \text{const} - T_A\left(\frac{1}{T_{ref}} - \frac{1}{T}\right)
```

Over the range of biological temperatures, where ``T - T_{ref} \ll T_{ref}``, this is close to linear in temperature,

```math
\ln t \approx \text{const} - \frac{T - T_{ref}}{E}, \qquad E = \frac{T_{ref}^2}{T_A}
```

which is the log-linear thermal death time model with a z-value, the rise in temperature for a ten-fold fall in time,
of

```math
z = E \ln 10 \approx 2.303\,E
```

``E`` is the thermal breadth of the universal thermal performance curve, so a narrow performance curve, a thermal
specialist, has a steep thermal death time curve. The relation is exact for the Arrhenius model, and approximate for the
Sharpe-Schoolfield models.

## Conversions

::: tabs

== Performance to TDT

[`tdt_from_tpc`](@ref) builds a [`LogLinearTDTModel`](@ref) from a [`UniversalTPCModel`](@ref), with ``z = E \ln 10``
and the critical thermal maximum of the performance curve as the static ``CT_{max}`` at the reference duration:

```@example bridge
tpc = utpc(optimal_temperature = 30.0u"°C", thermal_breadth = 4.0u"K", maximum_performance = 1.0)
tdt = tdt_from_tpc(tpc; reference_duration = 60.0u"minute", incipient_temperature = 30.0u"°C")
nothing # hide
```

```@eval
using ThermalPhysiology, Unitful, Markdown
tpc = utpc(optimal_temperature = 30.0u"°C", thermal_breadth = 4.0u"K", maximum_performance = 1.0)
tdt = tdt_from_tpc(tpc; reference_duration = 60.0u"minute", incipient_temperature = 30.0u"°C")
Markdown.parse("""
| Parameter | Value | From the performance curve |
|:----------|:------|:---------------------------|
| z-value | $(round(u"K", tdt.z_value; digits = 2)) | ``E \\ln 10`` |
| Static ``CT_{max}`` | $(round(u"°C", tdt.reference_ctmax; digits = 2)) | ``T_{opt} + E`` |
| Reference duration | $(tdt.reference_duration) | chosen |
| Incipient temperature | $(round(u"°C", tdt.incipient_temperature; digits = 2)) | chosen |
""")
```

== TDT to performance

[`thermal_breadth_from_tdt`](@ref) recovers the thermal breadth ``E = z / \ln 10`` from a thermal death time model:

```@example bridge
thermal_breadth_from_tdt(tdt)
```

== Arrhenius to z-value

[`z_value`](@ref) of an [`ArrheniusModel`](@ref) is ``T_{ref}^2 \ln 10 / T_A``:

```@example bridge
z_value(ArrheniusModel(T_A = 20000.0u"K", T_ref = 30.0u"°C"))
```

:::

An [`ArrheniusModel`](@ref) can also be used directly as a thermal death time model, with the survival time
``D_{crit} / f(T)``, where ``D_{crit}`` is given by the `critical_damage` keyword of [`survival_time`](@ref).

Thermal specialists and generalists, with narrow and broad performance curves, have steep and shallow thermal death time
curves:

```@example bridge
temperatures = collect(15.0:0.1:45.0) .* u"°C"
fig = Figure(size = (700, 600))
ax1 = Axis(fig[1, 1]; ylabel = "Relative performance")
ax2 = Axis(fig[2, 1]; xlabel = "Temperature (°C)", ylabel = "Survival time (minutes)", yscale = log10)
for (E, label) in ((2.0u"K", "Specialist"), (4.0u"K", "Intermediate"), (6.0u"K", "Generalist"))
    tpc = utpc(optimal_temperature = 30.0u"°C", thermal_breadth = E, maximum_performance = 1.0)
    tdt = tdt_from_tpc(tpc)
    lines!(ax1, ustrip.(temperatures), max.(0.0, tpc.(temperatures)); linewidth = 2, label = "$label, E = $E")
    hot = filter(>=(30.0u"°C"), temperatures)
    lines!(ax2, ustrip.(hot), ustrip.(u"minute", survival_time.(Ref(tdt), hot)); linewidth = 2)
end
axislegend(ax1; position = :lt)
linkxaxes!(ax1, ax2)
hidexdecorations!(ax1; grid = false)
fig
```
