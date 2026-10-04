# Thermal performance curves

A thermal performance curve (TPC) describes how the rate of a process, such as locomotion, feeding, growth or
development, depends on body temperature. Performance typically rises slowly with temperature to a peak at the optimal
temperature, then falls steeply to zero at the critical thermal maximum. The phenomenological models of this page are
evaluated with [`thermal_performance`](@ref), or by calling the model.

```@setup tpc
using Main.FigureHelpers
using CairoMakie, ThermalPhysiology, Unitful
```

## Models

In the equations, ``T`` is the temperature in °C unless stated otherwise.

::: tabs

== Universal

The universal thermal performance curve of Arnoldi et al. (2025), [`UniversalTPCModel`](@ref), built with
[`utpc`](@ref), has the same shape for all organisms and processes once temperature is scaled by the thermal breadth
``E`` and performance by its maximum ``P_{max}``:

```math
P(T) = P_{max}\, e^{x} (1 - x), \qquad x = \frac{T - T_{opt}}{E}
```

The peak is ``P_{max}`` at ``T_{opt}``, and performance falls to zero at the critical thermal maximum
``T_{opt} + E``. Below the optimum, performance approaches zero asymptotically.

```@example tpc
universal = utpc(optimal_temperature = 30.0u"°C", thermal_breadth = 8.0u"K", maximum_performance = 1.0)
universal.([20.0, 30.0, 35.0] .* u"°C")
```

== Deutsch

[`DeutschModel`](@ref), built with [`deutsch`](@ref), is the curve of Deutsch et al. (2008), a Gaussian rise to the
optimum and a quadratic decline to the critical thermal maximum ``CT_{max}``:

```math
P(T) = \begin{cases}
  P_{max} \exp\left(-\left(\frac{T - T_{opt}}{2\sigma}\right)^2\right) & T < T_{opt} \\
  P_{max} \left(1 - \left(\frac{T - T_{opt}}{T_{opt} - CT_{max}}\right)^2\right) & T \ge T_{opt}
\end{cases}
```

```@example tpc
deutsch_model = deutsch(maximum_rate = 1.0, optimal_temperature = 30.0u"°C",
                        critical_thermal_maximum = 38.0u"°C", width_parameter = 4.0u"K")
deutsch_model.([20.0, 30.0, 35.0] .* u"°C")
```

== Gaussian

[`GaussianModel`](@ref), built with [`gaussian`](@ref), is symmetric about the optimum:

```math
P(T) = P_{max} \exp\left(-\frac{1}{2}\left(\frac{T - T_{opt}}{\sigma}\right)^2\right)
```

```@example tpc
gaussian_model = gaussian(maximum_rate = 1.0, optimal_temperature = 30.0u"°C", width_parameter = 5.0u"K")
gaussian_model.([20.0, 30.0, 35.0] .* u"°C")
```

== Brière

The models of Brière et al. (1999) for the development rates of insects are zero below a minimum temperature
``T_{min}`` and above a maximum ``T_{max}``. [`Briere1Model`](@ref), built with [`briere`](@ref), is

```math
r(T) = c\, T (T - T_{min}) \sqrt{T_{max} - T}
```

and [`Briere2Model`](@ref) generalises the square root to a shape parameter ``d``:

```math
r(T) = c\, T (T - T_{min}) (T_{max} - T)^{1/d}
```

```@example tpc
briere1 = briere(rate_constant = 2.5e-4, minimum_temperature = 10.0u"°C", maximum_temperature = 38.0u"°C")
briere2 = Briere2Model(rate_constant = 2.5e-4, minimum_temperature = 10.0u"°C",
                       maximum_temperature = 38.0u"°C", shape_parameter = 3.0)
briere1.([20.0, 30.0, 35.0] .* u"°C"), briere2.([20.0, 30.0, 35.0] .* u"°C")
```

== Thomas 2012

[`Thomas2012Model`](@ref), built with [`thomas`](@ref), is the quadratic curve of Thomas et al. (2012), zero outside
``T_{opt} \pm b``:

```math
r(T) = c\, (T - T_{opt} + b)(b - (T - T_{opt}))
```

```@example tpc
thomas2012 = thomas(rate_constant = 1 / 144, shape_parameter = 12.0u"K", optimal_temperature = 28.0u"°C")
thomas2012.([20.0, 28.0, 35.0] .* u"°C")
```

== Thomas 2017

[`Thomas2017Model`](@ref) is a skewed Gaussian, after Thomas et al. (2017), with different widths below and above
the optimum:

```math
P(T) = P_{max} \exp\left(-\frac{(T - T_{opt})^2}{2\sigma(T)^2}\right), \qquad
\sigma(T) = \sigma \left(1 + s\,\mathrm{sign}(T - T_{opt})\right)
```

A negative skewness ``s`` gives the usual steep decline above the optimum.

```@example tpc
thomas2017 = Thomas2017Model(maximum_rate = 1.0, optimal_temperature = 30.0u"°C",
                             width_parameter = 5.0u"K", skewness = -0.5)
thomas2017.([20.0, 30.0, 35.0] .* u"°C")
```

== Pawar

[`PawarModel`](@ref), built with [`pawar`](@ref), is the metabolic curve of Pawar et al. (2016), an Arrhenius rise
with activation energy ``E`` and inactivation above a temperature ``T_{pk}`` with energy ``E_D``, with ``T`` in K:

```math
r(T) = r_{ref} \frac{\exp\left(\frac{E}{k_B}\left(\frac{1}{T_{ref}} - \frac{1}{T}\right)\right)}
       {1 + \exp\left(\frac{E_D}{k_B}\left(\frac{1}{T_{pk}} - \frac{1}{T}\right)\right)}
```

```@example tpc
pawar_model = pawar(rate_at_reference = 1.0, activation_energy = 0.65u"eV", deactivation_energy = 4.0u"eV",
                    peak_temperature = 33.0u"°C", reference_temperature = 20.0u"°C")
pawar_model.([20.0, 30.0, 35.0] .* u"°C")
```

== Lactin

[`Lactin2Model`](@ref), built with [`lactin2`](@ref), is the second model of Lactin et al. (1995) for insect
development:

```math
r(T) = e^{\rho T} - \exp\left(\rho T_{max} - \frac{T_{max} - T}{\Delta}\right) + \lambda
```

The coefficients ``\rho`` and ``\lambda`` are defined against the value of ``T`` in °C, and are bare numbers.

```@example tpc
lactin = lactin2(rate_constant = 0.12, maximum_temperature = 38.0u"°C", delta_temperature = 5.0u"K",
                 intercept = -1.1)
lactin.([20.0, 30.0, 35.0] .* u"°C")
```

:::

Scaled to their maxima, the curves are:

```@example tpc
models = (universal => "Universal", deutsch_model => "Deutsch", gaussian_model => "Gaussian",
          briere1 => "Brière 1", thomas2012 => "Thomas 2012", thomas2017 => "Thomas 2017",
          pawar_model => "Pawar", lactin => "Lactin 2")
temperatures = collect(5.0:0.1:42.0) .* u"°C"
fig, ax = figure_axis("Temperature (°C)", "Relative performance"; limits = (nothing, (0, 1.05)))
for (model, label) in models
    performance = model.(temperatures)
    lines!(ax, ustrip.(u"°C", temperatures), performance ./ maximum(performance); linewidth = 2, label)
end
axislegend(ax; position = :lt)
fig
```

## Properties of a curve

The optimal temperature, the critical thermal minimum and maximum, the thermal breadth between them, the maximum rate
and the ``Q_{10}`` are found analytically where there is a solution, and numerically otherwise:

```@example tpc
curves = (universal => "Universal", deutsch_model => "Deutsch", gaussian_model => "Gaussian",
          thomas2017 => "Thomas 2017", pawar_model => "Pawar", lactin => "Lactin 2")
nothing # hide
```

::: tabs

== Optimum

[`optimal_temperature`](@ref) and [`maximum_rate`](@ref):

```@example tpc
[(label, optimal_temperature(model), maximum_rate(model)) for (model, label) in curves]
```

== Critical maximum

[`critical_thermal_maximum`](@ref):

```@example tpc
[(label, critical_thermal_maximum(model)) for (model, label) in curves]
```

== Critical minimum

[`critical_thermal_minimum`](@ref):

```@example tpc
[(label, critical_thermal_minimum(model)) for (model, label) in curves]
```

== Q10

[`q10`](@ref) at 20 °C:

```@example tpc
[(label, q10(model, 20.0u"°C")) for (model, label) in curves]
```

:::

The critical thermal limits are the temperatures where performance falls to a `threshold`, zero by default. The
universal, Gaussian, Deutsch, Thomas 2017 and Pawar curves only approach zero on one or both sides of the optimum, so
some of their limits are ``-\infty`` or not found (`NaN`). A threshold gives limits, and a breadth, for any curve, for example the
range of temperatures with at least half of the maximum performance:

```@example tpc
half_maximum(model) = (; threshold = 0.5 * maximum_rate(model))
[(label, critical_thermal_minimum(model; half_maximum(model)...),
  critical_thermal_maximum(model; half_maximum(model)...),
  thermal_breadth(model; half_maximum(model)...)) for (model, label) in curves[2:end]]
```

The threshold is in the units of performance, so here it is set to half of each curve's own maximum rate.

## Units of performance

The maximum performance, or rate constant, can be a bare number for relative performance, or a quantity with units,
which are kept in the output:

```@example tpc
sprint_speed = utpc(optimal_temperature = 34.0u"°C", thermal_breadth = 7.0u"K",
                    maximum_performance = 2.5u"m/s")
sprint_speed(30.0u"°C")
```
