```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "ThermalPhysiology.jl"
  text: "Thermal performance and tolerance"
  tagline: "thermal performance curves, Arrhenius temperature corrections and thermal death time models for biophysical ecology, with units."
  actions:
    - theme: brand
      text: Get Started
      link: /get_started
    - theme: alt
      text: View on Github
      link: https://github.com/BiophysicalEcology/ThermalPhysiology.jl
    - theme: alt
      text: API Reference
      link: /api

features:
  - title: 🔥 Arrhenius corrections
    details: <a class="highlight-link">Arrhenius</a>, Sharpe-Schoolfield and Johnson-Lewin temperature corrections of rates, including the DEBtool normalisation used in Dynamic Energy Budget models.
    link: /manual/arrhenius
  - title: 📈 Performance curves
    details: <a class="highlight-link">Thermal performance curves</a> including the universal TPC of Arnoldi et al., Deutsch, Brière, Gaussian, Thomas, Pawar and Lactin models, with optima, critical limits and breadths.
    link: /manual/performance_curves
  - title: ☠️ Thermal death time
    details: <a class="highlight-link">Log-linear thermal death time</a> models and <a class="highlight-link">tolerance landscapes</a>, static and dynamic CTmax, and injury accumulation under fluctuating temperatures.
    link: /manual/thermal_death
  - title: 🌡️ Fluctuating temperatures
    details: The <a class="highlight-link">constant temperature equivalent</a> of a varying body temperature, and mean rates and corrections over a temperature series.
    link: /manual/fluctuating
  - title: 🧮 Fitting
    details: Fit performance curves, thermal death time curves and tolerance landscapes to <a class="highlight-link">experimental data</a> by least squares.
    link: /manual/fitting
  - title: 📏 Units
    details: Every temperature, time and energy is a <a class="highlight-link">Unitful.jl</a> quantity, so °C, K, minutes, hours and eV can be mixed freely, and bare numbers are rejected rather than guessed.
    link: /manual/introduction
---
```

## How to install ThermalPhysiology.jl?

ThermalPhysiology.jl can be installed from the Julia REPL:

```julia
julia> using Pkg
julia> Pkg.add("ThermalPhysiology")
# or
julia> ] # ']' should be pressed
pkg> add ThermalPhysiology
```

If you want to use the latest unreleased version, you can run the following command:

```julia
julia> using Pkg
julia> Pkg.add(url = "https://github.com/BiophysicalEcology/ThermalPhysiology.jl")
```

## Manual

ThermalPhysiology.jl is a general package for computing how temperature affects rates or survival: fitting and
comparing performance and thermal death time curves, correcting physiological rates for temperature, and estimating
performance, injury and survival from measured or modelled body temperatures.

It is part of the [BiophysicalEcology](https://github.com/BiophysicalEcology) ecosystem for mechanistic niche
modelling, to be brought together in [NicheMapper.jl](https://github.com/BiophysicalEcology/NicheMapper.jl) (in
development). It turns the body temperatures computed with
[HeatExchange.jl](https://github.com/BiophysicalEcology/HeatExchange.jl),
[BiophysicalBehaviour.jl](https://github.com/BiophysicalEcology/BiophysicalBehaviour.jl), 
[Microclimate.jl](https://github.com/BiophysicalEcology/Microclimate.jl) and
[MicroclimateMapper.jl](https://github.com/BiophysicalEcology/MicroclimateMapper.jl), into performance, temperature
corrections of physiological rates and heat injury. It can be used to adjust metabolic rates from
[BiologicalScaling.jl](https://github.com/BiophysicalEcology/BiologicalScaling.jl), and will be incorporated into 
Dynamic Energy Budget (DEB) models via [DEBtool_J.jl](https://github.com/add-my-pet/DEBtool_J.jl). See the
[Introduction](manual/introduction.md) for more about the design of this package and how it integrates with 
the [BiophysicalEcology](https://github.com/BiophysicalEcology) ecosystem.
