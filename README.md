# ThermalPhysiology.jl

[![CI](https://github.com/BiophysicalEcology/ThermalPhysiology.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/BiophysicalEcology/ThermalPhysiology.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/BiophysicalEcology/ThermalPhysiology.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/BiophysicalEcology/ThermalPhysiology.jl)
[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://BiophysicalEcology.github.io/ThermalPhysiology.jl/stable/)
[![Docs dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://BiophysicalEcology.github.io/ThermalPhysiology.jl/dev/)

Thermal performance curve (TPC) and thermal death time (TDT) models for biophysical ecology,
built for the [BiophysicalEcology.jl](https://github.com/BiophysicalEcology) ecosystem.

## Models

**Arrhenius family** (`temperature_correction`): `ArrheniusModel`, `SharpSchoolFullModel`,
`SharpSchoolDEBModel`, `SharpSchoolHighModel`, `SharpSchoolLowModel`, `JohnsonLewinModel`

**Phenomenological TPC** (`thermal_performance`): `UniversalTPCModel` (Arnoldi et al. 2025),
`DeutschModel`, `Briere1Model`, `Briere2Model`, `GaussianModel`, `Thomas2012Model`,
`Thomas2017Model`, `PawarModel`, `Lactin2Model`

**TDT** (`survival_time`): `LogLinearTDTModel` (Jørgensen et al. 2021),
`ToleranceLandscape` (Rezende et al. 2020)

## Quick start

```julia
using ThermalPhysiology, Unitful

# Universal TPC (Arnoldi et al. 2025)
m = utpc(optimal_temperature=30.0u"°C", thermal_breadth=10.0u"K")
thermal_performance(m, 30.0)   # → 1.0 (peak)
optimal_temperature(m)         # → 303.15 K
critical_thermal_maximum(m)    # → 40.0 °C

# Arrhenius temperature correction (DEBtool convention)
m_arr = arrhenius(activation_energy=0.65u"eV", reference_temperature=20.0u"°C")
temperature_correction(m_arr, 30.0)   # → correction factor at 30°C

# Sharpe-Schoolfield full model
m_ss = sharpe_schoolfield(activation_energy=0.65u"eV", reference_temperature=20.0u"°C",
                          low_temperature=0.0u"°C", low_deactivation_energy=1.5u"eV",
                          high_temperature=42.0u"°C", high_deactivation_energy=5.0u"eV")

# Thermal death time
m_tdt = log_linear_tdt(z_value=4.0, reference_ctmax=39.0, reference_duration=60.0)
survival_time(m_tdt, 41.0)   # → minutes to knockdown at 41°C

# Constant temperature equivalent (DEB/NicheMapR)
T_series = [20.0, 22.0, 25.0, 28.0, 25.0, 22.0, 20.0]
constant_temperature_equivalent(m_arr, T_series)

# Curve fitting (Schoolfield graphical initial estimates)
m_fit = fit_thermal_performance_curve(SharpSchoolFullModel, temps_C, rates;
                                      T_ref=298.15u"K", weights=1.0 ./ rates)
```

## References

- Schoolfield, Sharpe & Magnuson (1981) *J Theor Biol* 88:719
- Kooijman (2010) *Dynamic Energy Budget Theory*, §2.6
- Arnoldi, Jackson, Peralta-Maraver & Payne (2025) *PNAS*
- Jørgensen et al. (2021) *Sci Reports* 11:12962
- Rezende et al. (2020) *Science* 369:1242
