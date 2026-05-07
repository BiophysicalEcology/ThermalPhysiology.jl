# Models

## Arrhenius family

These models return a **dimensionless correction factor** (= 1 at `T_ref`).
Parameters are stored in Kelvin. Use `temperature_correction(m, T)` or call `m(T)`.

### `ArrheniusModel`

```math
f(T) = \exp\!\left(\frac{T_A}{T_{\text{ref}}} - \frac{T_A}{T}\right)
```

Single-parameter model (Kooijman 2010). `T_A` is the Arrhenius temperature (K);
`T_ref` is the reference temperature where `f = 1`.

Constructor: `arrhenius(; activation_energy=0.65u"eV", reference_temperature=293.15u"K")`

---

### `SharpSchoolHighModel`

Arrhenius with high-temperature enzyme deactivation:

```math
r(T) = \rho_{\text{ref}} \cdot
       \frac{e^{T_A/T_{\text{ref}} - T_A/T}}
       {1 + e^{T_{AH}/T_H - T_{AH}/T}}
```

Constructor: `sharpe_schoolfield_high(; activation_energy, reference_temperature, high_temperature, high_deactivation_energy, rate_at_reference)`

---

### `SharpSchoolLowModel`

Arrhenius with low-temperature enzyme suppression:

```math
r(T) = \rho_{\text{ref}} \cdot
       \frac{e^{T_A/T_{\text{ref}} - T_A/T}}
       {1 + e^{T_{AL}/T - T_{AL}/T_L}}
```

Constructor: `sharpe_schoolfield_low(; ..., low_temperature, low_deactivation_energy, ...)`

---

### `SharpSchoolFullModel`

Full Schoolfield (1981) model — both low and high inactivation,
with the ``T/T_{\text{ref}}`` pre-factor from the original paper (Eq. 4):

```math
r(T) = \rho_{\text{ref}} \cdot \frac{T}{T_{\text{ref}}} \cdot
       \frac{e^{T_A/T_{\text{ref}} - T_A/T}}
       {1 + e^{T_{AL}/T - T_{AL}/T_L} + e^{T_{AH}/T_H - T_{AH}/T}}
```

`rate_at_reference` is the rate **in the absence of inactivation** at `T_ref`
(T\_ref assumed to be in the central Arrhenius zone).

Constructor: `sharpe_schoolfield(; activation_energy, reference_temperature, low_temperature, low_deactivation_energy, high_temperature, high_deactivation_energy, rate_at_reference)`

---

### `SharpSchoolDEBModel`

DEBtool normalised Schoolfield model (`tempcorr` convention).
No ``T/T_{\text{ref}}`` factor; instead the numerator is corrected so that
`temperature_correction(m, T_ref) == rate_at_reference` exactly for any parameter values:

```math
r(T) = \rho_{\text{ref}} \cdot e^{T_A/T_{\text{ref}} - T_A/T} \cdot
       \frac{1 + e^{T_{AL}/T_{\text{ref}} - T_{AL}/T_L} + e^{T_{AH}/T_H - T_{AH}/T_{\text{ref}}}}
       {1 + e^{T_{AL}/T - T_{AL}/T_L} + e^{T_{AH}/T_H - T_{AH}/T}}
```

Use this variant when computing CTE for DEB parameter estimation (NicheMapR).

Constructor: `sharpe_schoolfield_deb(; ...)`

---

### `JohnsonLewinModel`

Original 1946 enzyme kinetics model — structurally identical to `SharpSchoolHighModel`:

```math
r(T) = \rho_{\text{ref}} \cdot
       \frac{e^{T_A/T_{\text{ref}} - T_A/T}}
       {1 + e^{T_{AH}/T_H - T_{AH}/T}}
```

Constructor: `johnson_lewin(; ...)`

---

## Phenomenological TPC family

These models use `thermal_performance(m, T)` (or `m(T)`) with temperatures in °C
(bare floats or Unitful).

### `UniversalTPCModel`

Scale-invariant two-parameter TPC (Arnoldi et al. 2025):

```math
y(x) = e^x(1 - x), \quad x = \frac{T - T_{\text{opt}}}{E}
```

Peak at `T_opt` (where x = 0, y = 1); zero at `T_opt + E` (CTmax); asymptotic
approach to zero at cold temperatures (CTmin = −∞).

Constructor: `utpc(; optimal_temperature, thermal_breadth, maximum_performance=1.0)`

Stored internally: `T_opt` in Kelvin (bare Float64), `E` in Kelvin.

---

### `DeutschModel`

Gaussian rise below optimum, quadratic decline above:

```math
r(T) = \begin{cases}
  r_{\max} \exp\!\left(-\left(\frac{T - T_{\text{opt}}}{2\sigma}\right)^2\right) & T < T_{\text{opt}} \\[4pt]
  r_{\max} \left(1 - \left(\frac{T - T_{\text{opt}}}{T_{\text{opt}} - T_{\max}}\right)^2\right) & T \geq T_{\text{opt}}
\end{cases}
```

Constructor: `deutsch(; maximum_rate, optimal_temperature, critical_thermal_maximum, width_parameter)`

---

### `Briere1Model` / `Briere2Model`

Brière (1999) models for insect development rates:

```math
r(T) = c \cdot T \cdot (T - T_{\min}) \cdot \sqrt{T_{\max} - T}
\quad \text{(Brière 1)}
```

```math
r(T) = c \cdot T \cdot (T - T_{\min}) \cdot (T_{\max} - T)^{1/d}
\quad \text{(Brière 2)}
```

Constructors: `briere(; rate_constant, minimum_temperature, maximum_temperature [, shape_parameter])`

---

### `GaussianModel`

Symmetric Gaussian:

```math
r(T) = r_{\max} \exp\!\left(-\frac{1}{2}\left(\frac{T - T_{\text{opt}}}{\sigma}\right)^2\right)
```

Constructor: `gaussian(; maximum_rate, optimal_temperature, width_parameter)`

---

### `Thomas2012Model`

Quadratic TPC (Thomas et al. 2012):

```math
r(T) = c \cdot (T - T_{\text{opt}} + b) \cdot (b - (T - T_{\text{opt}}))
```

Constructor: `thomas(; rate_constant, shape_parameter, optimal_temperature)` (2012 variant)

---

### `Thomas2017Model`

Asymmetric skewed-Gaussian (Thomas et al. 2017):

```math
r(T) = r_{\max} \exp\!\left(-\frac{(T-T_{\text{opt}})^2}{2\sigma(T)^2}\right)
```

with ``\sigma(T) = \sigma_L`` below `T_opt`, ``\sigma_R`` above (controlled by `skewness`).

Constructor: `thomas(; maximum_rate, optimal_temperature, width_parameter, skewness)` (2017 variant)

---

### `PawarModel`

Metabolic TPC with biologically-interpretable parameters (Pawar et al. 2018):
Arrhenius rise characterised by `activation_energy` (eV), with high-temperature
deactivation controlled by `deactivation_energy` and `peak_temperature`.

Constructor: `pawar(; rate_at_reference, activation_energy, deactivation_energy, peak_temperature, reference_temperature)`

---

### `Lactin2Model`

Modified Lactin (1995) model for insect development:

```math
r(T) = \exp(\rho T) - \exp\!\left(\rho T_{\max} - \frac{T_{\max} - T}{\delta}\right) + \lambda
```

Constructor: `lactin2(; rate_constant, maximum_temperature, delta_temperature, intercept)`

---

## TDT family

Use `survival_time(m, T)` (minutes to knockdown at constant temperature T °C).

### `LogLinearTDTModel`

Log-linear (Jørgensen et al. 2021):

```math
t(T) = t_{\text{ref}} \cdot 10^{(T_{\text{CTmax}} - T) / z}
```

**Parameters:**
- `z_value`: °C for a 10-fold change in knockdown time
- `reference_ctmax`: sCTmax at `reference_duration`
- `reference_duration`: exposure duration defining `reference_ctmax` (min)
- `incipient_temperature`: temperature below which injury is negligible

Constructor: `log_linear_tdt(; z_value, reference_ctmax, reference_duration=60.0, incipient_temperature=30.0)`

**Derived quantities:**
- `temperature_maximum(m)` — temperature where mean knockdown = 1 min (Rezende parameterisation)
- `ctmax_at_duration(m, t)` — sCTmax for any exposure duration
- `dynamic_ctmax(m, ramp_rate)` — predicted ramp CTmax (Jørgensen Eq. 7a)
- `static_ctmax_from_dynamic(m, dctmax, ramp_rate)` — recover sCTmax from dCTmax (Eq. 7b)

**Fluctuating exposures:**
- `accumulated_injury(m, T_series, dt)` — cumulative injury vector (Eq. 3)
- `time_to_failure(m, T_series, dt)` — minutes until injury = 1

---

### `ToleranceLandscape`

Semi-parametric distributional model (Rezende et al. 2014, 2020).
Stores the full empirical survival curve S(τ) z-shifted to a reference temperature.

**Fields:** `z_value`, `temperature_maximum`, `mean_assay_temperature`, `survival_curve` (n×2 matrix)

Constructor: `tolerance_landscape(; z_value, temperature_maximum, mean_assay_temperature, survival_curve)`
Usually created by `fit_tolerance_landscape(data::IndividualKnockdownData)`.

**Functions:**
- `dynamic_survival(tl, T_series; dt_minutes, recovery_model)` — iterative curve-shifting
- `daily_mortality(tl, T_series_24h)` — 1 − survival over one day
- `cumulative_survival(tl, days)` — product of daily survival fractions
