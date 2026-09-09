# Models

## Arrhenius family

These models return a **dimensionless correction factor** (= 1 at `T_ref`).
Parameters (`T_A`, `T_ref`, `T_L`, `T_AL`, `T_H`, `T_AH`) are stored as `Unitful`
Kelvin quantities and must be constructed with units (e.g. `9000.0u"K"`) — bare
numbers are rejected. Use `temperature_correction(m, T)` or call `m(T)`; `T` itself
still accepts a bare float (assumed °C) or a `Unitful` quantity.

Each struct has a direct field-based constructor (e.g. `SharpSchoolFullModel(;
T_A=9000.0u"K", T_ref=293.15u"K", ...)`). There's no generic factory that picks the
model type for you — `SharpSchoolFullModel` and `SharpSchoolDEBModel` take identical
keyword sets but compute different formulas, so the type has to be named explicitly
at the call site (`SomeType(; kwargs...)` also works when `SomeType` is a variable
holding the type, since Julia types are directly callable). No parameter has a
default except `rate_at_reference` (see below); each struct's docstring gives an example.

**`rate_at_reference`** (all models except `ArrheniusModel`) defaults to `nothing` —
not a stand-in physical value, but a mode switch: omit it and `temperature_correction`
returns the bare dimensionless correction; supply a rate (bare or Unitful, e.g.
`0.1u"d^-1"`) and it returns `rate_at_reference * correction` instead — the corrected
rate, in whatever units you supplied. This dispatches on the struct's type parameter,
so both forms are fully type-stable.

### `ArrheniusModel`

```math
f(T) = \exp\!\left(\frac{T_A}{T_{\text{ref}}} - \frac{T_A}{T}\right)
```

Single-parameter model (Kooijman 2010). `T_A` is the Arrhenius temperature (K);
`T_ref` is the reference temperature where `f = 1`.

Constructor: `ArrheniusModel(; T_A=..., T_ref=...)`, or `ArrheniusModel(E_A; T_ref=...)`
for activation energy (`0.65u"eV"`, converted via `ea_to_ta`).

---

### `SharpSchoolHighModel`

Arrhenius with high-temperature enzyme deactivation:

```math
r(T) = \rho_{\text{ref}} \cdot
       \frac{e^{T_A/T_{\text{ref}} - T_A/T}}
       {1 + e^{T_{AH}/T_H - T_{AH}/T}}
```

Constructor: `sharpe_schoolfield_high(; activation, reference_temperature, high_temperature, high_deactivation, rate_at_reference)`.
`activation`/`high_deactivation` each accept a Unitful temperature or energy.

---

### `SharpSchoolLowModel`

Arrhenius with low-temperature enzyme suppression:

```math
r(T) = \rho_{\text{ref}} \cdot
       \frac{e^{T_A/T_{\text{ref}} - T_A/T}}
       {1 + e^{T_{AL}/T - T_{AL}/T_L}}
```

Constructor: `sharpe_schoolfield_low(; ..., low_temperature, low_deactivation, ...)` —
`low_deactivation` accepts a Unitful temperature or energy.

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

Constructor: `sharpe_schoolfield(; activation, reference_temperature, low_temperature, low_deactivation, high_temperature, high_deactivation, rate_at_reference)`.
`activation`/`low_deactivation`/`high_deactivation` each accept a Unitful temperature or energy.

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
(bare floats or Unitful). No parameter has a default — every field is biologically
meaningful and must be supplied deliberately; each struct's docstring gives an example.

### `UniversalTPCModel`

Scale-invariant two-parameter TPC (Arnoldi et al. 2025):

```math
y(x) = e^x(1 - x), \quad x = \frac{T - T_{\text{opt}}}{E}
```

Peak at `T_opt` (where x = 0, y = 1); zero at `T_opt + E` (CTmax); asymptotic
approach to zero at cold temperatures (CTmin = −∞).

Constructor: `utpc(; optimal_temperature, thermal_breadth, maximum_performance)`

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

Constructors: `briere(; rate_constant, minimum_temperature, maximum_temperature)` (Brière 1);
`Briere2Model(; rate_constant, minimum_temperature, maximum_temperature, shape_parameter)` (Brière 2, no named wrapper)

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

Constructor: `Thomas2017Model(; maximum_rate, optimal_temperature, width_parameter, skewness)` (no named wrapper)

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

Use `survival_time(m, T)` (time to knockdown at constant temperature T).
No parameter has a default; see each struct's docstring for an example.

For `LogLinearTDTModel`, every scalar time and temperature quantity is Unitful and
required — bare numbers are rejected:
- `z_value`, `reference_ctmax`, `incipient_temperature` (struct fields); `z_value`
  accepts K or °C (reinterpreted directly as a K-sized difference — see the
  struct's docstring), the others accept any Unitful temperature
- `reference_duration`, `dt` (in `accumulated_injury`/`resettable_injury`/
  `time_to_failure`), and `duration` (in `ctmax_at_duration`/`lethal_temperature`)
- `ramp_rate` (in `dynamic_ctmax`/`static_ctmax_from_dynamic`) must be a
  temperature/time quantity such as `0.1u"K/minute"`; `°C/time` is rejected too
  since °C is an affine unit
- `T` in `survival_time(m, T)` directly

`survival_time`/`time_to_failure`/`temperature_maximum`/`ctmax_at_duration`/
`dynamic_ctmax`/`static_ctmax_from_dynamic`/`lethal_temperature` all return Unitful
quantities. **Bulk** temperature series (`T_series` in `accumulated_injury`,
`resettable_injury`, `time_to_failure`, `constant_temperature_equivalent`) and
`step_injury`'s single `temperature` argument (the per-element engine those bulk
functions call) stay bare-or-Unitful, matching the rest of the package — root-finding
internals (`lethal_temperature`, `z_value`, `constant_temperature_equivalent`) also
work in bare Kelvin/Celsius numerics, per `HeatExchange.jl`'s convention.

### `LogLinearTDTModel`

Log-linear (Jørgensen et al. 2021):

```math
t(T) = t_{\text{ref}} \cdot 10^{(T_{\text{CTmax}} - T) / z}
```

**Parameters:**
- `z_value`: temperature difference for a 10-fold change in knockdown time
- `reference_ctmax`: sCTmax at `reference_duration`
- `reference_duration`: exposure duration defining `reference_ctmax` — Unitful time
- `incipient_temperature`: temperature below which injury is negligible

Constructor: `log_linear_tdt(; z_value, reference_ctmax, reference_duration, incipient_temperature)`

**Derived quantities:**
- `temperature_maximum(m)` — temperature where mean knockdown = 1 min (Rezende parameterisation)
- `ctmax_at_duration(m, duration)` — sCTmax for any exposure duration
- `dynamic_ctmax(m, ramp_rate)` — predicted ramp CTmax (Jørgensen Eq. 7a)
- `static_ctmax_from_dynamic(m, dctmax, ramp_rate)` — recover sCTmax from dCTmax (Eq. 7b)

**Fluctuating exposures:**
- `accumulated_injury(m, T_series, dt)` — cumulative injury vector (Eq. 3)
- `time_to_failure(m, T_series, dt)` — time until injury = 1

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
