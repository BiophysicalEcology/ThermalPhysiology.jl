# Curve Fitting

All fitting functions use [LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl)
(nonlinear least squares) or ordinary least squares where appropriate.

## TPC fitting

### Generic models

```julia
fit_thermal_performance_curve(ModelType, temperatures, rates;
                              initial_parameters=nothing, weights=nothing)
```

Supported generic models: `UniversalTPCModel`, `GaussianModel`, `DeutschModel`.
Returns a fitted instance of `ModelType`. Initial parameters are estimated from
data heuristics if `initial_parameters` is `nothing`.

```julia
temps = [5, 10, 15, 20, 25, 30, 35, 40, 45]
rates = [0.1, 0.3, 0.6, 0.9, 1.0, 0.9, 0.5, 0.1, 0.0]

m_fit = fit_thermal_performance_curve(UniversalTPCModel, temps, rates)
```

### Sharpe-Schoolfield models

`SharpSchoolFullModel` and `SharpSchoolDEBModel` have a specialised overload with
initial parameter estimates from the Schoolfield et al. (1981) graphical method
(piecewise OLS on the Arrhenius plot):

```julia
fit_thermal_performance_curve(SharpSchoolFullModel, temperatures, rates;
                              T_ref=298.15u"K",
                              log_transform=true,
                              initial_parameters=nothing,
                              weights=nothing)
```

**`log_transform=true` (default):** minimises sum of squared log-residuals
``\sum(\log \hat{r}_i - \log r_i)^2``. Appropriate when rates span orders of
magnitude; matches the Schoolfield et al. (1981) Marquardt fit.

**`log_transform=false`:** minimises absolute residuals; use with `weights` for
weighted least squares (e.g. `weights = 1 ./ rates` for equal relative weighting
of each point).

**Graphical initial estimates (Schoolfield 1981, Figure 2):**

The data are split into cold / middle / hot fractions. OLS is fitted to
``\ln r \sim 1/T`` in each region. From the slopes:

```
T_A  ≈ −slope_middle
T_AL ≈ −slope_cold − T_A
T_AH ≈  slope_hot  + T_A
```

T_L and T_H are found from the intersections of the cold / hot OLS lines with the
half-Arrhenius line (middle line shifted down by ln 2).

## TDT fitting

### Static knockdown data

Provide group-summary mean knockdown times at each assay temperature:

```julia
data = StaticKnockdownData(temperatures=[36, 38, 40, 42, 44],
                           knockdown_times=[289, 72, 18, 4.5, 1.1])
m_tdt = fit_thermal_death_time_curve(data; reference_duration=60.0)
```

Fits ``\log_{10}(t) \sim T`` by OLS. Returns a `LogLinearTDTModel`.

### Dynamic CTmax data

Provide knockdown temperatures observed at multiple ramp rates:

```julia
data = DynamicKnockdownData(ramp_rates=[0.1, 0.25, 0.5],
                             dynamic_ctmax_values=[40.2, 41.8, 43.1],
                             start_temperature=20.0)
m_tdt = fit_thermal_death_time_curve(data; reference_duration=60.0)
```

Uses Jørgensen (2021) Eq. 7a in NLS (≥ 3 ramp rates) or a root-finding scan (2 rates).

### Tolerance landscape

Requires individual-level knockdown data (one row per organism):

```julia
data = IndividualKnockdownData(temperatures=assay_temps,
                               knockdown_times=times_to_knockdown)
tl = fit_tolerance_landscape(data; n_bins=1000)
```

**Algorithm (Rezende et al. 2014, `tolerance.landscape()`):**

1. OLS on ``\log_{10}(t) \sim T`` → z-value, T_max
2. Z-shift all knockdown times to the mean assay temperature:
   ``t_{\text{shifted}} = t \cdot 10^{(T - \bar{T}) / z}``
3. Compute empirical survival function S(τ) on shifted times
4. Store as `n_bins`-row matrix (interpolated to regular grid)

The resulting `ToleranceLandscape` can then be used for `dynamic_survival` predictions.

## Data container types

```julia
IndividualKnockdownData(; temperatures, knockdown_times)
StaticKnockdownData(; temperatures, knockdown_times)
DynamicKnockdownData(; ramp_rates, dynamic_ctmax_values, start_temperature)
```

All temperature arguments accept Unitful quantities or bare floats (°C assumed).
