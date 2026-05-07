# ThermalPhysiology.jl

A Julia package for thermal performance curve (TPC) and thermal death time (TDT) models,
built for the [BiophysicalEcology.jl](https://github.com/BiophysicalEcology) ecosystem.
Designed to be used alongside NicheMapR and DEBtool\_J for mechanistic niche modelling and
Dynamic Energy Budget (DEB) theory.

## Design

The package is built on three principles:

**Multiple dispatch on model structs.** Each model is a `@kwdef struct` containing its
parameters. The same function name (`thermal_performance`, `temperature_correction`,
`survival_time`) dispatches to the appropriate formula for each model type. Models are
also callable directly as functions:

```julia
m = utpc(optimal_temperature=30.0u"°C", thermal_breadth=10.0u"K")
m(25.0)                         # thermal_performance(m, 25.0)
thermal_performance(m, 25.0)    # same thing
```

**Unitful throughout.** Temperatures can be passed as bare floats (assumed °C for
phenomenological models, Kelvin for Arrhenius models) or as `Unitful` quantities.
Property functions return `Unitful` quantities. Internal arithmetic uses bare Kelvin.

**Two distinct model families** with a mathematical bridge between them:

| Family | Primary function | Output | Temperature convention |
|--------|-----------------|--------|----------------------|
| Arrhenius | `temperature_correction` | Dimensionless factor (= 1 at T\_ref) | Kelvin |
| Phenomenological | `thermal_performance` | Absolute or relative rate | °C |
| TDT | `survival_time` | Minutes to knockdown | °C |

## Mathematical background

### TPC ↔ TDT bridge

The Arrhenius rate law provides the shared mathematical core of both TPC correction
factors and TDT damage rates. If the damage rate scales as:

```
R(T) ∝ exp(T_A/T_ref - T_A/T)
```

then the time to a critical cumulative damage threshold is:

```
t_death = D_crit / R(T)  ⟹  log(t_death) = const - T_A·(1/T_ref - 1/T)
```

At biological temperatures where ``\Delta T \ll T^2``, this linearises to:

```
log(t_death) ≈ const - (T - T_ref) / E
```

where ``E = T_{\text{ref}}^2 / T_A`` (in Kelvin) is the **thermal breadth** — the same
parameter as in the Universal TPC. Consequently:

- TDT slope ``b = 1/E``; narrow TPC breadth → steep TDT → thermal specialist
- z-value (°C/decade) ``= E \times \ln(10) \approx 2.303\,E``
- These relationships are exact for `ArrheniusModel`; approximate for Schoolfield variants

### Universal TPC (Arnoldi et al. 2025)

A two-parameter, scale-invariant TPC:

```math
y(x) = e^x (1 - x), \quad x = (T - T_{\text{opt}}) / E
```

CTmax = T\_opt + E (where y = 0); CTmin → −∞ (asymptotic approach).
The peak is always 1 at x = 0 (scaled by `maximum_performance`).

### Sharpe-Schoolfield full model (Schoolfield et al. 1981)

Enzyme kinetics with both cold and heat inactivation:

```math
r(T) = \rho(T_{\text{ref}}) \cdot \frac{T}{T_{\text{ref}}} \cdot
       \frac{e^{T_A/T_{\text{ref}} - T_A/T}}
       {1 + e^{T_{AL}/T - T_{AL}/T_L} + e^{T_{AH}/T_H - T_{AH}/T}}
```

The ``T/T_{\text{ref}}`` pre-factor from collision-frequency theory is included
(following the original paper and rTPC). `SharpSchoolDEBModel` omits this factor and
normalises so that `temperature_correction(m, T_ref) == rate_at_reference` exactly
(DEBtool `tempcorr` convention).

### Log-linear TDT model (Jørgensen et al. 2021)

```math
t(T) = t_{\text{ref}} \cdot 10^{(T_{\text{CTmax}} - T) / z}
```

z (°C/decade) is the temperature increment for a 10-fold change in knockdown time.
Reference CTmax is the static knockdown temperature at duration `t_ref`.

### Constant temperature equivalent (CTE)

The single constant temperature that produces the same mean thermal response as a
fluctuating temperature series ``\{T_i\}``:

```
CTE: find T_eq such that f(T_eq) = mean_i[f(T_i)]
```

Used in DEBtool/NicheMapR to convert field body temperatures to a DEB-compatible
constant temperature. By Jensen's inequality, CTE > arithmetic mean for convex models
(Arrhenius). Analytic solution available for `ArrheniusModel`; numerical (Roots.jl)
otherwise.

## Quick start

```julia
using ThermalPhysiology, Unitful

# ── Thermal performance curves ──────────────────────────────────────────────

# Universal TPC (Arnoldi et al. 2025)
m = utpc(optimal_temperature=30.0u"°C", thermal_breadth=10.0u"K")
thermal_performance(m, 30.0)        # → 1.0 (peak)
optimal_temperature(m)              # → 303.15 K
critical_thermal_maximum(m)         # → 40.0 °C

# Arrhenius temperature correction (DEBtool convention, = 1 at T_ref)
m_arr = arrhenius(activation_energy=0.65u"eV", reference_temperature=20.0u"°C")
temperature_correction(m_arr, 30.0) # → correction factor at 30°C

# Sharpe-Schoolfield (full Schoolfield 1981 model)
m_ss = sharpe_schoolfield(
    activation_energy        = 0.65u"eV",
    reference_temperature    = 20.0u"°C",
    low_temperature          = 0.0u"°C",
    low_deactivation_energy  = 2.0u"eV",
    high_temperature         = 42.0u"°C",
    high_deactivation_energy = 5.0u"eV",
)

# ── TPC properties ───────────────────────────────────────────────────────────
optimal_temperature(m)         # Unitful temperature
critical_thermal_maximum(m)    # Unitful temperature
critical_thermal_minimum(m)    # Unitful temperature (−∞ for UTPC)
thermal_breadth(m)             # CTmax − CTmin (Unitful K)
maximum_rate(m)                # peak thermal performance
q10(m, 20.0)                   # temperature coefficient at 20°C

# ── Thermal death time ───────────────────────────────────────────────────────
m_tdt = log_linear_tdt(z_value=4.0, reference_ctmax=39.0, reference_duration=60.0)
survival_time(m_tdt, 41.0)          # → minutes to knockdown at 41°C
lethal_temperature(m_tdt, 120.0)    # → temperature lethal in 120 min
ctmax_at_duration(m_tdt, 10.0)      # → sCTmax at 10-min exposure

# Dynamic and fluctuating exposures
ramp = 0.25  # °C/min
dynamic_ctmax(m_tdt, ramp)                          # predicted ramp CTmax
static_ctmax_from_dynamic(m_tdt, 42.5, ramp)        # recover sCTmax from dCTmax

T_series = collect(28.0:0.5:43.0)
injury = accumulated_injury(m_tdt, T_series, 1.0)   # injury vs time
time_to_failure(m_tdt, T_series, 1.0)               # minutes until injury = 1

# ── ToleranceLandscape (Rezende et al. 2020) ─────────────────────────────────
# Build from individual knockdown data and predict dynamic survival
data = IndividualKnockdownData(temperatures=rand(100) .* 8 .+ 36,
                               knockdown_times=rand(100) .* 120)
tl = fit_tolerance_landscape(data)
surv = dynamic_survival(tl, 28.0 .+ 10.0 .* sin.(range(0, π, 60)))

# ── Constant temperature equivalent ──────────────────────────────────────────
T_body = 20.0 .+ 8.0 .* sin.(range(0, 2π, 24))
constant_temperature_equivalent(m_arr, T_body)   # CTE > mean(T_body)

# ── TPC ↔ TDT bridge ─────────────────────────────────────────────────────────
m_tpc = utpc(optimal_temperature=30.0u"°C", thermal_breadth=10.0u"K")
m_tdt2 = tdt_from_tpc(m_tpc)           # z ≈ E × log(10)
thermal_breadth_from_tdt(m_tdt2)        # recover E from TDT
```

## Model registry

All models are listed in `THERMAL_REGISTRY`:

```julia
model_names()                          # all model keys
model_names(family=:arrhenius)         # filter by family
model_names(family=:phenomenological)
model_names(family=:tdt)
thermal_models(family=:tdt)            # returns Dict of entries
```

## References

- Schoolfield, Sharpe & Magnuson (1981) *J Theor Biol* 88:719–731
- Kooijman (2010) *Dynamic Energy Budget Theory*, 3rd ed., §2.6
- Arnoldi, Jackson, Peralta-Maraver & Payne (2025) *PNAS* 122:e2416587122
- Deutsch et al. (2008) *PNAS* 105:6668–6672
- Pawar et al. (2018) *Funct Ecol* 32:1874–1886
- Rezende et al. (2014) *Funct Ecol* 28:799–809
- Rezende et al. (2020) *Science* 369:1242–1245
- Jørgensen et al. (2021) *Sci Reports* 11:12962
