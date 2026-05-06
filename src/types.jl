# ── Abstract type hierarchy ────────────────────────────────────────────────────

# TPC families
abstract type AbstractTPCModel end

# Temperature-correction models: return dimensionless factor (= 1 at T_ref).
# Parameters in Kelvin (T_A, T_ref). DEBtool convention.
abstract type AbstractArrheniusModel <: AbstractTPCModel end

# Performance-shape models: parameters in °C (T_opt, CTmax). Return absolute/relative performance.
abstract type AbstractPhenomenologicalModel <: AbstractTPCModel end

# TDT family — separate because correction factors and survival times are distinct quantities.
abstract type AbstractTDTModel end
