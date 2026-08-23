# ── Abstract type hierarchy ────────────────────────────────────────────────────

"""
    AbstractTPCModel

Root abstract type for all thermal performance curve models.
Subtypes dispatch on [`thermal_performance`](@ref) or [`temperature_correction`](@ref).
"""
abstract type AbstractTPCModel end

"""
    AbstractArrheniusModel <: AbstractTPCModel

Abstract type for Arrhenius-family temperature-correction models.
These return a dimensionless correction factor equal to 1 at `T_ref`.
Parameters are stored in Kelvin (T_A, T_ref). DEBtool convention.
"""
abstract type AbstractArrheniusModel <: AbstractTPCModel end

"""
    AbstractPhenomenologicalModel <: AbstractTPCModel

Abstract type for phenomenological TPC shape models.
Parameters are typically in °C (T_opt, CTmax). Return absolute or relative performance.
"""
abstract type AbstractPhenomenologicalModel <: AbstractTPCModel end

"""
    AbstractTDTModel

Abstract type for thermal death time models.
Subtypes dispatch on [`survival_time`](@ref).
Kept separate from TPC types because correction factors and survival times are distinct quantities.
"""
abstract type AbstractTDTModel end

"""
    AbstractRepairModel

Abstract type for models of recovery from accumulated thermal injury
(see [`step_injury`](@ref)), evaluated in parallel with damage accrual.
Subtypes dispatch on [`repair_rate`](@ref) (a
continuous per-step reduction) and/or [`resets_injury`](@ref) (an
instantaneous full reset, e.g. below some threshold temperature).
"""
abstract type AbstractRepairModel end
