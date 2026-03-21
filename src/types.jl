# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara

# --- Bond analysis types ---

"""
    AbstractBondAnalysis

Abstract supertype for bond analysis methods.
"""
abstract type AbstractBondAnalysis end

"""
    WOHP <: AbstractBondAnalysis
    WOHP()

Wannier Orbital Hamilton Population.
Positive values indicate bonding, negative values indicate antibonding.
Equivalent to `-pCOHP` in the Deringer/LOBSTER convention.
"""
struct WOHP <: AbstractBondAnalysis end

"""
    WOOP <: AbstractBondAnalysis
    WOOP()

Wannier Orbital Overlap Population.
Equivalent to projected density of states (PDOS).
Sum over all orbitals gives total DOS.
"""
struct WOOP <: AbstractBondAnalysis end

# --- Integration method types ---

"""
    IntegrationMethod

Abstract supertype for Brillouin zone integration methods.
"""
abstract type IntegrationMethod end

"""
    GaussianSmearing <: IntegrationMethod
    GaussianSmearing(sigma=0.15)

Gaussian smearing for delta function approximation.
`sigma` is the broadening width in eV.
"""
struct GaussianSmearing <: IntegrationMethod
    sigma::Float64
end
GaussianSmearing() = GaussianSmearing(0.15)

# --- Result types (GPU-ready: parametric over array type) ---

"""
    BondAnalysisResult{T, A, V}

Result of a bond analysis computation ([`compute`](@ref)).

# Fields

  - `matrix::A`: `(n_wan, n_wan, n_E)` orbital-pair resolved values
  - `total::V`: `(n_E,)` sum over all orbital pairs
  - `E_range::V`: `(n_E,)` energy grid in eV

Type parameters are generic over `AbstractArray` for GPU compatibility.
"""
struct BondAnalysisResult{T <: Real, A <: AbstractArray{T, 3}, V <: AbstractVector{T}}
    matrix::A
    total::V
    E_range::V
end

# --- Error fallbacks ---

function compute(analysis::AbstractBondAnalysis, args...; kwargs...)
    error("compute not implemented for $(typeof(analysis))")
end

# read_w90_hr: implemented in io.jl

# fourier_interpolate, diag_Hk: implemented in fourier.jl

# extract_bond, integrate: implemented in analysis.jl

export AbstractBondAnalysis, WOHP, WOOP
export IntegrationMethod, GaussianSmearing
export BondAnalysisResult
export compute, read_w90_hr, fourier_interpolate, diag_Hk, gaussian_delta
export extract_bond, integrate
