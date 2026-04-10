# API Reference

## Types

Bond analysis types and result containers.

```@docs
AbstractBondAnalysis
WOHP
WOOP
IntegrationMethod
GaussianSmearing
BondAnalysisResult
```

## Input

Read Wannier90 files and prepare k-space Hamiltonian.

```@docs
read_w90_hr
fourier_interpolate
diag_Hk
```

## Compute

Core computation of bond populations.

```@docs
compute(::WOHP, ::AbstractArray{<:Complex,3}, ::AbstractMatrix{<:Real}, ::AbstractArray{<:Complex,3})
compute(::WOOP, ::AbstractArray{<:Complex,3}, ::AbstractMatrix{<:Real}, ::AbstractArray{<:Complex,3})
gaussian_delta
```

## Output

Post-processing: bond extraction and integration.

```@docs
extract_bond
integrate
```
