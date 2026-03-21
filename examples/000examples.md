# CrystalBonds.jl Examples

## Naming convention

- `1xx`: Wannier90 file route (`_hr.dat` → WOHP). Uses CrystalBonds.jl only.
- `2xx`: Wannier.jl route (DFTK.jl → Wannier.jl → WOHP). Requires additional packages.
- `7xx`: Visualization utilities.

## Available examples

| File | System | Route | Description |
| ---- | ------ | ----- | ----------- |
| `110_diamond_w90.jl` | Diamond | `_hr.dat` | Minimal example using only CrystalBonds.jl |
| `120_gaas_w90.jl` | GaAs | `_hr.dat` | GaAs with orbital decomposition (s-s, s-p, p-s, p-p) |
| `210_diamond_dftk_pbe.jl` | Diamond | Wannier.jl | Full pipeline: DFT → Wannier → WOHP |

## Running

```sh
# From the CrystalBonds.jl root directory
julia --project=. examples/110_diamond_w90.jl
```

For `2xx` examples, install additional dependencies first:

```sh
julia --project=. -e 'using Pkg; Pkg.add(["DFTK", "PseudoPotentialData", "Wannier", "CairoMakie"])'
```
