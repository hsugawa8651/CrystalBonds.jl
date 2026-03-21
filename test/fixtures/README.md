# Test Fixtures

Wannier90-compatible `_hr.dat` files for unit tests.

## Diamond (`diamond_hr.dat`)

| Parameter | Value |
| --------- | ----- |
| DFT code | DFTK.jl v0.7.20 |
| Wannier code | Wannier.jl v0.3.6 |
| Pseudopotential | PseudoDojo NC SR PBE v0.4.1 standard (C, Z_ion=4) |
| Ecut | 30 Ha |
| DFT k-grid | 8 x 8 x 8 |
| Lattice constant | a = 6.7403 Bohr (3.567 A) |
| Structure | FCC, 2-atom basis, positions [1/8,1/8,1/8] and [-1/8,-1/8,-1/8] |
| n_bands = n_wannier | 8 (isolated bands, max_localize) |
| Projections | Hydrogenic 2s + 2p on each C atom (Z=6) |
| R-vectors | Wigner-Seitz, 617 vectors |
| Degeneracy | `.ndegen` file included |

## GaAs (`gaas_hr.dat`)

| Parameter | Value |
| --------- | ----- |
| DFT code | DFTK.jl v0.7.20 |
| Wannier code | Wannier.jl v0.3.6 |
| Pseudopotential | HGH PBE ga-q3 (Z_ion=3), as-q5 (Z_ion=5) |
| Ecut | 20 Ha |
| DFT k-grid | 8 x 8 x 8 |
| Lattice constant | a = 10.6829 Bohr (5.653 A) |
| Structure | Zincblende (FCC), Ga at [1/8,1/8,1/8], As at [-1/8,-1/8,-1/8] |
| n_bands = n_wannier | 8 (isolated bands, max_localize) |
| Projections | Hydrogenic 4s + 4p: Ga (Z=31), As (Z=33) |
| R-vectors | Wigner-Seitz, 617 vectors |
| Degeneracy | `.ndegen` file included |

## Regenerating

To regenerate these files, use the script in the development workspace:

```sh
julia crystalbonds-jl-dev/generate_hr_dat.jl
```

This requires DFTK.jl, PseudoPotentialData.jl, and Wannier.jl.
