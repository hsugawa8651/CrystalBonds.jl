# CrystalBonds.jl

Chemical bonding analysis for crystalline materials in Julia.

## Overview

CrystalBonds.jl computes Wannier Orbital Hamilton Population (WOHP) and related
quantities from Wannier90-compatible tight-binding Hamiltonians. It provides
energy-resolved bonding/antibonding analysis analogous to
[LOBSTER](http://www.cohp.de/)'s projected COHP (pCOHP), but using Wannier
functions as the local basis.

### Features (v0.1)

- **WOHP**: Wannier Orbital Hamilton Population (= -pCOHP in Deringer convention)
- **WOOP**: Wannier Orbital Overlap Population (= projected DOS)
- **IpCOHP**: Integrated WOHP for quantitative bond strength
- **Orbital decomposition**: s-s, s-p, p-p resolved bonding analysis
- **Minimal dependencies**: only `LinearAlgebra` (stdlib)

### Input

CrystalBonds.jl reads Wannier90 `_hr.dat` files directly. No dependency on
[DFTK.jl](https://github.com/JuliaMolSim/DFTK.jl),
[Wannier.jl](https://github.com/qiaojunfeng/Wannier.jl),
or any DFT package. Any DFT code that produces
[Wannier90](https://wannier.org/) output
(VASP, Quantum ESPRESSO, ABINIT, DFTK) can be used as input.

### References

- Kundu et al., "Population Analysis with Wannier Orbitals",
  J. Chem. Phys. (2020). DOI: [10.1063/5.0032605](https://doi.org/10.1063/5.0032605)
- Deringer et al., "Crystal Orbital Hamilton Population (COHP) Analysis
  As Projected from Plane-Wave Basis Sets",
  J. Phys. Chem. A 115, 5461 (2011). DOI: [10.1021/jp202489s](https://doi.org/10.1021/jp202489s)

See [API Reference](api.md) for the full list of exported functions and types.

### Planned

- Tetrahedron method (replacing Gaussian smearing)
- ICOBI (crystal orbital bond index)
- Mulliken / Löwdin charges
- GPU support
