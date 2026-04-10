# CrystalBonds.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://hsugawa8651.github.io/CrystalBonds.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://hsugawa8651.github.io/CrystalBonds.jl/dev/)
[![CI](https://github.com/hsugawa8651/CrystalBonds.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/hsugawa8651/CrystalBonds.jl/actions/workflows/CI.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19494353.svg)](https://doi.org/10.5281/zenodo.19494353)

Chemical bonding analysis for crystalline materials in Julia.

CrystalBonds.jl computes Wannier Orbital Hamilton Population (WOHP) from
Wannier90-compatible tight-binding Hamiltonians, providing energy-resolved
bonding/antibonding analysis for crystals.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/hsugawa8651/CrystalBonds.jl")
```

## Quick example

```julia
using CrystalBonds

# Read Wannier90 _hr.dat
HR, Rvectors, N = read_w90_hr("diamond_hr.dat")

# Interpolate H(R) -> H(k) on fine k-mesh
nk = 30
kpoints = zeros(Float64, 3, nk^3)
let idx = 1
    for i in 0:nk-1, j in 0:nk-1, kk in 0:nk-1
        kpoints[:, idx] = [i/nk, j/nk, kk/nk]
        idx += 1
    end
end
Hk = fourier_interpolate(HR, Rvectors, kpoints; N=N)
evals, evecs = diag_Hk(Hk)

# Compute WOHP
E_range = collect(range(minimum(evals) - 2, maximum(evals) + 2; length=500))
wohp = compute(WOHP(), Hk, evals, evecs; E_range=E_range)

# Extract bond and integrate
cc_wohp = extract_bond(wohp, 1:4, 5:8)
icohp = integrate(wohp, fermi_energy)
```

## Features (v0.1)

- **WOHP**: Wannier Orbital Hamilton Population
- **WOOP**: Wannier Orbital Overlap Population (projected DOS)
- **IpCOHP**: Integrated WOHP for bond strength
- **Orbital decomposition**: s-s, s-p, p-p resolved analysis
- **Minimal dependencies**: LinearAlgebra only

## Input

Any DFT code that produces Wannier90 `_hr.dat` files can be used:
VASP, Quantum ESPRESSO, ABINIT, DFTK.jl, etc.

See [`examples/`](examples/) for usage with DFTK.jl + Wannier.jl.

## Documentation

Documentation is available [here](https://hsugawa8651.github.io/CrystalBonds.jl/dev). Build locally:

```bash
cd docs
julia --project make.jl
```

## Contributing

Bug reports and feature requests are welcome via [GitHub Issues](https://github.com/hsugawa8651/CrystalBonds.jl/issues).
Before opening a pull request, start an issue or a discussion on the topic.
This project follows the [Julia Community Standards](https://julialang.org/community/standards/).

## References

- Kundu et al., "Population Analysis with Wannier Orbitals", J. Chem. Phys. (2020). DOI: [10.1063/5.0032605](https://doi.org/10.1063/5.0032605)
- Deringer et al., "Crystal Orbital Hamilton Population (COHP) Analysis As Projected from Plane-Wave Basis Sets", J. Phys. Chem. A 115, 5461 (2011). DOI: [10.1021/jp202489s](https://doi.org/10.1021/jp202489s)

## License

GPL-3.0-or-later. See [LICENSE](LICENSE) for details.
