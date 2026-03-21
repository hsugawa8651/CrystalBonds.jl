# Quick Start

## Installation

```julia
using Pkg
Pkg.add("CrystalBonds")
```

## Basic workflow

```julia
using CrystalBonds

# 1. Read Wannier90 _hr.dat
HR, Rvectors, N = read_w90_hr("diamond_hr.dat")

# 2. Build k-mesh and interpolate H(R) -> H(k)
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

# 3. Compute WOHP
E_range = collect(range(minimum(evals) - 2, maximum(evals) + 2; length=500))
wohp = compute(WOHP(), Hk, evals, evecs; E_range=E_range)

# 4. Extract bond and integrate
cc_wohp = extract_bond(wohp, 1:4, 5:8)  # atom1 orbitals, atom2 orbitals
icohp = integrate(wohp, fermi_energy)
```

## Sign convention

CrystalBonds.jl follows the Kundu convention:

- **Positive WOHP** = bonding interaction
- **Negative WOHP** = antibonding interaction

This is equivalent to `-pCOHP` in the Deringer/LOBSTER convention.

## Orbital indices

For a system with `s + 3p` orbitals per atom (e.g., Diamond with 2 C atoms):

| Index | Orbital |
| ----- | ------- |
| 1     | atom 1, s |
| 2     | atom 1, px |
| 3     | atom 1, py |
| 4     | atom 1, pz |
| 5     | atom 2, s |
| 6     | atom 2, px |
| 7     | atom 2, py |
| 8     | atom 2, pz |

Orbital decomposition:

```julia
wohp_ss = wohp.matrix[1, 5, :]                                            # s-s
wohp_sp = sum(wohp.matrix[1:1, 6:8, :]; dims=(1,2))[:] .+
          sum(wohp.matrix[2:4, 5:5, :]; dims=(1,2))[:]                    # s-p + p-s
wohp_pp = sum(wohp.matrix[2:4, 6:8, :]; dims=(1,2))[:]                   # p-p
```

## Integration methods

Currently supported:

```julia
GaussianSmearing(0.15)  # Gaussian broadening with sigma = 0.15 eV (default)
```

Tetrahedron method is planned for a future release.
