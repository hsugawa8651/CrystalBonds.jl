# Generating `_hr.dat` Files

CrystalBonds.jl reads Wannier90-compatible `_hr.dat` files as input.
This page describes how to generate these files from various DFT codes.

## From DFTK.jl + Wannier.jl (Julia)

The following example generates `diamond_hr.dat` for diamond (2 C atoms, sp3 bonding).

- **DFT**: DFTK.jl with PseudoDojo norm-conserving PBE pseudopotential, Ecut = 30 Ha, 8x8x8 k-grid
- **Wannier**: 8 maximally localized Wannier functions (2s + 2p per C atom, isolated bands)
- **Output**: `diamond_hr.dat` (real-space Hamiltonian) + `diamond_hr.dat.ndegen` (degeneracy weights)

```julia
using DFTK, PseudoPotentialData, Wannier

# 1. Run DFT
a = 6.7403  # Bohr
lattice = a / 2 * [[0 1 1.]; [1 0 1.]; [1 1 0.]]
family = PseudoFamily("dojo.nc.sr.pbe.v0_4_1.standard.upf")
C = ElementPsp(:C, family)
model = model_PBE(lattice, [C, C], [ones(3)/8, -ones(3)/8])
basis = PlaneWaveBasis(model; Ecut=30, kgrid=[8, 8, 8])
scfres = self_consistent_field(basis; tol=1e-8,
    nbandsalg=DFTK.FixedBands(; n_bands_converge=8))

# 2. Wannierize
positions = [ones(3)/8, -ones(3)/8]
projs = [DFTK.HydrogenicWannierProjection(positions[i], 2, l, m, 6.0)
         for i in 1:2 for (l, m) in [(0,0), (1,1), (1,-1), (1,0)]]
wm = Wannier.Model(scfres; n_bands=8, n_wannier=8, projections=projs,
    fileprefix="diamond")
U = Wannier.max_localize(wm; f_tol=1e-7, g_tol=1e-5, max_iter=200)

# 3. Build H(R) and write _hr.dat
U_arr = isa(U, Array{<:Complex, 3}) ? U : stack(U)
Hk = Wannier.get_Hk(wm.E, U_arr)
Rvecs_ws = Wannier.get_Rvectors_ws(wm.lattice, wm.kgrid)
kRvectors = Wannier.KRVectors(wm.lattice, wm.kgrid, wm.kpoints, Rvecs_ws)
HR = Wannier.fourier(kRvectors, Hk)

Wannier.WannierIO.write_HH_R("diamond_hr.dat", HR, Int.(Rvecs_ws.R);
    N=Int.(Rvecs_ws.N))
```

This produces `diamond_hr.dat` and `diamond_hr.dat.ndegen`.

## From VASP + Wannier90

```
1. Run VASP with LWANNIER90 = .TRUE.
2. Run wannier90.x to generate seedname_hr.dat
3. Use the resulting _hr.dat directly with CrystalBonds.jl
```

## From Quantum ESPRESSO + Wannier90

```
1. Run pw.x for SCF + NSCF
2. Run pw2wannier90.x to generate .amn, .mmn, .eig
3. Run wannier90.x to generate seedname_hr.dat
4. Use the resulting _hr.dat directly with CrystalBonds.jl
```

## File format

The `_hr.dat` format is a plain text file:

```
header line
n_wan
n_Rvecs
R1 R2 R3  mu  nu  Re(H)  Im(H)
R1 R2 R3  mu  nu  Re(H)  Im(H)
...
```

An optional `.ndegen` file contains Wigner-Seitz degeneracy weights
(one integer per R-vector). If present, CrystalBonds.jl reads it automatically.
If absent, all weights default to 1.

## Bundled test data

CrystalBonds.jl ships with `_hr.dat` files for Diamond and GaAs in `test/fixtures/`.
See `test/fixtures/README.md` for DFT parameters and regeneration instructions.
