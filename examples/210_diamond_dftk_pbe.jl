# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
#
# Diamond WOHP via DFTK.jl + Wannier.jl pipeline
# Requires: DFTK, PseudoPotentialData, Wannier, CairoMakie (optional)
#
# Usage: julia --project=.. examples/210_diamond_dftk_pbe.jl

using DFTK
using PseudoPotentialData
using Wannier
using CrystalBonds
using LinearAlgebra

# === Step 0: DFT ===
a = 6.7403  # Bohr
lattice = a / 2 * [[0 1 1.0]; [1 0 1.0]; [1 1 0.0]]
family = PseudoFamily("dojo.nc.sr.pbe.v0_4_1.standard.upf")
C = ElementPsp(:C, family)
atoms = [C, C]
positions = [ones(3)/8, -ones(3)/8]

model = model_PBE(lattice, atoms, positions)
basis = PlaneWaveBasis(model; Ecut = 30, kgrid = [8, 8, 8])
println("Running SCF...")
scfres = self_consistent_field(basis; tol = 1e-8,
    nbandsalg = DFTK.FixedBands(; n_bands_converge = 8))
eF_eV = scfres.εF * 27.211386
println("Fermi energy: $(round(eF_eV, digits=4)) eV")

# === Step 1: Wannierize ===
C_Z = 6.0
projs = [
    DFTK.HydrogenicWannierProjection(positions[i], 2, l, m, C_Z)
    for i in 1:2 for (l, m) in [(0, 0), (1, 1), (1, -1), (1, 0)]
]
workdir = joinpath(@__DIR__, "wannier_work")
mkpath(workdir)
wannier_model = Wannier.Model(scfres;
    n_bands = 8, n_wannier = 8, projections = projs,
    fileprefix = joinpath(workdir, "diamond"))
U = Wannier.max_localize(wannier_model; f_tol = 1e-7, g_tol = 1e-5, max_iter = 200)

# === Step 2: Build H(R) -> H(k) using Wannier.jl ===
U_arr = isa(U, Array{<:Complex, 3}) ? U : stack(U)
Hk_coarse = Wannier.get_Hk(wannier_model.E, U_arr)
Rvecs_ws = Wannier.get_Rvectors_ws(wannier_model.lattice, wannier_model.kgrid)
kRvectors = Wannier.KRVectors(wannier_model.lattice, wannier_model.kgrid,
    wannier_model.kpoints, Rvecs_ws)
HR_w = Wannier.fourier(kRvectors, Hk_coarse)

# === Step 3: Fine k-mesh interpolation using CrystalBonds ===
Rmat = Int.(Rvecs_ws.R)
N = Int.(Rvecs_ws.N)
nk = 30
kpoints = zeros(Float64, 3, nk^3)
let idx = 1
    for i in 0:(nk - 1), j in 0:(nk - 1), kk in 0:(nk - 1)
        kpoints[:, idx] = [i/nk, j/nk, kk/nk]
        idx += 1
    end
end

Hk = fourier_interpolate(HR_w, Rmat, kpoints; N = N)
evals, evecs = diag_Hk(Hk)
println(
    "Band range: $(round(minimum(evals), digits=2)) to $(round(maximum(evals), digits=2)) eV",
)

# === Step 4: WOHP using CrystalBonds ===
E_range = collect(range(minimum(evals) - 2, maximum(evals) + 2; length = 500))
wohp = compute(WOHP(), Hk, evals, evecs; E_range, method = GaussianSmearing(0.15))
woop = compute(WOOP(), Hk, evals, evecs; E_range, method = GaussianSmearing(0.15))

cc_wohp = extract_bond(wohp, 1:4, 5:8)
icohp = integrate(wohp, eF_eV)
println("IpCOHP at E_F: $(round(icohp, digits=4))")

# === Step 5: Plot ===
try
    using CairoMakie

    fig = Figure(size = (900, 500))
    ax1 = Axis(fig[1, 1]; xlabel = "DOS", ylabel = "E (eV)", title = "Diamond DOS")
    lines!(ax1, woop.total, E_range; color = :black, linewidth = 1.5)
    hlines!(ax1, [eF_eV]; color = :gray, linestyle = :dot)

    ax2 = Axis(fig[1, 2]; xlabel = "WOHP", ylabel = "E (eV)", title = "Diamond C-C WOHP")
    lines!(ax2, cc_wohp, E_range; color = :black, linewidth = 2)
    vlines!(ax2, [0]; color = :black, linestyle = :dash, linewidth = 0.5)
    hlines!(ax2, [eF_eV]; color = :gray, linestyle = :dot)
    linkyaxes!(ax1, ax2)

    outfile = joinpath(@__DIR__, "210_diamond_dftk_pbe.png")
    save(outfile, fig; px_per_unit = 2)
    println("Saved: $outfile")
catch e
    isa(e, ArgumentError) && println("CairoMakie not available, skipping plot.")
end
println("Done!")
