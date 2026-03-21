# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara
#
# GaAs WOHP from Wannier90 _hr.dat file
# No external dependencies beyond CrystalBonds.jl
#
# Usage: julia --project=.. examples/120_gaas_w90.jl

using CrystalBonds
using LinearAlgebra

# --- Read Wannier90 _hr.dat ---
hr_file = joinpath(@__DIR__, "..", "test", "fixtures", "gaas_hr.dat")
HR, Rvectors, N = read_w90_hr(hr_file)
n_wan = size(HR, 1)
println("Loaded: n_wan=$n_wan, n_Rvecs=$(size(HR,3))")

# --- Build fine k-mesh and interpolate ---
nk = 30
kpoints = zeros(Float64, 3, nk^3)
let idx = 1
    for i in 0:(nk - 1), j in 0:(nk - 1), kk in 0:(nk - 1)
        kpoints[:, idx] = [i/nk, j/nk, kk/nk]
        idx += 1
    end
end
println("k-mesh: $(nk)^3 = $(nk^3) points")

Hk = fourier_interpolate(HR, Rvectors, kpoints; N = N)
evals, evecs = diag_Hk(Hk)
println(
    "Band range: $(round(minimum(evals), digits=2)) to $(round(maximum(evals), digits=2)) eV",
)

# --- Compute WOHP ---
E_range = collect(range(minimum(evals) - 2, maximum(evals) + 2; length = 500))

wohp = compute(WOHP(), Hk, evals, evecs;
    E_range = E_range, method = GaussianSmearing(0.15))

woop = compute(WOOP(), Hk, evals, evecs;
    E_range = E_range, method = GaussianSmearing(0.15))

# --- Extract Ga-As bond (orbitals 1:4 = Ga, 5:8 = As) ---
bond_wohp = extract_bond(wohp, 1:4, 5:8)

# --- Orbital decomposition ---
wohp_ss = wohp.matrix[1, 5, :]
wohp_sp = dropdims(sum(wohp.matrix[1:1, 6:8, :]; dims = (1, 2)); dims = (1, 2))
wohp_ps = dropdims(sum(wohp.matrix[2:4, 5:5, :]; dims = (1, 2)); dims = (1, 2))
wohp_pp = dropdims(sum(wohp.matrix[2:4, 6:8, :]; dims = (1, 2)); dims = (1, 2))

# --- IpCOHP ---
gap_E = (maximum(evals[4, :]) + minimum(evals[5, :])) / 2
icohp = integrate(wohp, gap_E)
println("IpCOHP at mid-gap ($(round(gap_E, digits=2)) eV): $(round(icohp, digits=4))")

# --- DOS check ---
dE = E_range[2] - E_range[1]
dos_integral = sum(woop.total) * dE
println("DOS integral: $(round(dos_integral, digits=4)) (expected: $n_wan)")

# --- Plot ---
try
    using CairoMakie

    fig = Figure(size = (900, 500))

    ax1 = Axis(fig[1, 1];
        xlabel = "DOS (states/eV)", ylabel = "E (eV)", title = "GaAs DOS")
    lines!(ax1, woop.total, E_range; color = :black, linewidth = 1.5)
    hlines!(ax1, [0]; color = :gray, linestyle = :dot, linewidth = 0.5)

    ax2 = Axis(fig[1, 2];
        xlabel = "WOHP", ylabel = "E (eV)", title = "GaAs Ga-As WOHP")
    lines!(ax2, bond_wohp, E_range; color = :black, linewidth = 2, label = "total")
    lines!(ax2, wohp_ss, E_range; color = :red, linewidth = 1, label = "s-s")
    lines!(ax2, wohp_sp, E_range; color = :green, linewidth = 1, label = "s-p")
    lines!(ax2, wohp_ps, E_range; color = :orange, linewidth = 1, label = "p-s")
    lines!(ax2, wohp_pp, E_range; color = :blue, linewidth = 1, label = "p-p")
    vlines!(ax2, [0]; color = :black, linestyle = :dash, linewidth = 0.5)
    hlines!(ax2, [0]; color = :gray, linestyle = :dot, linewidth = 0.5)
    axislegend(ax2; position = :lb, labelsize = 10)

    linkyaxes!(ax1, ax2)

    outfile = joinpath(@__DIR__, "120_gaas_w90.png")
    save(outfile, fig; px_per_unit = 2)
    println("Saved: $outfile")
catch e
    if isa(e, ArgumentError)
        println("CairoMakie not available, skipping plot.")
    else
        rethrow(e)
    end
end

println("Done!")
