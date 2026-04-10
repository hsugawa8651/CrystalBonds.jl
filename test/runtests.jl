# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara

using CrystalBonds
using LinearAlgebra
using Test

@testset "CrystalBonds.jl" begin
    @testset "Type construction" begin
        @test WOHP() isa AbstractBondAnalysis
        @test WOOP() isa AbstractBondAnalysis
        @test GaussianSmearing(0.1) isa IntegrationMethod
        @test GaussianSmearing().sigma == 0.15  # default
    end

    @testset "Error fallbacks" begin
        Hk = zeros(ComplexF64, 2, 2, 1)
        # compute fallback fires for unimplemented AbstractBondAnalysis subtypes
        @test_throws ErrorException compute(WOHP(), Hk)
        @test_throws ErrorException compute(WOOP(), Hk)
    end

    @testset "read_w90_hr" begin
        fixture = joinpath(@__DIR__, "fixtures", "diamond_hr.dat")
        HR, Rvectors, N = read_w90_hr(fixture)

        n_wan = size(HR, 1)
        n_Rvecs = size(HR, 3)

        # Dimensions
        @test n_wan > 0
        @test n_Rvecs > 0
        @test size(HR) == (n_wan, n_wan, n_Rvecs)
        @test size(Rvectors) == (3, n_Rvecs)
        @test length(N) == n_Rvecs
        @test n_wan == 8       # Diamond sp3: 2 atoms x 4 orbitals
        @test n_Rvecs == 617   # Wigner-Seitz R-vectors for 8x8x8 kgrid

        # Degeneracy weights are positive integers
        @test all(N .>= 1)

        # R=0 must exist
        R0_idx = findfirst(i -> all(Rvectors[:, i] .== 0), 1:n_Rvecs)
        @test R0_idx !== nothing

        # H(R=0) should be approximately Hermitian
        H0 = HR[:, :, R0_idx]
        @test H0 ≈ H0' atol=1e-10
    end

    @testset "fourier_interpolate + diag_Hk" begin
        fixture = joinpath(@__DIR__, "fixtures", "diamond_hr.dat")
        HR, Rvectors, N = read_w90_hr(fixture)
        n_wan = size(HR, 1)

        # Small k-mesh for testing
        kpoints = Float64[0 0.5 0.25; 0 0.5 0.25; 0 0.5 0.25]  # Gamma, X, arbitrary
        Hk = fourier_interpolate(HR, Rvectors, kpoints; N = N)
        n_kpts = size(kpoints, 2)

        # Dimensions
        @test size(Hk) == (n_wan, n_wan, n_kpts)

        # H(k) is Hermitian at every k-point
        for ik in 1:n_kpts
            @test Hk[:, :, ik] ≈ Hk[:, :, ik]' atol=1e-10
        end

        # Diagonalize
        evals, evecs = diag_Hk(Hk)
        @test size(evals) == (n_wan, n_kpts)
        @test size(evecs) == (n_wan, n_wan, n_kpts)

        # Eigenvalues are real (guaranteed by Hermitian, but check sorted)
        for ik in 1:n_kpts
            @test issorted(evals[:, ik])
        end

        # Eigenvectors are orthonormal: V'V ≈ I
        for ik in 1:n_kpts
            V = evecs[:, :, ik]
            @test V' * V ≈ I(n_wan) atol=1e-10
        end

        # Self-consistency: [C†HC]_{jj} = ε_j
        for ik in 1:n_kpts
            V = evecs[:, :, ik]
            H = Hk[:, :, ik]
            reconstructed = real(diag(V' * H * V))
            @test reconstructed ≈ evals[:, ik] atol=1e-10
        end
    end

    @testset "gaussian_delta" begin
        sigma = 0.1

        # Non-negative for all inputs
        for eps in -10.0:0.5:10.0, E in -10.0:0.5:10.0
            @test gaussian_delta(eps, E, sigma) >= 0.0
        end

        # Peak at eps = E
        E0 = 5.0
        peak = gaussian_delta(E0, E0, sigma)
        @test peak > 0.0
        @test gaussian_delta(E0 + 0.5, E0, sigma) < peak
        @test gaussian_delta(E0 - 0.5, E0, sigma) < peak

        # Normalization: integral over E ≈ 1
        dE = 0.001
        E_grid = -5.0:dE:15.0
        integral = sum(gaussian_delta(5.0, E, sigma) for E in E_grid) * dE
        @test integral ≈ 1.0 atol=0.01

        # Cutoff: returns 0 beyond 6 sigma
        @test gaussian_delta(0.0, 10.0, sigma) == 0.0
    end

    @testset "compute WOHP/WOOP (Diamond)" begin
        fixture = joinpath(@__DIR__, "fixtures", "diamond_hr.dat")
        HR, Rvectors, N = read_w90_hr(fixture)
        n_wan = size(HR, 1)

        # Build k-mesh (small for testing)
        nk = 10
        kpoints = zeros(Float64, 3, nk^3)
        idx = 1
        for i in 0:(nk - 1), j in 0:(nk - 1), kk in 0:(nk - 1)
            kpoints[:, idx] = [i/nk, j/nk, kk/nk]
            idx += 1
        end

        Hk = fourier_interpolate(HR, Rvectors, kpoints; N = N)
        evals, evecs = diag_Hk(Hk)

        E_min = minimum(evals) - 2.0
        E_max = maximum(evals) + 2.0
        E_range = range(E_min, E_max; length = 200)
        sigma = 0.3

        # WOOP (= DOS)
        woop = compute(WOOP(), Hk, evals, evecs;
            E_range = collect(E_range), method = GaussianSmearing(sigma))

        # DOS >= 0 everywhere
        @test all(woop.total .>= -1e-15)

        # DOS integral = n_wannier (particle conservation)
        dE = woop.E_range[2] - woop.E_range[1]
        dos_integral = sum(woop.total) * dE
        @test dos_integral ≈ n_wan atol=0.1

        # Result dimensions
        @test size(woop.matrix) == (n_wan, n_wan, length(E_range))
        @test length(woop.total) == length(E_range)

        # WOHP
        wohp = compute(WOHP(), Hk, evals, evecs;
            E_range = collect(E_range), method = GaussianSmearing(sigma))

        @test size(wohp.matrix) == (n_wan, n_wan, length(E_range))

        # Self-consistency: total WOHP = -Σ ε_j δ(ε_j - E)
        # (NOT -E × DOS; that only holds in the delta function limit)
        # Verify by computing -Σ ε_j δ directly
        eps_dos = zeros(length(E_range))
        for ik in 1:size(evals, 2), j in 1:n_wan
            for iE in eachindex(E_range)
                d = gaussian_delta(evals[j, ik], E_range[iE], sigma)
                d == 0.0 && continue
                eps_dos[iE] += evals[j, ik] * d
            end
        end
        eps_dos ./= size(evals, 2)
        neg_eps_dos = -eps_dos
        nonzero = abs.(neg_eps_dos) .> 1e-10
        if any(nonzero)
            max_abs_err = maximum(abs.(wohp.total[nonzero] - neg_eps_dos[nonzero]))
            @test max_abs_err < 1e-10
        end

        # WOHP: C-C bond (orbitals 1:4 <-> 5:8) should be net bonding below gap
        cc_wohp = dropdims(sum(wohp.matrix[1:4, 5:8, :]; dims = (1, 2)); dims = (1, 2))
        gap_E = (maximum(evals[4, :]) + minimum(evals[5, :])) / 2  # mid-gap
        below_gap = collect(E_range) .< gap_E
        bonding_sum = sum(cc_wohp[below_gap]) * dE
        @test bonding_sum > 0  # net bonding below gap
    end

    @testset "extract_bond + integrate" begin
        fixture = joinpath(@__DIR__, "fixtures", "diamond_hr.dat")
        HR, Rvectors, N = read_w90_hr(fixture)
        n_wan = size(HR, 1)

        nk = 10
        kpoints = zeros(Float64, 3, nk^3)
        let idx = 1
            for i in 0:(nk - 1), j in 0:(nk - 1), kk in 0:(nk - 1)
                kpoints[:, idx] = [i/nk, j/nk, kk/nk]
                idx += 1
            end
        end

        Hk = fourier_interpolate(HR, Rvectors, kpoints; N = N)
        evals, evecs = diag_Hk(Hk)
        E_range = collect(range(minimum(evals)-2, maximum(evals)+2; length = 200))
        wohp = compute(WOHP(), Hk, evals, evecs;
            E_range = E_range, method = GaussianSmearing(0.3))

        # extract_bond
        cc = extract_bond(wohp, 1:4, 5:8)
        @test length(cc) == length(E_range)
        @test cc isa Vector{Float64}

        # integrate returns finite scalar
        gap_E = (maximum(evals[4, :]) + minimum(evals[5, :])) / 2
        icohp = integrate(wohp, gap_E)
        @test isfinite(icohp)

        # IpCOHP at energy below all bands should be ~0
        icohp_low = integrate(wohp, minimum(evals) - 10.0)
        @test abs(icohp_low) < 1e-10
    end

    @testset "read_w90_hr (GaAs)" begin
        fixture = joinpath(@__DIR__, "fixtures", "gaas_hr.dat")
        HR, Rvectors, N = read_w90_hr(fixture)

        n_wan = size(HR, 1)
        n_Rvecs = size(HR, 3)

        @test n_wan == 8        # GaAs sp3: 2 atoms x 4 orbitals
        @test n_Rvecs == 617
        @test size(HR) == (n_wan, n_wan, n_Rvecs)
        @test size(Rvectors) == (3, n_Rvecs)
        @test length(N) == n_Rvecs
        @test all(N .>= 1)

        # R=0 must exist
        R0_idx = findfirst(i -> all(Rvectors[:, i] .== 0), 1:n_Rvecs)
        @test R0_idx !== nothing

        # H(R=0) should be approximately Hermitian
        H0 = HR[:, :, R0_idx]
        @test H0 ≈ H0' atol=1e-10
    end

    @testset "compute WOHP/WOOP (GaAs)" begin
        fixture = joinpath(@__DIR__, "fixtures", "gaas_hr.dat")
        HR, Rvectors, N = read_w90_hr(fixture)
        n_wan = size(HR, 1)

        nk = 10
        kpoints = zeros(Float64, 3, nk^3)
        idx = 1
        for i in 0:(nk - 1), j in 0:(nk - 1), kk in 0:(nk - 1)
            kpoints[:, idx] = [i/nk, j/nk, kk/nk]
            idx += 1
        end

        Hk = fourier_interpolate(HR, Rvectors, kpoints; N = N)
        evals, evecs = diag_Hk(Hk)

        E_range = collect(range(minimum(evals) - 2, maximum(evals) + 2; length = 200))
        sigma = 0.3

        # WOOP (= DOS)
        woop = compute(WOOP(), Hk, evals, evecs;
            E_range = E_range, method = GaussianSmearing(sigma))

        # DOS >= 0 everywhere
        @test all(woop.total .>= -1e-15)

        # DOS integral = n_wannier (particle conservation)
        dE = E_range[2] - E_range[1]
        dos_integral = sum(woop.total) * dE
        @test dos_integral ≈ n_wan atol=0.1

        # WOHP
        wohp = compute(WOHP(), Hk, evals, evecs;
            E_range = E_range, method = GaussianSmearing(sigma))

        @test size(wohp.matrix) == (n_wan, n_wan, length(E_range))

        # Ga-As bond (orbitals 1:4 = Ga, 5:8 = As) should be net bonding below gap
        gaas_wohp = extract_bond(wohp, 1:4, 5:8)
        gap_E = (maximum(evals[4, :]) + minimum(evals[5, :])) / 2
        below_gap = E_range .< gap_E
        bonding_sum = sum(gaas_wohp[below_gap]) * dE
        @test bonding_sum > 0  # net bonding below gap

        # IpCOHP at mid-gap should be positive and finite
        icohp = integrate(wohp, gap_E)
        @test isfinite(icohp)
        @test icohp > 0
    end

    @testset "integrate boundary conditions (Diamond)" begin
        fixture = joinpath(@__DIR__, "fixtures", "diamond_hr.dat")
        HR, Rvectors, N = read_w90_hr(fixture)
        n_wan = size(HR, 1)

        nk = 10
        kpoints = zeros(Float64, 3, nk^3)
        let idx = 1
            for i in 0:(nk - 1), j in 0:(nk - 1), kk in 0:(nk - 1)
                kpoints[:, idx] = [i/nk, j/nk, kk/nk]
                idx += 1
            end
        end

        Hk = fourier_interpolate(HR, Rvectors, kpoints; N = N)
        evals, evecs = diag_Hk(Hk)
        E_range = collect(range(minimum(evals) - 2, maximum(evals) + 2; length = 200))
        wohp = compute(WOHP(), Hk, evals, evecs;
            E_range = E_range, method = GaussianSmearing(0.3))

        gap_E = (maximum(evals[4, :]) + minimum(evals[5, :])) / 2

        # IpCOHP below all bands should be ~0
        @test abs(integrate(wohp, minimum(evals) - 10.0)) < 1e-10

        # IpCOHP at mid-gap should be finite
        @test isfinite(integrate(wohp, gap_E))

        # IpCOHP above all bands should be finite
        @test isfinite(integrate(wohp, maximum(evals) + 10.0))

        # Bond-resolved IpCOHP (C-C) at mid-gap should be positive (net bonding)
        cc_wohp = extract_bond(wohp, 1:4, 5:8)
        dE = E_range[2] - E_range[1]
        gap_idx = findlast(E_range .<= gap_E)
        bond_icohp = sum(cc_wohp[1:gap_idx]) * dE
        @test bond_icohp > 0
    end
end
