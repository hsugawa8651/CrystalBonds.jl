# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara

"""
    compute(::WOHP, Hk, evals, evecs; E_range, method) -> BondAnalysisResult

Compute Wannier Orbital Hamilton Population (WOHP) in k-space.

    WOHP_{μν}(E) = -(1/Nk) Σ_{k,j} Re[H_{μν}(k) C*_{μj}(k) C_{νj}(k)] δ(ε_j(k) - E)

Positive WOHP = bonding, negative = antibonding (Kundu convention).

# Arguments

  - `Hk::AbstractArray{<:Complex, 3}`: `(n_wan, n_wan, n_kpts)` k-space Hamiltonian
  - `evals::AbstractMatrix{<:Real}`: `(n_wan, n_kpts)` eigenvalues
  - `evecs::AbstractArray{<:Complex, 3}`: `(n_wan, n_wan, n_kpts)` eigenvectors

# Keyword arguments

  - `E_range::AbstractVector{<:Real}`: energy grid (eV)
  - `method::IntegrationMethod`: integration method (default: `GaussianSmearing()`)
"""
function compute(
    ::WOHP,
    Hk::AbstractArray{<:Complex, 3},
    evals::AbstractMatrix{<:Real},
    evecs::AbstractArray{<:Complex, 3};
    E_range::AbstractVector{<:Real},
    method::IntegrationMethod = GaussianSmearing(),
)
    n_wan, n_kpts = size(evals)
    n_E = length(E_range)
    sigma = method.sigma

    wohp_matrix = zeros(Float64, n_wan, n_wan, n_E)

    for ik in 1:n_kpts
        H_k = @view Hk[:, :, ik]
        C = @view evecs[:, :, ik]
        eps_k = @view evals[:, ik]

        for j in 1:n_wan
            for iE in eachindex(E_range)
                d = gaussian_delta(eps_k[j], E_range[iE], sigma)
                d == 0.0 && continue
                for nu in 1:n_wan, mu in 1:n_wan
                    wohp_matrix[mu, nu, iE] +=
                        -real(H_k[mu, nu] *
                              conj(C[mu, j]) * C[nu, j]) * d
                end
            end
        end
    end

    wohp_matrix ./= n_kpts
    wohp_total = dropdims(sum(wohp_matrix; dims = (1, 2)); dims = (1, 2))
    return BondAnalysisResult(wohp_matrix, wohp_total, collect(Float64, E_range))
end

"""
    compute(::WOOP, Hk, evals, evecs; E_range, method) -> BondAnalysisResult

Compute Wannier Orbital Overlap Population (WOOP), equivalent to projected DOS.
Same as WOHP but without the Hamiltonian weight (weight = identity).
Sum over all orbitals gives total DOS.
"""
function compute(
    ::WOOP,
    Hk::AbstractArray{<:Complex, 3},
    evals::AbstractMatrix{<:Real},
    evecs::AbstractArray{<:Complex, 3};
    E_range::AbstractVector{<:Real},
    method::IntegrationMethod = GaussianSmearing(),
)
    n_wan, n_kpts = size(evals)
    n_E = length(E_range)
    sigma = method.sigma

    woop_matrix = zeros(Float64, n_wan, n_wan, n_E)

    for ik in 1:n_kpts
        C = @view evecs[:, :, ik]
        eps_k = @view evals[:, ik]

        for j in 1:n_wan
            for iE in eachindex(E_range)
                d = gaussian_delta(eps_k[j], E_range[iE], sigma)
                d == 0.0 && continue
                for nu in 1:n_wan, mu in 1:n_wan
                    woop_matrix[mu, nu, iE] += real(conj(C[mu, j]) * C[nu, j]) * d
                end
            end
        end
    end

    woop_matrix ./= n_kpts
    woop_total = dropdims(sum(woop_matrix; dims = (1, 2)); dims = (1, 2))
    return BondAnalysisResult(woop_matrix, woop_total, collect(Float64, E_range))
end
