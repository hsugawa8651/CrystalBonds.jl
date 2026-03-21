# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara

"""
    fourier_interpolate(HR, Rvectors, kpoints; N=ones(n_Rvecs)) -> Hk

Fourier interpolate real-space Hamiltonian H(R) to k-space H(k).

    H(k) = Σ_R (1/N_R) exp(i 2π k·R) H(R)

where N_R is the degeneracy weight of each R-vector.

# Arguments

  - `HR::AbstractArray{<:Complex, 3}`: `(n_wan, n_wan, n_Rvecs)` real-space Hamiltonian
  - `Rvectors::AbstractMatrix{<:Integer}`: `(3, n_Rvecs)` R-vector coordinates
  - `kpoints::AbstractMatrix{<:Real}`: `(3, n_kpts)` fractional k-point coordinates

# Keyword arguments

  - `N::AbstractVector{<:Integer}`: `(n_Rvecs,)` degeneracy weights (default: all ones)

# Returns

  - `Hk::Array{ComplexF64, 3}`: `(n_wan, n_wan, n_kpts)` k-space Hamiltonian
"""
function fourier_interpolate(
    HR::AbstractArray{<:Complex, 3},
    Rvectors::AbstractMatrix{<:Integer},
    kpoints::AbstractMatrix{<:Real};
    N::AbstractVector{<:Integer} = ones(Int, size(HR, 3)),
)
    n_wan = size(HR, 1)
    n_Rvecs = size(HR, 3)
    n_kpts = size(kpoints, 2)

    Hk = zeros(ComplexF64, n_wan, n_wan, n_kpts)

    for ik in 1:n_kpts
        k = @view kpoints[:, ik]
        for iR in 1:n_Rvecs
            R = @view Rvectors[:, iR]
            phase = cispi(2 * dot(k, R)) / N[iR]
            @views Hk[:, :, ik] .+= phase .* HR[:, :, iR]
        end
    end

    return Hk
end

"""
    diag_Hk(Hk) -> (eigenvalues, eigenvectors)

Diagonalize k-space Hamiltonian at each k-point.

# Arguments

  - `Hk::AbstractArray{<:Complex, 3}`: `(n_wan, n_wan, n_kpts)`

# Returns

  - `eigenvalues::Matrix{Float64}`: `(n_wan, n_kpts)` sorted eigenvalues
  - `eigenvectors::Array{ComplexF64, 3}`: `(n_wan, n_wan, n_kpts)` eigenvectors as columns
"""
function diag_Hk(Hk::AbstractArray{<:Complex, 3})
    n_wan, _, n_kpts = size(Hk)

    eigenvalues = zeros(Float64, n_wan, n_kpts)
    eigenvectors = zeros(ComplexF64, n_wan, n_wan, n_kpts)

    for ik in 1:n_kpts
        F = eigen(Hermitian(Hk[:, :, ik]))
        eigenvalues[:, ik] = F.values
        eigenvectors[:, :, ik] = F.vectors
    end

    return eigenvalues, eigenvectors
end
