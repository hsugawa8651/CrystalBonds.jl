# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara

"""
    extract_bond(result, orb1, orb2) -> Vector{Float64}

Extract bond-resolved WOHP/WOOP between two sets of orbitals.

# Arguments

  - `result::BondAnalysisResult`: output from `compute`
  - `orb1`: orbital indices for atom 1 (e.g., `1:4`)
  - `orb2`: orbital indices for atom 2 (e.g., `5:8`)

# Returns

  - Energy-resolved bond WOHP/WOOP vector of length `n_E`
"""
function extract_bond(result::BondAnalysisResult, orb1, orb2)
    return dropdims(sum(result.matrix[orb1, orb2, :]; dims = (1, 2)); dims = (1, 2))
end

"""
    integrate(result, E_fermi) -> Float64

Compute integrated WOHP (IpCOHP) up to the Fermi energy.

    IpCOHP = ∫_{-∞}^{E_F} WOHP(E) dE

# Arguments

  - `result::BondAnalysisResult`: output from `compute`
  - `E_fermi::Real`: Fermi energy (eV)

# Returns

  - Integrated value (scalar). Positive = net bonding, negative = net antibonding.
"""
function integrate(result::BondAnalysisResult, E_fermi::Real)
    E = result.E_range
    dE = length(E) > 1 ? E[2] - E[1] : 0.0
    below = E .<= E_fermi
    return sum(result.total[below]) * dE
end
