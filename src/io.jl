# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara

"""
    read_w90_hr(filename) -> (HR, Rvectors, N)

Read a Wannier90-compatible `_hr.dat` or `_HH_R.dat` file.

If a `.ndegen` file exists alongside, degeneracy weights are read from it.
Otherwise, all weights default to 1.

# Returns

  - `HR::Array{ComplexF64, 3}`: `(n_wan, n_wan, n_Rvecs)` real-space Hamiltonian
  - `Rvectors::Matrix{Int}`: `(3, n_Rvecs)` integer R-vector coordinates
  - `N::Vector{Int}`: `(n_Rvecs,)` degeneracy weights for each R-vector
"""
function read_w90_hr(filename::AbstractString)
    lines = readlines(filename)
    idx = 1

    # Line 1: header (skip)
    idx += 1

    # Line 2: n_wan
    n_wan = parse(Int, strip(lines[idx]))
    idx += 1

    # Line 3: n_Rvecs
    n_Rvecs = parse(Int, strip(lines[idx]))
    idx += 1

    # Allocate arrays
    HR = zeros(ComplexF64, n_wan, n_wan, n_Rvecs)
    Rvectors = zeros(Int, 3, n_Rvecs)

    # Read H(R) data: n_wan * n_wan entries per R-vector
    for iR in 1:n_Rvecs
        for j in 1:n_wan
            for i in 1:n_wan
                tokens = split(strip(lines[idx]))
                R1 = parse(Int, tokens[1])
                R2 = parse(Int, tokens[2])
                R3 = parse(Int, tokens[3])
                mu = parse(Int, tokens[4])
                nu = parse(Int, tokens[5])
                re = parse(Float64, tokens[6])
                im_val = parse(Float64, tokens[7])

                Rvectors[:, iR] = [R1, R2, R3]
                HR[mu, nu, iR] = complex(re, im_val)
                idx += 1
            end
        end
    end

    # Read degeneracy weights from .ndegen file if it exists
    ndegen_file = filename * ".ndegen"
    if isfile(ndegen_file)
        ndegen_text = read(ndegen_file, String)
        N = parse.(Int, split(strip(ndegen_text)))
        length(N) == n_Rvecs ||
            error("ndegen file has $(length(N)) entries, expected $(n_Rvecs)")
    else
        N = ones(Int, n_Rvecs)
    end

    return HR, Rvectors, N
end
