# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Hiroharu Sugawara

"""
    gaussian_delta(eps, E, sigma)

Gaussian approximation to the Dirac delta function.

    δ(ε - E) ≈ exp(-(ε-E)²/(2σ²)) / (σ√(2π))

Returns 0.0 when |ε - E| > 6σ for efficiency.
"""
function gaussian_delta(eps::Real, E::Real, sigma::Real)
    x = (eps - E) / sigma
    abs(x) > 6 && return 0.0
    return exp(-x^2 / 2) / (sigma * sqrt(2π))
end
