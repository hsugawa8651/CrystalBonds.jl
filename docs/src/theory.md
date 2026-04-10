# Theory

## Wannier Orbital Hamilton Population (WOHP)

The WOHP provides energy-resolved bonding analysis using Wannier functions
as the local basis. For orbital pair ``(\mu, \nu)`` at energy ``E``:

```math
\text{WOHP}_{\mu\nu}(E) = -\frac{1}{N_k} \sum_{k,j}
\text{Re}\!\left[ H_{\mu\nu}(\mathbf{k})\, C^*_{\mu j}(\mathbf{k})\, C_{\nu j}(\mathbf{k}) \right]
\tilde{\delta}(\varepsilon_j(\mathbf{k}) - E)
```

where:
- ``H_{\mu\nu}(\mathbf{k})`` is the Hamiltonian matrix element in the Wannier basis
- ``C_{\mu j}(\mathbf{k})`` are eigenvector components
- ``\varepsilon_j(\mathbf{k})`` are eigenvalues
- ``\tilde{\delta}`` is a broadened delta function (e.g., Gaussian smearing)

## Wannier Orbital Overlap Population (WOOP)

The WOOP is the same as WOHP but without the Hamiltonian weight:

```math
\text{WOOP}_{\mu\nu}(E) = \frac{1}{N_k} \sum_{k,j}
\text{Re}\!\left[ C^*_{\mu j}(\mathbf{k})\, C_{\nu j}(\mathbf{k}) \right]
\tilde{\delta}(\varepsilon_j(\mathbf{k}) - E)
```

The sum over all orbital pairs gives the total density of states (DOS):

```math
\text{DOS}(E) = \sum_{\mu} \text{WOOP}_{\mu\mu}(E)
```

## Integrated WOHP (IpCOHP)

The integrated WOHP up to the Fermi energy quantifies
net bond strength:

```math
\text{IpCOHP} = \int_{-\infty}^{E_F} \text{WOHP}(E)\, dE
```

Positive IpCOHP indicates net bonding; negative indicates net antibonding.

## Sign convention

CrystalBonds.jl follows the Kundu convention:

| Sign | Meaning | Equivalent in LOBSTER |
| ---- | ------- | --------------------- |
| WOHP > 0 | bonding | -pCOHP > 0 |
| WOHP < 0 | antibonding | -pCOHP < 0 |

The relationship to the Deringer/LOBSTER convention is:

```math
\text{WOHP} = -\text{pCOHP}
```

When comparing with LOBSTER output, note that LOBSTER plots `-pCOHP`
(positive = bonding), which matches WOHP directly.

## Fourier interpolation

The real-space Hamiltonian ``H(\mathbf{R})`` from Wannier90 is interpolated
to an arbitrary k-point:

```math
H_{\mu\nu}(\mathbf{k}) = \sum_{\mathbf{R}} \frac{1}{N_{\mathbf{R}}}
e^{i 2\pi \mathbf{k} \cdot \mathbf{R}}\, H_{\mu\nu}(\mathbf{R})
```

where ``N_{\mathbf{R}}`` are the Wigner-Seitz degeneracy weights.

## Gaussian smearing

The delta function is approximated by a Gaussian:

```math
\tilde{\delta}(\varepsilon - E) = \frac{1}{\sigma\sqrt{2\pi}}
\exp\!\left( -\frac{(\varepsilon - E)^2}{2\sigma^2} \right)
```

with a cutoff at ``|{\varepsilon - E}| > 6\sigma`` for efficiency.

## Self-consistency identities

For Gaussian smearing, the following identity holds exactly:

```math
\text{total WOHP}(E) = -\frac{1}{N_k} \sum_{k,j} \varepsilon_j(\mathbf{k})\,
\tilde{\delta}(\varepsilon_j(\mathbf{k}) - E)
```

!!! note
    The simpler identity ``\text{total WOHP}(E) = -E \times \text{DOS}(E)``
    only holds in the exact delta function limit, not with Gaussian smearing.

## References

- Kundu et al., "Population Analysis with Wannier Orbitals",
  J. Chem. Phys. (2020). DOI: [10.1063/5.0032605](https://doi.org/10.1063/5.0032605)
- Deringer et al., "Crystal Orbital Hamilton Population (COHP) Analysis
  As Projected from Plane-Wave Basis Sets",
  J. Phys. Chem. A 115, 5461 (2011). DOI: [10.1021/jp202489s](https://doi.org/10.1021/jp202489s)
- Dronskowski and Bloechl, "Crystal Orbital Hamilton Populations (COHP).
  Energy-Resolved Visualization of Chemical Bonding in Solids Based on
  Density-Functional Calculations",
  J. Phys. Chem. 97, 8617 (1993). DOI: [10.1021/j100135a014](https://doi.org/10.1021/j100135a014)
