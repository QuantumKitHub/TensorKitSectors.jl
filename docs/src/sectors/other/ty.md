# The $\mathbb Z_N$-Tambara-Yamagami categories: `TambaraYamagami`

`TambaraYamagami{N, K}` represents the Tambara-Yamagami fusion category ${\rm TY}(\mathbb Z_N, K)$ based on the cyclic group of order $N$.

The $N+1$ simple objects coincide with the group elements of $\mathbb Z_N$ supplemented with one non-invertible object $m$. The type parameter `K` specifies the Frobenius-Schur indicator of $m$ as $\varkappa_m = (-1)^K = \pm 1$, i.e. `K = false` corresponds to $\varkappa_m = 1$ and `K = true` to $\varkappa_m = -1$.

## Sector type

```@docs; canonical = false
TambaraYamagami
```
Here, the type parameters `N` and `K` correspond respectively to the order of the underlying cyclic group and the Frobenius-Schur sign $(-1)^K$ of the non-invertible object.

## Fusion Rules

The non-trivial fusion rules read

```math
g ⊗ h = g + h \mod N, \qquad g ⊗ m = m ⊗ g = m, \qquad m ⊗ m = \bigoplus_{g ∈ ℤ_N} g,
```
for all group elements $g, h \in \mathbb Z_N$.

Hence `FusionStyle(::Type{<:TambaraYamagami}) = SimpleFusion()`.

The quantum dimensions are

```math
d_g = 1, \quad \forall g\in\mathbb Z_N, \quad {\rm and} \quad d_m = \sqrt{N}.
```

## Topological Data

We write $χ(g, h) = \exp(2π i g h / N)$ for a normalised non-degenerate symmetric bicharacter on $\mathbb Z_N$. The nontrivial F-symbols are then given by

```math
F^{g\,m\,h}_{m} = χ(g, h), \qquad
F^{m\,g\,m}_{h} = χ(g, h), \qquad
\left[F^{m\,m\,m}_{m}\right]_g^h = \frac{\varkappa_m}{\sqrt{N}}\,\overline{χ(g, h)},
```

for all $g, h \in \mathbb Z_N$.

Crucially, there exists no braiding on this fusion category, except when $N=2$, in which case it coincides with `IsingAnyon`.

## Iteration and basis conventions
`values(TambaraYamagami{N,K})` iterates the labels `0, 1, …, N-1,:m` in increasing order.

## References
[1] D. Tambara and S. Yamagami, *Tensor categories with fusion rules of self-duality for
    finite abelian groups*, J. Algebra **209**, 692-707 (1998).

[2] M. Barkeshli, P. Bonderson, M. Cheng and Z. Wang, *Symmetry Fractionalization, Defects, and Gauging of Topological Phases*, Phys. Rev. B **100**, 115147 (2019), [arXiv:1410.4540](https://arxiv.org/abs/1410.4540).