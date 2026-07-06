# ``CU(1)`` Representations: `CU1Irrep`

`CU1Irrep` represents irreducible representations of ``U(1) ⋊ C``, where `C` acts by charge conjugation.
This group is also known as ``O(2)``.

## Sector type

```@docs; canonical = false
CU1Irrep
```

Labels are pairs `(j, s)`, where `j` is a `U1Irrep` charge and `s` is an integer in `0:2`.
For `j == 0`, there are two one-dimensional irreps: `CU1Irrep(0, 0)` and `CU1Irrep(0, 1)`.
For `j > 0`, only `s == 2` is valid and the representation is two-dimensional.
The unit is `CU1Irrep(0, 0)`, and every irrep is self-dual.

## Fusion Rules

The two one-dimensional irreps fuse by XOR of the `s` label.
Fusing either one-dimensional irrep with a two-dimensional irrep returns the same two-dimensional label.
For two positive charges, the outputs are governed by sum and absolute difference:

```math
j_1 \otimes j_2 =
|j_1-j_2| \oplus (j_1+j_2),
```

with the special case `j1 == j2`, where the zero-charge difference splits into the two one-dimensional irreps.
The category has `FusionStyle(CU1Irrep) = SimpleFusion()`.

Quantum dimensions are

```math
d_{(0,0)} = d_{(0,1)} = 1,\qquad d_{(j,2)} = 2\quad (j>0).
```

## Topological Data

`CU1Irrep` is a bosonic representation category, so all twists are `1`.
The [`Fsymbol`](@ref) values are real and encode the chosen Clebsch-Gordan basis.
The [`Rsymbol`](@ref) is real and equals the N-symbol, with an additional negative sign on channels with `c.s == 1` and positive left overbraiding input charge.

## Fusion Tensor and Basis Conventions

`fusiontensor(a, b, c)` returns a rank-4 array of size ``d_a \times d_b \times d_c \times N_c^{ab}``.
Since `CU1Irrep` has simple fusion, the final multiplicity axis has length `0` or `1`.
For allowed channels, the tensor entries are real Clebsch-Gordan coefficients in the convention explained below.

For a two-dimensional irrep `(j, 2)` with `j > 0`, the basis is ordered as a pair of charge-conjugate ``U(1)`` weights, which we can denote by

```math
\ket{+j},\ \ket{-j}.
```

The zero-charge irreps `(0, 0)` and `(0, 1)` are one-dimensional. The label `(0, 0)` is even under charge conjugation, while `(0, 1)` is odd.

When two equal positive-charge irreps fuse to a zero-charge irrep, the fusion tensors pick the symmetric and antisymmetric combinations:

```math
\ket{(0,0)} =
\frac{1}{\sqrt{2}}\left(\ket{+j}\otimes\ket{-j} + \ket{-j}\otimes\ket{+j}\right),
```

```math
\ket{(0,1)} =
\frac{1}{\sqrt{2}}\left(\ket{+j}\otimes\ket{-j} - \ket{-j}\otimes\ket{+j}\right).
```

In array form, these are the entries

```math
C_{1,2,1} = \frac{1}{\sqrt{2}},\qquad
C_{2,1,1} = \pm\frac{1}{\sqrt{2}},
```

with the plus sign for `(0, 0)` and the minus sign for `(0, 1)`.

Fusing a zero-charge irrep with a two-dimensional irrep leaves the charge label unchanged.
The odd zero-charge irrep contributes a sign on the second basis vector:

```math
(0,s) \otimes (j,2) \to (j,2):\qquad
\ket{+j}\mapsto\ket{+j},\quad
\ket{-j}\mapsto (-1)^s\ket{-j}.
```

The same convention is used for `(j,2) ⊗ (0,s)`, with the sign attached to the second basis vector of the two-dimensional input.

For two positive charges, the sum channel is diagonal:

```math
(j_a,2)\otimes(j_b,2)\to(j_a+j_b,2):\qquad
\ket{+j_a,+j_b}\mapsto\ket{+(j_a+j_b)},\quad
\ket{-j_a,-j_b}\mapsto\ket{-(j_a+j_b)}.
```

The difference channel pairs opposite weights. If `j_a > j_b`,

```math
\ket{+j_a,-j_b}\mapsto\ket{+(j_a-j_b)},\qquad
\ket{-j_a,+j_b}\mapsto\ket{-(j_a-j_b)}.
```

If `j_b > j_a`, the output basis is ordered by the positive charge `j_b - j_a`, so the two nonzero entries are swapped accordingly:

```math
\ket{-j_a,+j_b}\mapsto\ket{+(j_b-j_a)},\qquad
\ket{+j_a,-j_b}\mapsto\ket{-(j_b-j_a)}.
```

All omitted entries are zero. These conventions determine the real [`Fsymbol`](@ref) values.

## Iteration

`values(CU1Irrep)` is infinite.
It starts with the two zero-charge irreps and then lists positive half-integer charges:

```jldoctest
julia> using TensorKitSectors

julia> values(CU1Irrep)[1]
Irrep[CU₁](0, 0)

julia> values(CU1Irrep)[2]
Irrep[CU₁](0, 1)

julia> values(CU1Irrep)[5]
Irrep[CU₁](3/2, 2)
```

## References

- [Orthogonal group in two dimensions](https://en.wikipedia.org/wiki/Orthogonal_group_in_two_dimensions)
- [Semidirect product](https://en.wikipedia.org/wiki/Semidirect_product)