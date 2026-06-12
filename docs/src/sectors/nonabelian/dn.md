# Dihedral Group Representations: `DNIrrep`

`DNIrrep{N}` represents irreducible representations of the dihedral group ``D_N = \mathbb{Z}_N ⋊ C``.
The aliases `D3Irrep` and `D4Irrep` are provided for `DNIrrep{3}` and `DNIrrep{4}`.

## Sector type

```@docs; canonical = false
DNIrrep
```

Labels are constructed as `DNIrrep{N}(j, isodd=false)`.
The integer `j` satisfies `0 <= j <= N ÷ 2`.
The `isodd` flag is valid only for one-dimensional irreps: always at `j == 0`, and also at `j == N/2` when `N` is even.

For odd `N`, the labels are

```math
(0,\text{false}),\ (0,\text{true}),\ 1,2,\ldots,(\frac{N-1}{2},\text{false}).
```

For even `N`, the labels are

```math
(0,\text{false}),\ (0,\text{true}),\ 1,2,\ldots,\frac{N}{2}-1,\ (N/2,\text{false}),\ (N/2,\text{true}).
```

The unit is `DNIrrep{N}(0, false)`, and all irreps are self-dual.

## Fusion Rules

The one-dimensional irreps fuse by XOR of the parity label.
The fusion product of a one-dimensional irrep with a two-dimensional irrep preserves the two-dimensional label, up to the special one-dimensional labels at `j == N/2` for even `N`.

For ordinary two-dimensional labels, fusion follows the dihedral character rule:

```math
\rho_i \otimes \rho_j = \rho_{|i-j|} \oplus \rho_{i+j},
```

where labels are folded back into the range `0:N÷2`.
When a folded output lands on a one-dimensional point (`0`, or `N/2` for even `N`), it splits into the corresponding even and odd one-dimensional irreps.

For `N < 3`, `FusionStyle(DNIrrep{N}) = UniqueFusion()`.
Otherwise, `FusionStyle(DNIrrep{N}) = SimpleFusion()`.

Quantum dimensions are

```math
d_{(0,\pm)} = 1,\qquad d_{\rho_j} = 2,\qquad
d_{(N/2,\pm)} = 1\quad (N \text{ even}).
```

## Topological Data

`DNIrrep{N}` is a bosonic representation category, so all twists are `1`.
The [`Fsymbol`](@ref) is computed from the selected [`fusiontensor`](@ref) basis.
The [`Rsymbol`](@ref) is the N-symbol with a possible sign on odd one-dimensional output channels from two non-trivial two-dimensional inputs.

## Fusion tensor and basis conventions

`fusiontensor(a, b, c)` returns a rank-4 array of size ``d_a \times d_b \times d_c \times N_c^{ab}``.
The implementation uses real Clebsch-Gordan coefficients.
#TODO: source for the convention used here

## Iteration

The `DNIrrep{N}` labels are iterable, with the same ordering as described above.

## References

- [Dihedral group](https://en.wikipedia.org/wiki/Dihedral_group)