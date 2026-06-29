# Cyclic Group Representations: `ZNIrrep`

`ZNIrrep{N}` and `LargeZNIrrep{N}` represent irreducible representations of the cyclic group ``\mathbb{Z}_N``.
Use `ZNIrrep{N}` when possible; it selects the compact storage type automatically.

## Sector types

```@docs; canonical = false
ZNIrrep
LargeZNIrrep
```

`Z2Irrep`, `Z3Irrep`, and `Z4Irrep` are aliases for `ZNIrrep{2}`, `ZNIrrep{3}`, and `ZNIrrep{4}`.
For small `N`, labels are stored as `UInt8`; larger `N` uses `LargeZNIrrep{N}`.

## Fusion Rules

Labels are charges modulo `N`.
The unit is charge `0`, duality negates the charge modulo `N`, and fusion adds charges:

```math
a \otimes b = (a + b) \bmod N,\qquad a^* = -a \bmod N.
```

The category has `FusionStyle(ZNIrrep{N}) = UniqueFusion()`.

## Topological Data

All quantum dimensions are `1`, as well as all twist.
F-symbols, R-symbols and fusion tensors are trivial in allowed channels.
The braiding style is `Bosonic()`.

## Iteration and basis conventions

`values(ZNIrrep{N})` iterates charges in increasing order from `0` to `N - 1`.

```julia
using TensorKitSectors

values(ZNIrrep{5})[4] # output is Irrep[ℤ{5}](3)
values(ZNIrrep{34})[34] # output is Irrep[ℤ{34}](33)
```

## References

- [Cyclic group](https://en.wikipedia.org/wiki/Cyclic_group)
- [Representation theory of finite groups](https://en.wikipedia.org/wiki/Representation_theory_of_finite_groups)
