# Planar Trivial Anyon: `PlanarTrivial`

`PlanarTrivial` is the one-object fusion category without braiding.
It is mostly useful as the smallest planar, non-braided sector in tests and examples.

## Sector type

```@docs; canonical = false
PlanarTrivial
```

There is only one label, `PlanarTrivial()`.
It is its own unit and dual.

## Fusion Rules

Fusion is uniquely trivial:

```math
1 \otimes 1 = 1.
```

`FusionStyle(PlanarTrivial) = UniqueFusion()` and `Nsymbol` is `1` for the only allowed channel.

## Topological Data

The [`Fsymbol`](@ref) is `1`.
The quantum dimension is `1`.

Unlike [`Trivial`](@ref), this sector has `BraidingStyle(PlanarTrivial) = NoBraiding()`.
Consequently, braiding data such as [`Rsymbol`](@ref) and twists are not part of this sector.
