```@meta
CollapsedDocStrings = true
```

# Traits and Styles

Traits define compile-time properties of sector types that affect how operations are specialized and optimized.

## FusionStyle

The `FusionStyle` trait is arguably the most important characteristic of a sector, determining how many sectors appear when fusing two sectors, and how many times a unique output can appear.
Various optimizations become available whenever we are not dealing with the `GenericFusion` case.
Firstly, since the shape (size of the arrays) of the topological data is determined by combinations of the [`Nsymbol`](@ref), we can avoid allocating arrays and use scalar quantities for `UniqueFusion` and `SimpleFusion`.
Furthermore, in the `UniqueFusion` case, there is only a single channel for `a ⊗ b ⊗ c ⊗ ...`, paving the way for various optimizations when dealing with fusion trees.

```@docs; canonical = false
FusionStyle
UniqueFusion
SimpleFusion
GenericFusion
```

It is additionally possible to combine fusion styles through the `&` operator, which returns the style with the least assumptions.
For example:

```julia
UniqueFusion() & SimpleFusion() # SimpleFusion()
GenericFusion() & UniqueFusion() # GenericFusion()
```

Finally, some predefined combinations that appear often have dedicated names:

```@docs; canonical = false
MultipleFusion
MultiplicityFreeFusion
```

## BraidingStyle

The `BraidingStyle` describes how sectors behave under exchange (braiding) operations.
In other words, this trait defines the behavior of [`Rsymbol`](@ref) and [`twist`](@ref).
Different braiding styles not only enable different optimizations, but also dictate the allowed operations.
For example, symmetric braiding allows for permutation group statistics, while anyonic systems require full braid group representations.

```@docs; canonical = false
BraidingStyle
NoBraiding
Bosonic
Fermionic
Anyonic
```

It is also possible to combine braiding styles through the `&` operator, which returns the style with the least assumptions.
For example:

```julia
Bosonic() & Fermionic() # Fermionic()
Fermionic() & Anyonic() # Anyonic()
Bosonic() & NoBraiding() # NoBraiding()
```

Finally, some predefined combinations that appear often have dedicated names:

```@docs; canonical = false
HasBraiding
SymmetricBraiding
```

## UnitStyle

The `UnitStyle` describes whether the fusion category has a simple or semisimple unit object.
In other words, this trait determines when we can define a unique value for [`unit`](@ref), or multiple units exist and we have to resort to [`leftunit`](@ref) and [`rightunit`](@ref).
By default, this is derived from `length(allunits(I))`.

```@docs; canonical = false
UnitStyle
SimpleUnit
GenericUnit
```
