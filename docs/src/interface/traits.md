```@meta
CollapsedDocStrings = true

DocTestSetup = quote
    using TensorKitSectors
end
```

# Traits and Styles

Traits define compile-time properties that can be assumed about a sector type.
They control behavior and enable optimizations.

## FusionStyle

The `FusionStyle` trait indicates how many outputs to expect when fusing two sectors.

```@docs; canonical = false
FusionStyle
UniqueFusion
SimpleFusion
GenericFusion
```

This enables various optimizations for different cases.
Firstly, since the shape (size of the arrays) of the topological data is determined by combinations of the [`Nsymbol`](@ref), for `UniqueFusion` and `SimpleFusion` we can use scalar quantities instead of arrays.
Secondly, in the `UniqueFusion` case, there is only a single channel for ``a \otimes b \otimes c \otimes \ldots``, avoiding the need to iterate through all options.

It is additionally possible to combine fusion styles through the `&` operator, which returns the style with the least assumptions.
For example:

```jldoctest
julia> UniqueFusion() & SimpleFusion()
SimpleFusion()

julia> GenericFusion() & UniqueFusion()
GenericFusion()
```

Finally, some predefined combinations that appear often have dedicated names:

```@docs; canonical = false
MultipleFusion
MultiplicityFreeFusion
```

## BraidingStyle

The `BraidingStyle` describes whether and how exchange of sectors is defined.
It determines how TensorKit interprets [`Rsymbol`](@ref) and [`twist`](@ref).

```@docs; canonical = false
BraidingStyle
NoBraiding
Bosonic
Fermionic
Anyonic
```

Additionally, this dictates whether or not permutations are sufficient to specify generic exchanges, or if a full braid group representation is needed.

It is also possible to combine braiding styles through the `&` operator, which returns the style with the least assumptions.
For example:

```jldoctest
julia> Bosonic() & Fermionic()
Fermionic()

julia> Fermionic() & Anyonic()
Anyonic()

julia> Bosonic() & NoBraiding()
NoBraiding()
```

Finally, some predefined combinations that appear often have dedicated names:

```@docs; canonical = false
HasBraiding
SymmetricBraiding
```

## UnitStyle

The `UnitStyle` tells whether there is a single identity label or multiple units.
By default, it is derived from `length(allunits(I))`.

```@docs; canonical = false
UnitStyle
SimpleUnit
GenericUnit
```

Whenever the style is `SimpleUnit`, a unique value of [`unit`](@ref) can be defined and there is no distinction between [`leftunit`](@ref) and [`rightunit`](@ref).
For `GenericUnit`, this is no longer the case and special care has to be taken to use the *correct* unit for various fusion diagrams.
