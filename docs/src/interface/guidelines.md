```@meta
CollapsedDocStrings = true
```

# Implementation Guidelines

This section provides practical advice for implementing new sector types efficiently and correctly.

## Helper Type: SectorProductIterator

Instead of materializing all fusion outputs in a tuple or array, one can use `SectorProductIterator` for type-stable, lazy iteration.

```@docs; canonical = false
TensorKitSectors.SectorProductIterator
```

This is enabled by default, and new sectors should provide implementations for the following functions:

```julia
Base.iterate(ab::SectorProductIterator{I}, [state]) = ...
# optional optimizations:
Base.IteratorSize(::Type{SectorProductIterator{I}}) = HasLength()
Base.length(ab::SectorProductIterator{I}) = ...
```

This has the benefit of helping with type stability, and has a pretty-printing overload:

```@example
using TensorKitSectors # hide
SU2Irrep(1) ⊗ SU2Irrep(1)
```

## Shape of Topological Data

The dimensionality of [`Fsymbol`](@ref) and [`Rsymbol`](@ref) (and other derived data) can be expressed in terms of the [`Nsymbol`](@ref).
Therefore, depending on the [`FusionStyle`](@ref), we can avoid the allocation of arrays and simply return scalars.

### [`MultiplicityFreeFusion`](@ref)

- `Nsymbol(a, b, c)`: Returns a `Bool`.
- `Fsymbol(a, b, c, d, e, f)`: Returns a scalar of type `sectorscalartype(I)`.
- `Rsymbol(a, b, c)`: Returns a scalar of type `sectorscalartype(I)`.

Additionally, if the [`Fsymbol`](@ref) and [`Rsymbol`](@ref) do not correspond to valid fusion channels, the result is ``0``.

In other words, we have:

```math
\left(N^{ab}_e N^{ec}_d = 0 \lor N^{af}_d N^{bc}_f = 0\right) \implies (F_{abc}^d)_e^f = 0
```

```math
N^{ab}_c = 0 \implies R_{ab}^c = 0
```

### [`GenericFusion`](@ref)

- `Nsymbol(a, b, c)`: Returns a positive `Integer`.
- `Fsymbol(a, b, c, d, e, f)`: Returns a ``N^{ab}_e \times N^{ec}_d \times N^{af}_d \times N^{bc}_f`` array of `sectorscalartype(I)` elements.
- `Rsymbol(a, b, c)`: Returns a ``N^{ab}_c \times N^{ba}_c`` array of `sectorscalartype(I)` elements.

Here invalid fusion channels will automatically lead to empty arrays.
