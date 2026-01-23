```@meta
CollapsedDocStrings = true
```

# Sector Interface

## Introduction

A `Sector` represents the isomorphism classes of simple objects in unitary and pivotal fusion categories.
In the context of tensor networks and representation theory, sectors label the different *symmetry charges* or *quantum numbers* that organize graded vector spaces.

Mathematically, sectors encode the topological data of fusion categories:
- **Simple objects**: Irreducible representations or anyonic excitations
- **Fusion rules**: How sectors combine (``a ⊗ b → \bigoplus_c N^{ab}_c c``)
- **Associativity**: F-symbols (6j symbols) encoding recoupling transformations
- **Braiding**: R-symbols encoding exchange statistics (for braided categories)

Sectors can be used as labels for graded vector spaces, which represent vector spaces decomposed into different sectors:
```math
V = \bigoplus_{a \in \text{Sectors}} \mathbb{C}^{n_a} \otimes V_a
```
where ``n_a`` is the multiplicity of sector ``a`` and ``V_a`` is the associated vector space.

This page explains the required and optional methods needed to create a new sector type.
Once this interface is fulfilled, TensorKit.jl should be able to create symmetric tensors that correspond to the implemented sector.

## Required Methods

The following methods **must** be implemented for any new sector type `I <: Sector`.
These methods are grouped by functionality to help understand their purpose.

### Defining the Set of Sectors

A sector type `I <: Sector` represents the set of all simple objects in a fusion category — equivalently, all irreducible representations of a group.
The first requirement is defining this set by making it **enumerable** through the iterator interface.

The set of all sector values is obtained via `values(I)`, which returns a `SectorValues{I}()` singleton type by default.
This `SectorValues{I}()` must be iterable, enabling enumeration of all sectors of type `I`.
To do so, one needs to implement the [Iteration interface](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-iteration).
In particular, we require the following methods to be defined:

- `Base.iterate(::SectorValues{I}, state...)` — Iterate over all possible values of sector type `I`.
- `Base.IteratorSize(::Type{SectorValues{I}})` — Specify whether the number of sector values is known, finite, or infinite.

Here the `IteratorSize` is either `HasLength()`, `SizeUnknown()` or `IsInfinite()`.
If the length is known (`HasLength()`), we additionally require:

- `Base.length(::SectorValues{I})`: Return the number of sectors

TODO: this is technically optional!
- `Base.getindex(::SectorValues{I}, i::Int)`: Access the i-th sector value
- `findindex(::SectorValues{I}, c::I)`: Find the index of sector `c`

!!! note
    The choice of `IteratorSize` determines how associative containers are constructed in e.g. a `GradedSpace`.
    In the case of `HasLength()`, an implicit mapping between the values and the position is used to enable storage through `Tuple`s or `Vector`s.
    In the other cases one has to resort to `AbstractDict`-like containers.
    If the set of simple objects is sufficiently large, it might be beneficial to register its length as `SizeUnknown()` to avoid putting too much pressure on the compiler.

**Example: `Z3Irrep`**
```julia
# ℤ₃ has exactly 3 irreps labeled by 0, 1, 2
Base.iterate(::SectorValues{Z3Irrep}, i = 0) = i > 2 ? nothing : (Z2Irrep(i), i+1)
Base.IteratorSize(::Type{SectorValues{Z3Irrep}}) = HasLength()
Base.length(::SectorValues{Z3Irrep}) = 3
Base.getindex(::SectorValues{Z3Irrep}, i::Int) = Z3Irrep(i-1)
findindex(::SectorValues{Z3Irrep}, c::Z3Irrep) = c.n + 1
```

**Example: `U1Irrep`**
```julia
# U₁ an irrep for each integer: 0, 1, -1, 2, -2, ...
Base.iterate(::SectorValues{U1Irrep}, i = 0) = (U1Irrep(i), i <= 0 ? (-i + 1) : -1)
Base.IteratorSize(::Type{SectorValues{U1Irrep}}) = IsInfinite()
```

### Fusion Structure

The fusion structure describes how sectors combine when taking tensor products.
Mathematically, this is the tensor product decomposition of irreps or simple objects:
```math
a ⊗ b = \bigoplus_c N^{ab}_c \, c
```
where ``N^{ab}_c`` is the fusion multiplicity indicating how many times sector ``c`` appears in the fusion of ``a`` and ``b``.

The fusion structure is defined through three related methods:

```@docs; canonical = false
⊗
Nsymbol
FusionStyle
```

**Example: U₁ (Unique Fusion)**
```julia
⊗(c1::U1Irrep, c2::U1Irrep) = (U1Irrep(charge(c1) + charge(c2)),)
Nsymbol(c1::U1Irrep, c2::U1Irrep, c3::U1Irrep) = charge(c1) + charge(c2) == charge(c3)
FusionStyle(::Type{U1Irrep}) = UniqueFusion()
```

**Example: SU₂ (Simple Fusion)**
TODO: finish this example
```julia
function ⊗(s1::SU2Irrep, s2::SU2Irrep)
    return SectorProductIterator(s1, s2)  # Yields SU2Irrep(j) for j ∈ [|j1-j2|, j1+j2]
end

function Nsymbol(s1::SU2Irrep, s2::SU2Irrep, s3::SU2Irrep)
    j1, j2, j3 = spin(s1), spin(s2), spin(s3)
    return abs(j1 - j2) <= j3 <= j1 + j2 && isinteger(j1 + j2 + j3)
end

FusionStyle(::Type{SU2Irrep}) = SimpleFusion()
```

### Topological Data

### Utility Methods

## Optional Methods

## Traits and Styles

### FusionStyle

### BraidingStyle

### UnitStyle

## Implementation Guidelines

### Helper Type: SectorProductIterator

### Shape of Topological Data (determined by FusionStyle)

### Key Considerations

### Registering with IrrepTable

### Common Patterns
