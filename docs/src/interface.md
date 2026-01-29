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

### Identity and Duality

Every fusion category has distinguished objects: the unit (identity) and duals (conjugates).
The unit ``\mathbb{1}`` acts as the identity under fusion:
```math
\mathbb{1} ⊗ a ≅ a ≅ a ⊗ \mathbb{1}
```
The dual of a sector ``a`` is the unique sector ``\bar{a}`` such that their fusion contains the unit:
```math
N^{a\bar{a}}_{\mathbb{1}} = 1 \quad \text{and} \quad N^{\bar{a}a}_{\mathbb{1}} = 1
```

```@docs; canonical = false
unit
dual
```

**Example: U₁**
```julia
unit(::Type{U1Irrep}) = U1Irrep(0)         # charge 0 is the unit
dual(c::U1Irrep) = U1Irrep(-charge(c))    # opposite charge
```

**Example: SU₂**
```julia
unit(::Type{SU2Irrep}) = SU2Irrep(0)      # spin 0 is the unit
dual(s::SU2Irrep) = s                      # self-dual
```

**Multifusion categories** can have multiple units and may distinguish between left and right units.
For such cases, additional methods are available:

```@docs; canonical = false
allunits
leftunit
rightunit
```

For regular fusion categories the unit object is unique, such that `unit`, `leftunit` and `rightunit` all coincide.

### Associativity

The associativity of the fusion tensor product is encoded in the F-symbols, which relate different ways of fusing three sectors.
Formally, the F-symbol ``F^{abc}_d`` with intermediate sectors ``e`` and ``f`` is a linear transformation between the two different parenthesizations:
```math
[F^{abc}_d]^f_e : (a ⊗ b → e) ⊗ c → d \quad \longrightarrow \quad a ⊗ (b ⊗ c → f) → d
```
The basis states ``|(a ⊗ b → e) ⊗ c → d\rangle`` are linearly transformed into ``|a ⊗ (b ⊗ c → f) → d\rangle``.
For sectors with `UniqueFusion` or `SimpleFusion`, the F-symbol is a scalar (complex number).
For `GenericFusion`, it is a rank-4 tensor with indices corresponding to the multiplicity labels of each fusion vertex.

The F-symbols must satisfy the **pentagon equation** for every choice of sectors:
```math
\sum_n F^{bcd}_{gn} F^{abn}_{fe} = \sum_m F^{abc}_{em} F^{amc}_{fg} F^{bcd}_{gf}
```
This ensures that all ways of reassociating four tensor factors ``(((a ⊗ b) ⊗ c) ⊗ d)`` to ``(a ⊗ (b ⊗ (c ⊗ d)))`` give the same result, regardless of the sequence of reassociations.

```@docs; canonical = false
Fsymbol
```

**Examples:**
TODO: fix these examples
```julia
# Trivial category: all F-symbols are 1
Fsymbol(::Trivial, ::Trivial, ::Trivial, ::Trivial, ::Trivial, ::Trivial) = 1

# U₁ irreps: all F-symbols are 1 (canonical gauge choice)
Fsymbol(a::U1Irrep, b::U1Irrep, c::U1Irrep, d::U1Irrep, e::U1Irrep, f::U1Irrep) = 1

# SU₂ irreps: computed from 6j-symbols
Fsymbol(a::SU2Irrep, b::SU2Irrep, c::SU2Irrep, d::SU2Irrep, e::SU2Irrep, f::SU2Irrep) = ...
```

### Braiding

Sectors can have a braiding structure that describes the effect of exchanging two tensor factors.
The braiding is encoded in the R-symbol ``R^{ab}_c``, which is a linear transformation between the fusion channels ``a ⊗ b → c`` and ``b ⊗ a → c``.
For sectors with `UniqueFusion` or `SimpleFusion`, the R-symbol is a complex phase.
For `GenericFusion`, it is a square matrix relating the multiplicity spaces of the two fusion orders.

The R-symbols must satisfy the **hexagon equations** together with the F-symbols:
```math
\sum_g R^{ab}_g F^{abx}_{cg} R^{ax}_c = \sum_{f,h} F^{bax}_{cf} R^{af}_h F^{axb}_{ch}
```
and the analogous equation with ``a`` and ``b`` swapped.
These ensure that the braiding is compatible with the associativity encoded by F-symbols.

The type of braiding behavior is declared by the `BraidingStyle` trait, which categorizes sectors into four classes:
- `NoBraiding()` for planar categories where braiding is undefined
- `Bosonic()` for symmetric braiding with trivial twist (all R-symbols square to identity, all twists equal +1)
- `Fermionic()` for symmetric braiding with fermion parity (twists can be ±1)
- `Anyonic()` for general braiding with arbitrary phases or non-symmetric exchange

```@docs; canonical = false
Rsymbol
BraidingStyle
```

**Examples:**
TODO: fix these examples
```julia
# Trivial category: bosonic braiding, all R-symbols are 1
BraidingStyle(::Type{Trivial}) = Bosonic()
Rsymbol(::Trivial, ::Trivial, ::Trivial) = 1

# Fermion parity: fermionic braiding
BraidingStyle(::Type{FermionParity}) = Fermionic()
Rsymbol(a::FermionParity, b::FermionParity, c::FermionParity) = 
    iseven(a) || iseven(b) ? 1 : -1

# Fibonacci anyons: anyonic braiding
BraidingStyle(::Type{FibonacciAnyon}) = Anyonic()
Rsymbol(::FibonacciAnyon, ::FibonacciAnyon, ::FibonacciAnyon) = exp(4π*im/5)

# Planar trivial: no braiding
BraidingStyle(::Type{PlanarTrivial}) = NoBraiding()
```

### Utility Methods

Sectors must support a deterministic ordering and hashing so they can be used as dictionary keys, sorted collections, and canonical fusion outputs.
The ordering should be a strict total order that is consistent with enumerating the sector values, and the hash must be stable with respect to equality.
To achieve this, we must have

- `Base.isless(::Sector, ::Sector)` — Define an order on the sectors.
- `Base.hash(::Sector, h::UInt)` — Associate a hash value with a sector.

**Example: U₁**
```julia
# U₁: order by absolute charge, then prefer positive over negative to comply with `values(U1Irrep)`
function Base.isless(c1::U1Irrep, c2::U1Irrep)
    q1, q2 = charge(c1), charge(c2)
    return abs(q1) < abs(q2) || (abs(q1) == abs(q2) && q1 > q2 > 0)
end

# Hash consistent with equality
Base.hash(c::U1Irrep, h::UInt) = hash(c.charge, h)
```

## Optional Methods

The following methods have default implementations but can be overridden for performance or to provide additional functionality.

### Quantum Dimensions

The quantum dimension of a sector is a fundamental invariant that determines the behavior of fusions and braiding.
The default implementation extracts the dimension from [`Fsymbol`](@ref) via the quantum dimension formula.
For many common sectors, however, the dimension is known directly from representation theory.
In these cases, it can be beneficial to overload `dim` to bypass computing F-symbols, either for performance reasons or to enforce tighter output types.
For example, dimensions of irreducible representations of groups are always integers.

```math
d_a = \left| \frac{1}{(F_{a \bar{a} a}^a)^{1_a}_{_a1} } \right|
```

```@docs; canonical = false
dim
sqrtdim
invsqrtdim
```

### Frobenius-Schur Indicators

The Frobenius-Schur indicator and phase characterize the self-duality properties of sectors.
The indicator distinguishes real, complex, and quaternionic representations.
The phase is the category-theoretic version that appears in line bending operations.

```math
\kappa_a = \text{sign}\left( (F_{a \bar{a} a}^a)^{1_a}_{_a1} \right)
```

```@docs; canonical = false
frobenius_schur_indicator
frobenius_schur_phase
```

### Scalar Type

The scalar type declares what number type is used in the F-symbols and R-symbols for a given sector.
This is automatically inferred from the return types of those symbols, but can be explicitly declared for efficiency.

```@docs; canonical = false
sectorscalartype
```

### Topological Data Symbols

The [`Asymbol`](@ref), [`Bsymbol`](@ref) and [`twist`](@ref) are derived from F- and R-symbols but are often occurring combinations.
The A-symbol and B-symbol relate different ways of bending strands, while the twist is the topological spin (quantum dimension phase) of a sector.

The A-symbol ``A^{ab}_c`` relates splitting and fusion vertices:
```math
A^{ab}_c = \sqrt{\frac{d_c}{d_b}} \overline{\kappa_a} F^{\bar{a} a b}_b
```

The B-symbol ``B^{ab}_c`` relates splitting and fusion vertices:
```math
B^{ab}_c = \sqrt{\frac{d_c}{d_a}} F^{a b \bar{b}}_a R^{ab}_a F^{a \bar{b} b}_a
```

The twist ``\theta_a`` of a sector is the topological spin phase, computed as the trace of the R-matrix for braiding a sector with itself:
```math
\theta_a = \frac{1}{d_a} \sum_{b \in a \otimes a} d_b \text{tr}(R^{aa}_b)
```

```@docs; canonical = false
Asymbol
Bsymbol
twist
```

### Fusion Basis

The fusion tensor provides explicit matrix elements for the tensor product of representations.
It is a rank-4 array whose components are the Clebsch-Gordan coefficients for fusing sector `a` and `b` into `c`.
The fusion tensor is not uniquely determined by topological data alone, instead the topological data can be extracted from it when it is available.
However, there is not always a concrete representation of these tensors in terms of simple `Array` objects, and the exact representation is not important for TensorKit.jl.
For this reason, it is optional: TensorKit can work with the topological data alone.
Note however that they are still required whenever we want to convert symmetric tensors to and from dense arrays.

```@docs; canonical = false
fusiontensor
```

## Traits and Styles

Traits define compile-time properties of sector types that affect how operations are specialized and optimized.

### FusionStyle

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

### BraidingStyle

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

### UnitStyle

The `UnitStyle` describes whether the fusion category has a simple or semisimple unit object.
In other words, this trait determines when we can define a unique value for [`unit`](@ref), or multiple units exist and we have to resort to [`leftunit`](@ref) and [`rightunit`](@ref).
By default, this is derived from `length(allunits(I))`.

```@docs; canonical = false
UnitStyle
SimpleUnit
GenericUnit
```

## Implementation Guidelines

This section provides practical advice for implementing new sector types efficiently and correctly.

### Helper Type: SectorProductIterator

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

### Shape of Topological Data

The dimensionality of [`Fsymbol`](@ref) and [`Rsymbol`](@ref) (and other derived data) can be expressed in terms of the [`Nsymbol`](@ref).
Therefore, depending on the [`FusionStyle`](@ref), we can avoid the allocation of arrays and simply return scalars.

#### [`MultiplicityFreeFusion`](@ref)

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

#### [`GenericFusion`](@ref)

- `Nsymbol(a, b, c)`: Returns a positive `Integer`.
- `Fsymbol(a, b, c, d, e, f)`: Returns a ``N^{ab}_e \times N^{ec}_d \times N^{af}_d \times N^{bc}_f`` array of `sectorscalartype(I)` elements.
- `Rsymbol(a, b, c)`: Returns a ``N^{ab}_c \times N^{ba}_c`` array of `sectorscalartype(I)` elements.

Here invalid fusion channels will automatically lead to empty arrays.
