```@meta
CollapsedDocStrings = true
```

# Required Methods

The following methods **must** be implemented for any new sector type `I <: Sector`.
These methods are grouped by functionality to help understand their purpose.

## Defining the Set of Sectors

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
Base.iterate(::SectorValues{U1Irrep}, i = 0) = (U1Irrep(i), i <= 0 ? (-i + 1) : -i - 1)
Base.IteratorSize(::Type{SectorValues{U1Irrep}}) = IsInfinite()
```

## Fusion Structure

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

## Identity and Duality

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

## Associativity

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

## Braiding

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

## Utility Methods

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
