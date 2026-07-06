```@meta
CollapsedDocStrings = true
```

# Required Methods

The following methods **must** be implemented for any new sector type `I <: Sector`.
These methods are grouped by functionality to help understand their purpose.

## Defining the Set of Sectors

A sector type `I <: Sector` represents the set of all labels that can be used to grade a vector space.
This corresponds to all irreducible representations of a group, or all simple objects in a fusion category.
The first requirement is making this set **enumerable** through the iterator interface.

The set of all sector values is obtained via `values(I)`, which returns a `SectorValues{I}()` singleton.
This `SectorValues{I}()` must be iterable, enabling enumeration of all sectors of type `I`.
To do so, one needs to implement the [Iteration interface](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-iteration).
In particular, we require the following methods to be defined:

- `Base.iterate(::SectorValues{I}, state...)` — Iterate over all possible values of sector type `I`.
- `Base.IteratorSize(::Type{SectorValues{I}})` — Specify whether the number of sector values is known, finite, or infinite.

Here the `IteratorSize` is either `HasLength()`, `SizeUnknown()` or `IsInfinite()`.
If the length is known (`HasLength()`), sectors can be indexed by position.
This requires `Base.length`, while `Base.getindex` and `findindex` have generic fallbacks that are typically overridden for performance:

- `Base.length(::SectorValues{I})` — Return the number of sectors. **(required for `HasLength()`)**
- `Base.getindex(::SectorValues{I}, i::Int)` — Access the `i`-th sector value. A fallback linearly iterates the values.
- `findindex(::SectorValues{I}, c::I)` — Find the index of sector `c`. A fallback linearly searches the values.

!!! note
    The choice of `IteratorSize` determines how associative containers are constructed in e.g. a `GradedSpace`.
    In the case of `HasLength()`, an implicit mapping between the values and the position is used to enable storage through `Tuple`s or `Vector`s.
    In the other cases one has to resort to `AbstractDict`-like containers.
    If the set of simple objects is sufficiently large, it might be beneficial to register its length as `SizeUnknown()` to avoid putting too much pressure on the compiler.


## Fusion Structure

The fusion structure describes how sectors combine when taking tensor products.
In code, `a ⊗ b` returns the allowed output labels, and `Nsymbol(a, b, c)` tells you whether (or how many times) `c` appears.

The formal decomposition is:
```math
a ⊗ b = \bigoplus_c N^{ab}_c \, c
```
where ``N^{ab}_c`` is the fusion multiplicity.

The fusion structure is defined through three related methods:

```@docs; canonical = false
⊗
Nsymbol
FusionStyle
```


## Identity and Duality

Every sector type has a unit (identity) label and a notion of dual (conjugate).
`unit` returns the identity label, and `dual` returns the label that fuses with `a` to give the unit.

The unit ``\mathbb{1}`` acts as the identity under fusion:
```math
\mathbb{1} ⊗ a ≅ a ≅ a ⊗ \mathbb{1}
```
The dual of a sector ``a`` is the unique sector ``\bar{a}`` such that:
```math
N^{a\bar{a}}_{\mathbb{1}} = 1 \quad \text{and} \quad N^{\bar{a}a}_{\mathbb{1}} = 1
```

```@docs; canonical = false
unit
dual
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

The associativity of the fusion tensor product tells us how to relate the basis states ``|(a ⊗ b → e) ⊗ c → d\rangle`` to the states ``|a ⊗ (b ⊗ c → f) → d\rangle``.
This is encoded in the F-symbols, which give the coefficients to transform the different ways of fusing three sectors to one.

```@docs; canonical = false
Fsymbol
```

Formally, the F-symbol ``F^{abc}_d`` with intermediate sectors ``e`` and ``f`` is a linear transformation between the two different parenthesizations:

```math
(F_{abc}^d)^e_f : (a ⊗ b → e) ⊗ c → d \quad \longrightarrow \quad a ⊗ (b ⊗ c → f) → d
```

For sectors with `UniqueFusion` or `SimpleFusion`, the F-symbol is a scalar `<:Number`.
For `GenericFusion`, it is a rank-4 tensor with indices corresponding to the multiplicity labels of each fusion vertex.

The F-symbols must satisfy the **pentagon equation** for every choice of sectors:

```math
(F_{fcd}^e)^g_h (F_{abh}^e)^f_i = (F_{abc}^g)^f_j (F_{ajd}^e)^g_i (F_{bcd}^i)^j_h
```

This ensures that all ways of reassociating four tensor factors ``(((a ⊗ b) ⊗ c) ⊗ d)`` to ``(a ⊗ (b ⊗ (c ⊗ d)))`` give the same result, regardless of the sequence of reassociations.


## Braiding

Sectors can have a braiding structure that describes the effect of exchanging two tensor factors.
The braiding is encoded in the R-symbol ``R^{ab}_c``, which is a linear transformation between the fusion channels ``a ⊗ b → c`` and ``b ⊗ a → c``.
For sectors with `UniqueFusion` or `SimpleFusion`, the R-symbol is a complex phase.
For `GenericFusion`, it is a square matrix relating the multiplicity spaces of the two fusion orders.

The R-symbols must satisfy the **hexagon equations** together with the F-symbols:

```math
R^{cd}_e (\overline{F}_{dab}^e)^g_c \overline{R}^{da}_g = (F_{abd}^e)^c_f R^{bd}_f (\overline{F}_{adb}^e)^g_f
```

and the analogous equation with ``a`` and ``b`` swapped.
These ensure that the braiding is compatible with the associativity encoded by F-symbols.

The `BraidingStyle` trait categorizes behavior into four classes:
- `NoBraiding()` for planar categories where braiding is undefined
- `Bosonic()` for symmetric braiding with trivial twist (all R-symbols square to identity, all twists equal +1)
- `Fermionic()` for symmetric braiding with fermion parity (twists can be ±1)
- `Anyonic()` for general braiding with arbitrary phases or non-symmetric exchange

```@docs; canonical = false
Rsymbol
BraidingStyle
```


## Utility Methods

Sectors must support a deterministic ordering and hashing so they can be used as dictionary keys, sorted collections, and canonical fusion outputs.
One should keep the order consistent with `values(I)`, i.e. the enumeration of objects should happen in a sorted fashion.
To achieve this, we must have

- `Base.isless(::Sector, ::Sector)` — Define an order on the sectors.
- `Base.hash(::Sector, h::UInt)` — Associate a hash value with a sector.
