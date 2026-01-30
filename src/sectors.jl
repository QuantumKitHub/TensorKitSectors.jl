"""
    abstract type Sector

Abstract type for representing the (isomorphism classes of) simple objects in (unitary
and pivotal) (pre-)fusion categories, e.g. the irreducible representations of a finite or
compact group. Subtypes `I <: Sector` as the set of labels of a `GradedSpace`.

Every new `I <: Sector` should implement the following methods:
*   `unit(::Type{I})`: unit element of `I`. If there are multiple, implement `allunits(::Type{I})`
    instead.
*   `dual(a::I)`: ``a̅``, conjugate or dual label of ``a``
*   `⊗(a::I, b::I)`: iterable with unique fusion outputs of ``a ⊗ b``
    (i.e. don't repeat in case of multiplicities)
*   `Nsymbol(a::I, b::I, c::I)`: number of times `c` appears in `a ⊗ b`, i.e. the
    multiplicity
*   `FusionStyle(::Type{I})`: `UniqueFusion()`, `SimpleFusion()` or
    `GenericFusion()`
*   `BraidingStyle(::Type{I})`: `Bosonic()`, `Fermionic()`, `Anyonic()`, ...
*   `Fsymbol(a::I, b::I, c::I, d::I, e::I, f::I)`: F-symbol: scalar (in case of
    `UniqueFusion`/`SimpleFusion`) or rank-4 array (in case of `GenericFusion`)
*   `Rsymbol(a::I, b::I, c::I)`: R-symbol: scalar (in case of
    `UniqueFusion`/`SimpleFusion`) or matrix (in case of `GenericFusion`)
*   `isless(a::I, b::I)`: defines a canonical ordering of sectors
*   `hash(a::I)`: hash function for sectors
and optionally
*   `dim(a::I)`: quantum dimension of sector `a`
*   `frobenius_schur_indicator(a::I)`: Frobenius-Schur indicator of `a` (1, 0, -1)
*   `frobenius_schur_phase(a::I)`: Frobenius-Schur phase of `a` (±1)
*   `sectorscalartype(::Type{I})`: scalar type of F- and R-symbols
*   `Bsymbol(a::I, b::I, c::I)`: B-symbol: scalar (in case of
    `UniqueFusion`/`SimpleFusion`) or matrix (in case of `GenericFusion`)
*   `twist(a::I)` -> twist of sector `a`

Furthermore, `iterate` and `Base.IteratorSize` should be made to work for the singleton type
[`SectorValues{I}`](@ref).

To help with the implementation of `⊗(a::I, b::I)` as an iterator, the provided `struct` type
[`SectorProductIterator{I}`](@ref) can be used, which stores `a` and `b` and requires the
implementation of `Base.iterate(::SectorProductIterator{I}, state...)`.
"""
abstract type Sector end

"""
    type_repr(T::Type)

Return a string representation of the type `T`, which is used to modify the default
way in which `Sector` subtypes are displayed in other objects that depend on them.
"""
type_repr(T::Type) = repr(T)

# iterator over the values (i.e., elements of representative set of simple objects)
# in the sector
"""
    struct SectorValues{I <: Sector}

Singleton type to represent an iterator over the possible values of type `I`, whose
instance is obtained as `values(I)`. For a new `I::Sector`, the following should be defined
*   `Base.iterate(::SectorValues{I}[, state])`: iterate over the values
*   `Base.IteratorSize(::Type{SectorValues{I}})`: `HasLength()`, `SizeUnknown()`
    or `IsInfinite()` depending on whether the number of values of type `I` is finite
    (and sufficiently small) or infinite; for a large number of values, `SizeUnknown()` is
    recommended because this will trigger the use of `GenericGradedSpace`.
If `IteratorSize(I) == HasLength()`, also the following must be implemented:
*   `Base.length(::SectorValues{I})`: the number of different values
*   `Base.getindex(::SectorValues{I}, i::Int)`: a mapping between an index `i` and an
    instance of `I`. A fallback implementation exists that returns the `i`th value
    of the `SectorValues` iterator.
*   `findindex(::SectorValues{I}, c::I)`: reverse mapping between a value `c::I` and an
    index `i::Integer ∈ 1:length(values(I))`. A fallback implementation exists that
    linearly searches through the `SectorValues` iterator.
"""
struct SectorValues{I <: Sector} end
Base.IteratorEltype(::Type{<:SectorValues}) = HasEltype()
Base.eltype(::Type{SectorValues{I}}) where {I <: Sector} = I
Base.values(::Type{I}) where {I <: Sector} = SectorValues{I}()

Base.@propagate_inbounds function Base.getindex(
        v::SectorValues{I}, i::Int
    ) where {I <: Sector}
    @boundscheck begin
        if Base.IteratorSize(v) === HasLength()
            1 ≤ i ≤ length(v) || throw(BoundsError(v, i))
        else
            1 ≤ i || throw(BoundsError(v, i))
        end
    end
    for (j, c) in enumerate(v)
        j == i && return c
    end
    throw(BoundsError(v, i))
end

"""
    findindex(v::SectorValues{I}, c::I)

Reverse mapping between a value `c::I` and an index `i::Integer ∈ 1:length(values(I))`.
"""
function findindex(v::SectorValues{I}, c::I) where {I <: Sector}
    for (i, cc) in enumerate(v)
        cc == c && return i
    end
    throw(ArgumentError(lazy"Cannot locate sector $c"))
end


"""
    allunits(I::Type{<:Sector}) -> Tuple{I}

Return a tuple with all units of the sector type `I`.
For fusion categories, this will contain only one element.
"""
allunits(::Type{I}) where {I <: Sector} = (unit(I),)

"""
    unit(::Sector) -> Sector
    unit(::Type{<:Sector}) -> Sector

Return the unit element of this type of sector, provided it is unique.
"""
unit(a::Sector) = unit(typeof(a))
Base.one(a::Sector) = unit(a)
Base.one(::Type{I}) where {I <: Sector} = unit(I)

"""
    isunit(a::Sector) -> Bool

Return whether sector `a` is a unit element.
"""
function isunit(a::Sector)
    return if UnitStyle(a) === SimpleUnit()
        a == unit(a)
    else
        a in allunits(typeof(a))
    end
end
Base.isone(a::Sector) = isunit(a)

"""
    leftunit(a::Sector) -> Sector

Return the left unit element corresponding to `a`;
this is necessary for multifusion categories, where the unit may not be unique.
See also [`rightunit`](@ref) and [`unit`](@ref).
"""
leftunit(a::Sector) = unit(a)

"""
    rightunit(a::Sector) -> Sector

Return the right unit element corresponding to `a`;
this is necessary for multifusion categories, where the unit may not be unique.
See also [`leftunit`](@ref) and [`unit`](@ref).
"""
rightunit(a::Sector) = unit(a)

@doc """
    dual(a::Sector) -> Sector

Return the dual label of `a`, i.e. the unique label `ā = dual(a)` such that 
`Nsymbol(a, ā, leftunit(a)) == 1` and `Nsymbol(ā, a, rightunit(a)) == 1`.
""" dual(::Sector)
Base.conj(a::Sector) = dual(a)

"""
    sectorscalartype(I::Type{<:Sector}) -> Type{<:Number}

Return the scalar type of the topological data of the sector `I`.
In particular, this is a combination of the scalar type of both the [`Fsymbol`](@ref) and [`Rsymbol`](@ref),
and determines the scalar type of the [`fusiontensor`](@ref) whenever it is defined.

See also [`fusionscalartype`](@ref) and [`braidingscalartype`](@ref).
"""
function sectorscalartype(::Type{I}) where {I <: Sector}
    if BraidingStyle(I) === NoBraiding()
        return fusionscalartype(I)
    else
        return Base.promote_op(*, fusionscalartype(I), braidingscalartype(I))
    end
end

"""
    fusionscalartype(I::Type{<:Sector}) -> Type{<:Number}

Return the scalar type of the topological data associated to fusion of the sector `I`.
In particular, this is the scalar type of [`Fsymbol`](@ref).

See also [`braidingscalartype`](@ref) and [`sectorscalartype`](@ref).
"""
function fusionscalartype(::Type{I}) where {I <: Sector}
    Ftype = Core.Compiler.return_type(Fsymbol, NTuple{6, I})
    return FusionStyle(I) === UniqueFusion() ? Ftype : eltype(Ftype)
end

"""
    braidingscalartype(I::Type{<:Sector}) -> Type{<:Number}

Return the scalar type of the topological data associated to braiding of the sector `I`.
In particular, this is the scalar type of [`Rsymbol`](@ref).

See also [`fusionscalartype`](@ref) and [`sectorscalartype`](@ref).
"""
function braidingscalartype(::Type{I}) where {I <: Sector}
    BraidingStyle(I) === NoBraiding() && throw(ArgumentError("No braiding for sector $I"))
    Rtype = Core.Compiler.return_type(Rsymbol, NTuple{3, I})
    return FusionStyle(I) === UniqueFusion() ? Rtype : eltype(Rtype)
end

"""
    isreal(::Type{<:Sector}) -> Bool

Return whether the topological data (Fsymbol, Rsymbol) of the sector is real or not
(in which case it is complex).
"""
Base.isreal(I::Type{<:Sector}) = sectorscalartype(I) <: Real

# FusionStyle: the most important aspect of Sector
#---------------------------------------------
"""
    ⊗(a::I, b::I...) where {I <: Sector}
    otimes(a::I, b::I...) where {I <: Sector}

Return an iterable of elements of `c::I` that appear in the fusion product `a ⊗ b`.

Note that every element `c` should appear at most once, fusion degeneracies (if
`FusionStyle(I) == GenericFusion()`) should be accessed via `Nsymbol(a, b, c)`.
"""
function ⊗ end
const otimes = ⊗

⊗(I::Sector) = (I,)

# NOTE: the following inline is extremely important for performance, especially
# in the case of UniqueFusion, because ⊗(...) is computed very often
@inline function ⊗(a::I, b::I, c::I, rest::Vararg{I}) where {I <: Sector}
    if FusionStyle(I) isa UniqueFusion
        return a ⊗ first(⊗(b, c, rest...))
    else
        s = Set{I}()
        for d in ⊗(b, c, rest...)
            for e in a ⊗ d
                push!(s, e)
            end
        end
        return s
    end
end

"""
    struct SectorProductIterator{I <: Sector}
    SectorProductIterator(a::I, b::I) where {I <: Sector}

Custom iterator to represent the (unique) fusion outputs of ``a ⊗ b``.

Custom sectors that aim to use this have to provide the following functionality:

* `Base.iterate(::SectorProductIterator{I}, state...) where {I <: Sector}`: iterate over
    the fusion outputs of `a ⊗ b`

If desired and it is possible to easily compute the number of unique fusion outputs, it is also
possible to define `Base.IteratorSize(::Type{SectorProductIterator{I}}) = Base.HasLength()`, in which
case `Base.length(::SectorProductIterator{I})` has to be implemented.

See also [`⊗`](@ref).
"""
struct SectorProductIterator{I <: Sector}
    a::I
    b::I
end

⊗(a::I, b::I) where {I <: Sector} = SectorProductIterator(a, b)

Base.IteratorSize(::Type{SectorProductIterator{I}}) where {I} = Base.SizeUnknown()
Base.IteratorEltype(::Type{SectorProductIterator{I}}) where {I} = Base.HasEltype()
Base.eltype(::Type{SectorProductIterator{I}}) where {I} = I

function Base.show(io::IO, it::SectorProductIterator)
    show(io, it.a)
    print(io, " ⊗ ")
    show(io, it.b)
    return nothing
end

function Base.show(io::IO, mime::MIME"text/plain", ab::SectorProductIterator)
    show(io, ab.a)
    print(io, " ⊗ ")
    show(io, ab.b)
    get(io, :compact, false) && return nothing
    print(io, ":")
    ioc = IOContext(io, :typeinfo => eltype(ab))
    for c in ab
        print(io, "\n ")
        show(ioc, mime, c)
    end
    return nothing
end

"""
    Nsymbol(a::I, b::I, c::I) where {I <: Sector} -> Integer

Return an `Integer` representing the number of times `c` appears in the fusion product
`a ⊗ b`. Could be a `Bool` if `FusionStyle(I) == UniqueFusion()` or `SimpleFusion()`.
"""
function Nsymbol end

# trait to describe the fusion of superselection sectors
"""
    abstract type FusionStyle
    FusionStyle(::Sector)
    FusionStyle(I::Type{<:Sector})

Trait to describe the fusion behavior of sectors of type `I`, which can be either
*   `UniqueFusion()`: single fusion output when fusing two sectors;
*   `SimpleFusion()`: multiple outputs, but every output occurs at most one,
    also known as multiplicity-free (e.g. irreps of ``SU(2)``);
*   `GenericFusion()`: multiple outputs that can occur more than once (e.g. irreps
    of ``SU(3)``).

There is an abstract supertype `MultipleFusion` of which both `SimpleFusion` and
`GenericFusion` are subtypes. Furthermore, there is a type alias `MultiplicityFreeFusion`
for those fusion types which do not require muliplicity labels, i.e.
`MultiplicityFreeFusion = Union{UniqueFusion,SimpleFusion}`.
"""
abstract type FusionStyle end
FusionStyle(a::Sector) = FusionStyle(typeof(a))

struct UniqueFusion <: FusionStyle end # unique fusion output when fusing two sectors
abstract type MultipleFusion <: FusionStyle end
struct SimpleFusion <: MultipleFusion end # multiple fusion but multiplicity free
struct GenericFusion <: MultipleFusion end # multiple fusion with multiplicities
const MultiplicityFreeFusion = Union{UniqueFusion, SimpleFusion}

@doc (@doc FusionStyle) UniqueFusion
@doc (@doc FusionStyle) SimpleFusion
@doc (@doc FusionStyle) GenericFusion

# combine fusion properties of tensor products of sectors
Base.:&(f::F, ::F) where {F <: FusionStyle} = f
Base.:&(f₁::FusionStyle, f₂::FusionStyle) = f₂ & f₁

Base.:&(::SimpleFusion, ::UniqueFusion) = SimpleFusion()
Base.:&(::GenericFusion, ::UniqueFusion) = GenericFusion()
Base.:&(::GenericFusion, ::SimpleFusion) = GenericFusion()

# similar, but for multifusion categories
"""
    abstract type UnitStyle
    UnitStyle(::Sector)
    UnitStyle(I::Type{<:Sector})

Trait to describe the semisimplicity of the unit sector of type `I`.
This can be either
*   `SimpleUnit()`: the unit is simple (e.g. fusion categories);
*   `GenericUnit()`: the unit is semisimple.
"""
abstract type UnitStyle end
UnitStyle(a::Sector) = UnitStyle(typeof(a))

struct SimpleUnit <: UnitStyle end
struct GenericUnit <: UnitStyle end

@doc (@doc UnitStyle) SimpleUnit
@doc (@doc UnitStyle) GenericUnit

UnitStyle(::Type{I}) where {I <: Sector} = length(allunits(I)) == 1 ? SimpleUnit() : GenericUnit()

@noinline function throw_genericunit_error(I)
    throw(DomainError(I, "Sector has multiple units, use `allunits` instead of `unit`"))
end

# combine fusion properties of tensor products of multifusion sectors
Base.:&(f::F, ::F) where {F <: UnitStyle} = f
Base.:&(f₁::UnitStyle, f₂::UnitStyle) = f₂ & f₁

Base.:&(::GenericUnit, ::SimpleUnit) = GenericUnit()

@doc """
    fusiontensor(a::I, b::I, c::I) where {I <: Sector} -> AbstractArray{T, 4}

Return the fusion tensor for the fusion `a ⊗ b -> c`. The dimensions of the returned array
are `(dim(a), dim(b), dim(c), Nsymbol(a, b, c))`. The components of the fusion tensor are
simply the Clebsch-Gordan coefficients, describing the unitary basis change from the tensor
product of irreps `a` and `b` to the coupled irrep `c`.
""" fusiontensor(::I, ::I, ::I) where {I <: Sector}

"""
    Fsymbol(a::I, b::I, c::I, d::I, e::I, f::I) where {I <: Sector}

Return the F-symbol ``F^{abc}_d`` that associates the two different fusion orders of sectors
`a`, `b` and `c` into an ouput sector `d`, using either an intermediate sector ``a ⊗ b → e``
or ``b ⊗ c → f``:
```
a-<-μ-<-e-<-ν-<-d                                     a-<-λ-<-d
    ∨       ∨       -> Fsymbol(a,b,c,d,e,f)[μ,ν,κ,λ]      ∨
    b       c                                             f
                                                          v
                                                      b-<-κ
                                                          ∨
                                                          c
```
If `FusionStyle(I)` is `UniqueFusion` or `SimpleFusion`, the F-symbol is a number. Otherwise
it is a rank 4 array of size
`(Nsymbol(a, b, e), Nsymbol(e, c, d), Nsymbol(b, c, f), Nsymbol(a, f, d))`.
"""
function Fsymbol end

# properties that can be determined in terms of the F symbol
# TODO: find mechanism for returning these numbers with custom type T <: AbstractFloat
"""
    dim(a::Sector)

Return the (quantum) dimension of the sector `a`.
"""
function dim(a::Sector)
    return if FusionStyle(a) isa UniqueFusion
        1
    elseif FusionStyle(a) isa SimpleFusion
        abs(1 / Fsymbol(a, dual(a), a, a, leftunit(a), rightunit(a)))
    else
        abs(1 / Fsymbol(a, dual(a), a, a, leftunit(a), rightunit(a))[1])
    end
end
sqrtdim(a::Sector) = (FusionStyle(a) isa UniqueFusion) ? 1 : sqrt(dim(a))
invsqrtdim(a::Sector) = (FusionStyle(a) isa UniqueFusion) ? 1 : inv(sqrt(dim(a)))

"""
    frobenius_schur_phase(a::Sector)

Return the Frobenius-Schur phase ``κₐ`` of a sector ``a``, which is a complex phase that
appears in the context of bending lines and is obtained from ``F^{a a̅ a}_a``.
When `a == dual(a)`, it is restricted to ``κₐ ∈ \\{1, -1\\}`` and coincides with
the group-theoretic version [`frobenius_schur_indicator`](@ref).
When `a != dual(a)`, the value of ``κₐ`` can be gauged to be `1`, though is not required to be.
"""
function frobenius_schur_phase(a::Sector)
    return if FusionStyle(a) isa UniqueFusion || FusionStyle(a) isa SimpleFusion
        sign(Fsymbol(a, dual(a), a, a, leftunit(a), rightunit(a)))
    else
        sign(Fsymbol(a, dual(a), a, a, leftunit(a), rightunit(a))[1])
    end
end

"""
    frobenius_schur_indicator(a::Sector)

Return the Frobenius-Schur indicator of a sector ``νₐ ∈ \\{1, 0, -1\\}``, which distinguishes
between real, complex and quaternionic representations.

See also [`frobenius_schur_phase`](@ref) for the category-theoretic version that appears in the
context of line bending.
"""
function frobenius_schur_indicator(a::Sector)
    ν = frobenius_schur_phase(a)
    return a == conj(a) ? ν : zero(ν)
end

# Not necessary
"""
    Asymbol(a::I, b::I, c::I) where {I <: Sector}

Return the value of ``A^{ab}_c`` which appears in transforming a splitting vertex
into a fusion vertex using the transformation
```
a -<-μ-<- c                                                    b -<-ν-<- dual(a)
     ∨       -> √(dim(c) / dim(b)) * Asymbol(a, b, c)[μ, ν]         ∧
     b                                                              c
```
If `FusionStyle(I)` is `UniqueFusion()` or `SimpleFusion()`, the A-symbol is a
number. Otherwise it is a square matrix with row and column size
`Nsymbol(a, b, c) == Nsymbol(dual(a), c, b)`.
"""
function Asymbol(a::I, b::I, c::I) where {I <: Sector}
    return if FusionStyle(I) isa UniqueFusion || FusionStyle(I) isa SimpleFusion
        (sqrtdim(a) * sqrtdim(b) * invsqrtdim(c)) *
            conj(frobenius_schur_phase(a) * Fsymbol(dual(a), a, b, b, leftunit(a), c))
    else
        reshape(
            (sqrtdim(a) * sqrtdim(b) * invsqrtdim(c)) *
                conj(frobenius_schur_phase(a) * Fsymbol(dual(a), a, b, b, leftunit(a), c)),
            (Nsymbol(a, b, c), Nsymbol(dual(a), c, b))
        )
    end
end

"""
    Bsymbol(a::I, b::I, c::I) where {I <: Sector}

Return the value of ``B^{ab}_c`` which appears in transforming a splitting vertex
into a fusion vertex using the transformation
```
a -<-μ-<- c                                                    a -<-ν-<- c
     ∨       -> √(dim(c) / dim(a)) * Bsymbol(a, b, c)[μ, ν]         ∧
     b                                                            dual(b)
```
If `FusionStyle(I)` is `UniqueFusion()` or `SimpleFusion()`, the B-symbol is a
number. Otherwise it is a square matrix with row and column size
`Nsymbol(a, b, c) == Nsymbol(c, dual(b), a)`.
"""
function Bsymbol(a::I, b::I, c::I) where {I <: Sector}
    return if FusionStyle(I) isa UniqueFusion || FusionStyle(I) isa SimpleFusion
        (sqrtdim(a) * sqrtdim(b) * invsqrtdim(c)) * Fsymbol(a, b, dual(b), a, c, rightunit(a))
    else
        reshape(
            (sqrtdim(a) * sqrtdim(b) * invsqrtdim(c)) * Fsymbol(a, b, dual(b), a, c, rightunit(a)),
            (Nsymbol(a, b, c), Nsymbol(c, dual(b), a))
        )
    end
end

# Braiding:
#-------------------------------------------------
# trait to describe type to denote how the elementary spaces in a tensor product space
# interact under permutations or actions of the braid group
"""
    abstract type BradingStyle
    BraidingStyle(::Sector) -> ::BraidingStyle
    BraidingStyle(I::Type{<:Sector}) -> ::BraidingStyle

Return the type of braiding and twist behavior of sectors of type `I`, which can be either
*   `NoBraiding()`: no braiding structure
*   `Bosonic()`: symmetric braiding with trivial twist (i.e. identity)
*   `Fermionic()`: symmetric braiding with non-trivial twist (squares to identity)
*   `Anyonic()`: general ``R^{ab}_c`` phase or matrix (depending on `SimpleFusion` or
    `GenericFusion` fusion) and arbitrary twists

Note that `Bosonic` and `Fermionic` are subtypes of `SymmetricBraiding`, which means that
braids are in fact equivalent to crossings (i.e. braiding twice is an identity:
`isone(Rsymbol(b,a,c)*Rsymbol(a,b,c)) == true`) and permutations are uniquely defined.
"""
abstract type BraidingStyle end
BraidingStyle(a::Sector) = BraidingStyle(typeof(a))

abstract type HasBraiding <: BraidingStyle end
struct NoBraiding <: BraidingStyle end
abstract type SymmetricBraiding <: HasBraiding end # symmetric braiding => actions of permutation group are well defined
struct Bosonic <: SymmetricBraiding end # all twists are one
struct Fermionic <: SymmetricBraiding end # twists one and minus one
struct Anyonic <: HasBraiding end

@doc (@doc BraidingStyle) NoBraiding
@doc (@doc BraidingStyle) Bosonic
@doc (@doc BraidingStyle) Fermionic
@doc (@doc BraidingStyle) Anyonic

Base.:&(b::B, ::B) where {B <: BraidingStyle} = b
Base.:&(B1::BraidingStyle, B2::BraidingStyle) = B2 & B1
Base.:&(::Bosonic, ::Fermionic) = Fermionic()
Base.:&(::Bosonic, ::Anyonic) = Anyonic()
Base.:&(::Fermionic, ::Anyonic) = Anyonic()
Base.:&(::Bosonic, ::NoBraiding) = NoBraiding()
Base.:&(::Fermionic, ::NoBraiding) = NoBraiding()
Base.:&(::Anyonic, ::NoBraiding) = NoBraiding()

"""
    Rsymbol(a::I, b::I, c::I) where {I <: Sector}

Returns the R-symbol ``R^{ab}_c`` that maps between ``c → a ⊗ b`` and ``c → b ⊗ a`` as in
```
a -<-μ-<- c                                 b -<-ν-<- c
     ∨        -> Rsymbol(a, b, c)[μ, ν]          v
     b                                           a
```
If `FusionStyle(I)` is `UniqueFusion()` or `SimpleFusion()`, the R-symbol is a
number. Otherwise it is a square matrix with row and column size
`Nsymbol(a, b, c) == Nsymbol(b, a, c)`.
"""
function Rsymbol end

# properties that can be determined in terms of the R symbol

"""
    twist(a::Sector)

Return the twist of a sector `a`.
"""
twist(a::Sector) = sum(dim(b) / dim(a) * tr(Rsymbol(a, a, b)) for b in a ⊗ a)

# Triangle equation
#-------------------------------------------------------------------------------
# requirement that certain F-moves involving unit objects are trivial
function triangle_equation(a::I, b::I; kwargs...) where {I <: Sector}
    for c in ⊗(a, b)
        F1 = Fsymbol(leftunit(a), a, b, c, a, c)
        F2 = Fsymbol(a, rightunit(a), b, c, a, b)
        F3 = Fsymbol(a, b, rightunit(b), c, c, b)

        isapproxone(F) = isapprox(F, one(F); kwargs...)
        if FusionStyle(I) isa MultiplicityFreeFusion
            all(isapproxone, (F1, F2, F3)) || return false
        else
            N = Nsymbol(a, b, c)
            tomatrix = Base.Fix2(reshape, (N, N))
            all(isapproxone ∘ tomatrix, (F1, F2, F3)) || return false
        end
    end
    return true
end

# Pentagon and Hexagon equations
#-------------------------------------------------------------------------------
# Consistency equations for F- and R-symbols
"""
    pentagon_equation(a::I, b::I, c::I, d::I; kwargs...) where {I <: Sector} -> Bool

Check whether the pentagon equation holds for fusing the sectors `a`, `b`, `c` and `d`
to each of their sectors along the two different fusion paths.

If `kwargs` are provided, they are forwarded to `isapprox` when comparing the two sides
of the pentagon equation.
"""
function pentagon_equation(a::I, b::I, c::I, d::I; kwargs...) where {I <: Sector}
    for f in ⊗(a, b), h in ⊗(c, d)
        for g in ⊗(f, c), i in ⊗(b, h)
            for e in intersect(⊗(g, d), ⊗(a, i))
                if FusionStyle(I) isa MultiplicityFreeFusion
                    p1 = Fsymbol(f, c, d, e, g, h) * Fsymbol(a, b, h, e, f, i)
                    p2 = zero(p1)
                    for j in ⊗(b, c)
                        p2 += Fsymbol(a, b, c, g, f, j) *
                            Fsymbol(a, j, d, e, g, i) *
                            Fsymbol(b, c, d, i, j, h)
                    end
                else
                    @tensor p1[λ, μ, ν, κ, ρ, σ] := Fsymbol(f, c, d, e, g, h)[λ, μ, ν, τ] *
                        Fsymbol(a, b, h, e, f, i)[κ, τ, ρ, σ]
                    p2 = zero(p1)
                    for j in ⊗(b, c)
                        @tensor p2[λ, μ, ν, κ, ρ, σ] += Fsymbol(a, b, c, g, f, j)[
                            κ, λ, α,
                            β,
                        ] *
                            Fsymbol(a, j, d, e, g, i)[
                            β, μ, τ,
                            σ,
                        ] *
                            Fsymbol(b, c, d, i, j, h)[
                            α, τ, ν,
                            ρ,
                        ]
                    end
                end
                isapprox(p1, p2; kwargs...) || return false
            end
        end
    end
    return true
end

"""
    hexagon_equation(a::I, b::I, c::I; kwargs...) where {I <: Sector} -> Bool

Check whether the hexagon equation holds for braiding the sector `a` around the fusion
product of `b` and `c` along the two different paths.

If `kwargs` are provided, they are forwarded to `isapprox` when comparing the two sides
of the hexagon equation.
"""
function hexagon_equation(a::I, b::I, c::I; kwargs...) where {I <: Sector}
    BraidingStyle(I) isa NoBraiding &&
        throw(ArgumentError("Hexagon equation only defined for sectors with braiding"))
    for e in ⊗(c, a), f in ⊗(c, b)
        for d in intersect(⊗(e, b), ⊗(a, f))
            if FusionStyle(I) isa MultiplicityFreeFusion
                p1 = Rsymbol(c, a, e) * Fsymbol(a, c, b, d, e, f) * Rsymbol(c, b, f)
                p2 = zero(p1)
                for g in ⊗(a, b)
                    p2 += Fsymbol(c, a, b, d, e, g) * Rsymbol(c, g, d) *
                        Fsymbol(a, b, c, d, g, f)
                end
            else
                @tensor p1[α, β, μ, ν] := Rsymbol(c, a, e)[α, λ] *
                    Fsymbol(a, c, b, d, e, f)[λ, β, γ, ν] * Rsymbol(c, b, f)[γ, μ]
                p2 = zero(p1)
                for g in ⊗(a, b)
                    @tensor p2[α, β, μ, ν] += Fsymbol(c, a, b, d, e, g)[α, β, δ, σ] *
                        Rsymbol(c, g, d)[σ, ψ] * Fsymbol(a, b, c, d, g, f)[δ, ψ, μ, ν]
                end
            end
            isapprox(p1, p2; kwargs...) || return false
        end
    end
    return true
end

# SectorSet:
#-------------------------------------------------------------------------------------------
# Custom generator to represent sets of sectors with type inference
struct SectorSet{I <: Sector, F, S}
    f::F
    set::S
end
SectorSet{I}(::Type{F}, set::S) where {I <: Sector, F, S} = SectorSet{I, Type{F}, S}(F, set)
SectorSet{I}(f::F, set::S) where {I <: Sector, F, S} = SectorSet{I, F, S}(f, set)
SectorSet{I}(set) where {I <: Sector} = SectorSet{I}(identity, set)

Base.IteratorEltype(::Type{<:SectorSet}) = HasEltype()
Base.IteratorSize(::Type{SectorSet{I, F, S}}) where {I <: Sector, F, S} = Base.IteratorSize(S)

Base.eltype(::SectorSet{I}) where {I <: Sector} = I
Base.length(s::SectorSet) = length(s.set)
Base.size(s::SectorSet) = size(s.set)

function Base.iterate(s::SectorSet{I}, args...) where {I <: Sector}
    next = iterate(s.set, args...)
    next === nothing && return nothing
    val, state = next
    return convert(I, s.f(val)), state
end
