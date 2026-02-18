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
    return if BraidingStyle(I) === NoBraiding()
        fusionscalartype(I)
    else
        typeof(zero(fusionscalartype(I)) * zero(braidingscalartype(I)))
    end
end

"""
    fusionscalartype(I::Type{<:Sector}) -> Type{<:Number}

Return the scalar type of the topological data associated to fusion of the sector `I`.
In particular, this is the scalar type of [`Fsymbol`](@ref).

See also [`braidingscalartype`](@ref) and [`sectorscalartype`](@ref).
"""
function fusionscalartype(::Type{I}) where {I <: Sector}
    u = first(allunits(I))
    return eltype(Fsymbol(u, u, u, u, u, u))
end

"""
    braidingscalartype(I::Type{<:Sector}) -> Type{<:Number}

Return the scalar type of the topological data associated to braiding of the sector `I`.
In particular, this is the scalar type of [`Rsymbol`](@ref).

See also [`fusionscalartype`](@ref) and [`sectorscalartype`](@ref).
"""
function braidingscalartype(::Type{I}) where {I <: Sector}
    BraidingStyle(I) === NoBraiding() && throw(ArgumentError("No braiding for sector $I"))
    u = first(allunits(I))
    return eltype(Rsymbol(u, u, u))
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
Each sector `c` should appear at most once in this iteration, even if the multiplicity ``N_c^{ab} > 1``.
The actual multiplicities are accessed separately through [`Nsymbol`](@ref).

The return type is typically [`SectorProductIterator{I}`](@ref) which provides a type-stable iterable that supports pretty-printing, but could also be any custom iterable.

See also [`FusionStyle`](@ref) for the trait associated to the fusion behavior of a given sector type.
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

The fusion multiplicity ``N_c^{ab}``, indicating how many times sector `c` appears in the fusion product `a ⊗ b`.

The return type depends on the [`FusionStyle]`(@ref), where [`UniqueFusion`](@ref) and [`SimpleFusion`](@ref) return `Bool` values, while [`GenericFusion`] returns `Int`.

See also [`⊗`](@ref) to obtain the set of sectors `c` that appear in `a ⊗ b`.
"""
function Nsymbol end

# trait to describe the fusion of superselection sectors
"""
    abstract type FusionStyle
    FusionStyle(::Sector)
    FusionStyle(I::Type{<:Sector})

Trait to describe the fusion behavior of sectors of type `I`, which can be either
* [`UniqueFusion`](@ref): each fusion `a ⊗ b` has exactly one output `c`.
* [`SimpleFusion`](@ref): fusing `a ⊗ b` can lead to multiple values `c`, but each appears at most once.
* [`GenericFusion`](@ref): fusing `a ⊗ b` can lead to multiple values `c` that could appear multiple times.

There is an abstract supertype [`MultipleFusion`](@ref) of which both `SimpleFusion` and `GenericFusion` are subtypes.
Furthermore, there is a type alias [`MultiplicityFreeFusion`](@ref) for those fusion types which do not require muliplicity labels.
"""
abstract type FusionStyle end
FusionStyle(a::Sector) = FusionStyle(typeof(a))

"""
    struct UniqueFusion <: FusionStyle

Fusion style where every product `a ⊗ b` has exactly one output `c`.
As a result, ``N_c^{ab} ≤ 1`` and no multiplicity labels are needed.

See also [`FusionStyle`](@ref).
"""
struct UniqueFusion <: FusionStyle end

"""
    abstract type MultipleFusion <: FusionStyle

Fusion styles that allow more than one fusion output for `a ⊗ b`.

See also [`SimpleFusion`](@ref), [`GenericFusion`](@ref) and [`FusionStyle`](@ref).
"""
abstract type MultipleFusion <: FusionStyle end

"""
    struct SimpleFusion <: MultipleFusion

Fusion style where multiple outputs `c` can appear in `a ⊗ b`, but each appears at most once.
As a result, ``N_c^{ab} ≤ 1`` and no multiplicity labels are needed.

See also [`FusionStyle`](@ref).
"""
struct SimpleFusion <: MultipleFusion end

"""
    struct GenericFusion <: MultipleFusion

Fusion style with potentially multiple outputs `c` and nontrivial multiplicities.
Here ``N_c^{ab}`` can exceed 1, and multiplicity labels are required.

See also [`FusionStyle`](@ref).
"""
struct GenericFusion <: MultipleFusion end

"""
    const MultiplicityFreeFusion = Union{UniqueFusion, SimpleFusion}

Convenience alias for fusion styles that can assume `Nsymbol(a, b, c)::Bool`, and therefore never require multiplicity labels.

See also [`UniqueFusion`](@ref), [`SimpleFusion`](@ref) and [`FusionStyle`](@ref).
"""
const MultiplicityFreeFusion = Union{UniqueFusion, SimpleFusion}

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
* [`SimpleUnit`](@ref): the unit is simple (e.g. fusion categories).
* [`GenericUnit`](@ref): the unit is semisimple (e.g. multifusion categories).
"""
abstract type UnitStyle end
UnitStyle(a::Sector) = UnitStyle(typeof(a))

"""
    struct SimpleUnit <: UnitStyle

Unit style for fusion categories with a unique unit (identity) object.
The unit satisfies ``\\mathbb{1} ⊗ a ≅ a ≅ a ⊗ \\mathbb{1}`` for all sectors.

See also [`UnitStyle`](@ref).
"""
struct SimpleUnit <: UnitStyle end

"""
    struct GenericUnit <: UnitStyle

Unit style for multifusion categories with multiple unit objects (semisimple unit).
Requires implementation of `allunits(::Type{I})`, `leftunit(a)`, and `rightunit(a)`.

See also [`UnitStyle`](@ref).
"""
struct GenericUnit <: UnitStyle end

UnitStyle(::Type{I}) where {I <: Sector} = length(allunits(I)) == 1 ? SimpleUnit() : GenericUnit()

@noinline function throw_genericunit_error(I)
    throw(DomainError(I, "Sector has multiple units, use `allunits` instead of `unit`"))
end

# combine unitstyle properties of tensor products of multifusion sectors
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

function Fsymbol_from_fusiontensor(a::I, b::I, c::I, d::I, e::I, f::I) where {I <: Sector}
    T = fusionscalartype(I)
    Nabe, Necd, Nbcd, Nafd = Nsymbol(a, b, e), Nsymbol(e, c, d), Nsymbol(b, c, f), Nsymbol(a, f, d)
    if iszero(Nabe * Necd * Nbcd * Nafd)
        return FusionStyle(I) isa MultiplicityFreeFusion ? zero(T) : zeros(T, Nabe, Necd, Nbcd, Nafd)
    else
        A = fusiontensor(a, b, e)
        B = @view fusiontensor(e, c, d)[:, :, 1, :]
        C = fusiontensor(b, c, f)
        D = @view fusiontensor(a, f, d)[:, :, 1, :]

        @tensor F[-1, -2, -3, -4] := conj(D[1, 5, -4]) * conj(C[2, 4, 5, -3]) * A[1, 2, 3, -1] * B[3, 4, -2]
        return FusionStyle(I) isa MultiplicityFreeFusion ? only(F) : F
    end
end

# properties that can be determined in terms of the F symbol
# TODO: find mechanism for returning these numbers with custom type T <: AbstractFloat
"""
    dim(a::Sector)

Return the (quantum) dimension of the sector `a`.
"""
dim(a::Sector) = dim_from_Fsymbol(a)

function dim_from_Fsymbol(a::Sector)
    return if FusionStyle(a) isa UniqueFusion
        1
    elseif FusionStyle(a) isa SimpleFusion
        abs(1 / Fsymbol(a, dual(a), a, a, leftunit(a), rightunit(a)))
    else
        abs(1 / Fsymbol(a, dual(a), a, a, leftunit(a), rightunit(a))[1])
    end
end

"""
    sqrtdim(a::Sector)

Return the square root of the (quantum) dimension of sector `a`.

This is a performance specialization that avoids computing `sqrt(1)` for sectors with 
`UniqueFusion`, preserving the number type (returning `1::Int` instead of `1.0::Float64`).
For other sectors, it is equivalent to `sqrt(dim(a))`.
"""
sqrtdim(a::Sector) = (FusionStyle(a) isa UniqueFusion) ? 1 : sqrt(dim(a))

"""
    invsqrtdim(a::Sector)

Return the inverse square root of the (quantum) dimension of sector `a`.

This is a performance specialization that avoids computing `inv(sqrt(1))` for sectors with 
`UniqueFusion`, preserving the number type (returning `1::Int` instead of `1.0::Float64`).
For other sectors, it is equivalent to `inv(sqrt(dim(a)))`.
"""
invsqrtdim(a::Sector) = (FusionStyle(a) isa UniqueFusion) ? 1 : inv(sqrt(dim(a)))

"""
    dimscalartype(::Type{<:Sector}) -> Type{<:Number}

Return the scalar type of the quantum dimensions associated to sectors of type `I`.
In particular, this is the scalar type of [`dim`](@ref).
"""
dimscalartype(::Type{I}) where {I <: Sector} =
    FusionStyle(I) isa UniqueFusion ? Int : typeof(dim(first(allunits(I))))

"""
    frobenius_schur_phase(a::Sector)

Return the Frobenius-Schur phase ``κₐ`` of a sector ``a``, which is a complex phase that
appears in the context of bending lines and is obtained from ``F^{a a̅ a}_a``.
When `a == dual(a)`, it is restricted to ``κₐ ∈ \\{1, -1\\}`` and coincides with
the group-theoretic version [`frobenius_schur_indicator`](@ref).
When `a != dual(a)`, the value of ``κₐ`` can be gauged to be `1`, though is not required to be.
"""
frobenius_schur_phase(a::Sector) = frobenius_schur_phase_from_Fsymbol(a)

function frobenius_schur_phase_from_Fsymbol(a::Sector)
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
    return a == dual(a) ? ν : zero(ν)
end

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
Asymbol(a::I, b::I, c::I) where {I <: Sector} = Asymbol_from_Fsymbol(a, b, c)

function Asymbol_from_Fsymbol(a::I, b::I, c::I) where {I <: Sector}
    return if FusionStyle(I) isa MultiplicityFreeFusion
        (sqrtdim(a) * sqrtdim(b) * invsqrtdim(c)) *
            conj(frobenius_schur_phase(a) * Fsymbol(dual(a), a, b, b, rightunit(a), c))
    else
        reshape(
            (sqrtdim(a) * sqrtdim(b) * invsqrtdim(c)) *
                conj(frobenius_schur_phase(a) * Fsymbol(dual(a), a, b, b, rightunit(a), c)),
            (Nsymbol(a, b, c), Nsymbol(dual(a), c, b))
        )
    end
end
function Asymbol_from_fusiontensor(a::I, b::I, c::I) where {I <: Sector}
    Nabc = Nsymbol(a, b, c)
    T = fusionscalartype(I)
    if Nabc == 0
        return FusionStyle(I) isa MultiplicityFreeFusion ? zero(T) : zeros(T, 0, 0)
    else
        C1 = view(fusiontensor(a, b, c), :, 1, :, :)
        C2 = view(fusiontensor(dual(a), c, b), :, :, 1, :)
        Za = sqrtdim(a) * view(fusiontensor(a, dual(a), leftunit(a)), :, :, 1, 1)
        @tensor A[-1, -2] := sqrtdim(b) / sqrtdim(c) * conj(Za[1, 2]) * C1[1, 3, -1] * C2[2, 3, -2]
        return FusionStyle(I) isa MultiplicityFreeFusion ? only(A) : A
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
Bsymbol(a::I, b::I, c::I) where {I <: Sector} = Bsymbol_from_Fsymbol(a, b, c)

function Bsymbol_from_Fsymbol(a::I, b::I, c::I) where {I <: Sector}
    return if FusionStyle(I) isa MultiplicityFreeFusion
        (sqrtdim(a) * sqrtdim(b) * invsqrtdim(c)) * Fsymbol(a, b, dual(b), a, c, rightunit(a))
    else
        reshape(
            (sqrtdim(a) * sqrtdim(b) * invsqrtdim(c)) * Fsymbol(a, b, dual(b), a, c, rightunit(a)),
            (Nsymbol(a, b, c), Nsymbol(c, dual(b), a))
        )
    end
end
function Bsymbol_from_fusiontensor(a::I, b::I, c::I) where {I <: Sector}
    Nabc = Nsymbol(a, b, c)
    T = fusionscalartype(I)
    if Nabc == 0
        return FusionStyle(I) isa MultiplicityFreeFusion ? zero(T) : zeros(T, 0, 0)
    else
        C1 = view(fusiontensor(a, b, c), 1, :, :, :)
        C2 = view(fusiontensor(c, dual(b), a), :, :, 1, :)
        Zb = sqrtdim(b) * view(fusiontensor(b, dual(b), leftunit(b)), :, :, 1, 1)
        @tensor B[-1, -2] := sqrtdim(a) / sqrtdim(c) * conj(Zb[1, 2]) * C1[1, 3, -1] * C2[3, 2, -2]
        return FusionStyle(I) isa MultiplicityFreeFusion ? only(B) : B
    end
end

# Braiding:
#-------------------------------------------------
# trait to describe type to denote how the elementary spaces in a tensor product space
# interact under permutations or actions of the braid group
"""
    abstract type BraidingStyle
    BraidingStyle(::Sector) -> ::BraidingStyle
    BraidingStyle(I::Type{<:Sector}) -> ::BraidingStyle

Trait to describe the braiding behavior of sectors of type `I`, which can be either
* [`NoBraiding`](@ref): no braiding structure defined.
* [`Bosonic`](@ref): symmetric braiding structure with a trivial twist.
* [`Fermionic`](@ref): symmetric braiding structure with a non-trivial twist that squares to identity.
* [`Anyonic`](@ref): general braiding structure and arbitrary twists.

There is an abstract supertype [`HasBraiding`](@ref) that includes all styles that define [`Rsymbol`](@ref) (everything but `NoBraiding`).
Furthermore, the abstract supertype [`SymmetricBraiding`](@ref) denotes the cases where braidings are equivalent to crossings, i.e. braiding twice is an identity operation.
This includes the `Bosonic` and `Fermionic` styles, for which we can uniquely define permutations.
"""
abstract type BraidingStyle end
BraidingStyle(a::Sector) = BraidingStyle(typeof(a))

"""
    abstract type HasBraiding <: BraidingStyle

Supertype for all braiding styles where an [`Rsymbol`](@ref) is defined.
This includes all current `BraidingStyle`s except `NoBraiding`.
"""
abstract type HasBraiding <: BraidingStyle end

"""
    struct NoBraiding <: BraidingStyle

Braiding style for categories without a braiding structure.
Except for braiding with the unit sector, only planar diagrams are meaningful; [`Rsymbol`](@ref) is undefined.

See also [`BraidingStyle`](@ref).
"""
struct NoBraiding <: BraidingStyle end

"""
    abstract type SymmetricBraiding <: HasBraiding

Supertype for braiding styles with symmetric braiding, where braiding twice is the identity operation.
Subtypes include [`Bosonic`](@ref) (trivial twist) and [`Fermionic`](@ref) (nontrivial twist ±1).
Supports permutation group statistics.

See also [`BraidingStyle`](@ref).
"""
abstract type SymmetricBraiding <: HasBraiding end

"""
    struct Bosonic <: SymmetricBraiding

Braiding style with symmetric braiding and trivial twist.
This is characterized by ``R^{ab}_c R^{ba}_c = 1`` and ``\\theta_a = 1`` for all sectors.

See also [`BraidingStyle`](@ref).
"""
struct Bosonic <: SymmetricBraiding end

"""
    struct Fermionic <: SymmetricBraiding

Braiding style with symmetric braiding and nontrivial (symmetric) twist.
This is characterized by ``R^{ab}_c R^{ba}_c = 1`` and ``\\theta_a = \\pm 1`` for all sectors.

See also [`BraidingStyle`](@ref).
"""
struct Fermionic <: SymmetricBraiding end

"""
    struct Anyonic <: HasBraiding

Braiding style with general (non-symmetric) braiding and arbitrary twists.
Characterized by nontrivial braid group representations where ``R^{ab}_c R^{ba}_c ≠ 1`` in general.

See also [`BraidingStyle`](@ref).
"""
struct Anyonic <: HasBraiding end

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

function Rsymbol_from_fusiontensor(a::I, b::I, c::I) where {I <: Sector}
    Nabc = Nsymbol(a, b, c)
    T = braidingscalartype(I)
    if Nabc == 0
        return FusionStyle(I) isa MultiplicityFreeFusion ? zero(T) : zeros(T, 0, 0)
    else
        A = view(fusiontensor(a, b, c), :, :, 1, :)
        B = view(fusiontensor(b, a, c), :, :, 1, :)
        @tensor R[-1 -2] := conj(B[1 2 -2]) * A[2 1 -1]
        return FusionStyle(I) isa MultiplicityFreeFusion ? only(R) : R
    end
end

# properties that can be determined in terms of the R symbol

"""
    twist(a::Sector)

Return the twist of a sector `a`.
"""
twist(a::Sector) = twist_from_Rsymbol(a)
twist_from_Rsymbol(a::Sector) = sum(dim(b) / dim(a) * tr(Rsymbol(a, a, b)) for b in a ⊗ a)

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
