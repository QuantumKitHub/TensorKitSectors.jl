"""
    abstract type Sector end

Abstract type for representing the (isomorphism classes of) simple objects in (unitary
and pivotal) (pre-)fusion categories, e.g. the irreducible representations of a finite or
compact group. Subtypes `I<:Sector` as the set of labels of a `GradedSpace`.

Every new `I<:Sector` should implement the following methods:
*   `one(::Type{I})`: unit element of `I`
*   `conj(a::I)`: ``a̅``, conjugate or dual label of ``a``
*   `⊗(a::I, b::I)`: iterable with unique fusion outputs of ``a ⊗ b``
    (i.e. don't repeat in case of multiplicities)
*   `Nsymbol(a::I, b::I, c::I)`: number of times `c` appears in `a ⊗ b`, i.e. the
    multiplicity
*   `FusionStyle(::Type{I})`: `UniqueFusion()`, `SimpleFusion()` or
    `GenericFusion()`
*   `BraidingStyle(::Type{I})`: `Bosonic()`, `Fermionic()`, `Anyonic()`, ...
*   `Fsymbol(a::I, b::I, c::I, d::I, e::I, f::I)`: F-symbol: scalar (in case of
    `UniqueFusion`/`SimpleFusion`) or matrix (in case of `GenericFusion`)
*   `Rsymbol(a::I, b::I, c::I)`: R-symbol: scalar (in case of
    `UniqueFusion`/`SimpleFusion`) or matrix (in case of `GenericFusion`)
and optionally
*   `dim(a::I)`: quantum dimension of sector `a`
*   `frobeniusschur(a::I)`: Frobenius-Schur indicator of `a`
*   `Bsymbol(a::I, b::I, c::I)`: B-symbol: scalar (in case of
    `UniqueFusion`/`SimpleFusion`) or matrix (in case of `GenericFusion`)
*   `twist(a::I)` -> twist of sector `a`

Furthermore, `iterate` and `Base.IteratorSize` should be made to work for the singleton type
[`SectorValues{I}`](@ref).
"""
abstract type Sector end

# iterator over the values (i.e., elements of representative set of simple objects)
# in the sector
"""
    struct SectorValues{I<:Sector}

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
    one(::Sector) -> Sector
    one(::Type{<:Sector}) -> Sector

Return the unit element within this type of sector.
"""
Base.one(a::Sector) = one(typeof(a))

"""
    leftone(a::Sector) -> Sector

Return the left unit element within this type of sector.
See also [`rightone`](@ref) and [`Base.one`](@ref).
"""
leftone(a::Sector) = one(a)

"""
    rightone(a::Sector) -> Sector

Return the right unit element within this type of sector.
See also [`leftone`](@ref) and [`Base.one`](@ref).
"""
rightone(a::Sector) = one(a)

"""
    allones(I::Type{<:Sector}) -> Tuple{I}

Return a tuple with all units of the sector type `I`.
For fusion categories, this will contain only one element.
"""
allones(I::Type{<:Sector}) = (one(I),)

"""
    dual(a::Sector) -> Sector

Return the conjugate label `conj(a)`.
"""
dual(a::Sector) = conj(a)

"""
    sectorscalartype(I::Type{<:Sector}) -> Type

Return the scalar type of the topological data (Fsymbol, Rsymbol) of the sector `I`.
"""
function sectorscalartype(::Type{I}) where {I <: Sector}
    if BraidingStyle(I) === NoBraiding()
        return _Fscalartype(I)
    else
        return Base.promote_op(*, _Fscalartype(I), _Rscalartype(I))
    end
end
function _Fscalartype(::Type{I}) where {I <: Sector}
    Ftype = Core.Compiler.return_type(Fsymbol, NTuple{6, I})
    return FusionStyle(I) === UniqueFusion() ? Ftype : eltype(Ftype)
end
function _Rscalartype(::Type{I}) where {I <: Sector}
    BraidingStyle(I) === NoBraiding() && throw(ArgumentError("No braiding for sector $I"))
    Rtype = Core.Compiler.return_type(Rsymbol, NTuple{3, I})
    return FusionStyle(I) === UniqueFusion() ? Rtype : eltype(Rtype)
end

"""
    isreal(::Type{<:Sector}) -> Bool

Return whether the topological data (Fsymbol, Rsymbol) of the sector is real or not (in
which case it is complex).
"""
Base.isreal(I::Type{<:Sector}) = sectorscalartype(I) <: Real

# FusionStyle: the most important aspect of Sector
#---------------------------------------------
"""
    ⊗(a::I, b::I...) where {I<:Sector}
    otimes(a::I, b::I...) where {I<:Sector}

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
    Nsymbol(a::I, b::I, c::I) where {I<:Sector} -> Integer

Return an `Integer` representing the number of times `c` appears in the fusion product
`a ⊗ b`. Could be a `Bool` if `FusionStyle(I) == UniqueFusion()` or `SimpleFusion()`.
"""
function Nsymbol end

# trait to describe the fusion of superselection sectors
"""
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

# combine fusion properties of tensor products of sectors
Base.:&(f::F, ::F) where {F <: FusionStyle} = f
Base.:&(f₁::FusionStyle, f₂::FusionStyle) = f₂ & f₁

Base.:&(::SimpleFusion, ::UniqueFusion) = SimpleFusion()
Base.:&(::GenericFusion, ::UniqueFusion) = GenericFusion()
Base.:&(::GenericFusion, ::SimpleFusion) = GenericFusion()

# similar, but for multifusion categories
"""
    MultiFusionStyle(::Sector)
    MultiFusionStyle(I::Type{<:Sector})

Trait to describe the semisimplicity of the unit sector of type `I`.
This can be either
*   `SimpleMultiFusion()`: the unit is simple (e.g. fusion categories);
*   `GenericMultiFusion()`: the unit is semisimple.
"""
abstract type MultiFusionStyle end
MultiFusionStyle(a::Sector) = MultiFusionStyle(typeof(a))

struct SimpleMultiFusion <: MultiFusionStyle end
struct GenericMultiFusion <: MultiFusionStyle end

MultiFusionStyle(::Type{<:Sector}) = SimpleMultiFusion() # default is simple unit

# combine fusion properties of tensor products of multifusion sectors
Base.:&(f::F, ::F) where {F <: MultiFusionStyle} = f
Base.:&(f₁::MultiFusionStyle, f₂::MultiFusionStyle) = f₂ & f₁

Base.:&(::GenericMultiFusion, ::SimpleMultiFusion) = GenericMultiFusion()

"""
    Fsymbol(a::I, b::I, c::I, d::I, e::I, f::I) where {I<:Sector}

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
# TODO: find mechanism for returning these numbers with custom type T<:AbstractFloat
"""
    dim(a::Sector)

Return the (quantum) dimension of the sector `a`.
"""
function dim(a::Sector)
    return if FusionStyle(a) isa UniqueFusion
        1
    elseif FusionStyle(a) isa SimpleFusion
        abs(1 / Fsymbol(a, conj(a), a, a, leftone(a), rightone(a)))
    else
        abs(1 / Fsymbol(a, conj(a), a, a, leftone(a), rightone(a))[1])
    end
end
sqrtdim(a::Sector) = (FusionStyle(a) isa UniqueFusion) ? 1 : sqrt(dim(a))
invsqrtdim(a::Sector) = (FusionStyle(a) isa UniqueFusion) ? 1 : inv(sqrt(dim(a)))

"""
    frobeniusschur(a::Sector)

Return the Frobenius-Schur indicator of a sector `a`.
"""
function frobeniusschur(a::Sector)
    return if FusionStyle(a) isa UniqueFusion || FusionStyle(a) isa SimpleFusion
        sign(Fsymbol(a, conj(a), a, a, leftone(a), rightone(a)))
    else
        sign(Fsymbol(a, conj(a), a, a, leftone(a), rightone(a))[1])
    end
end

# Not necessary
function Asymbol(a::I, b::I, c::I) where {I <: Sector}
    return if FusionStyle(I) isa UniqueFusion || FusionStyle(I) isa SimpleFusion
        (sqrtdim(a) * sqrtdim(b) * invsqrtdim(c)) *
            conj(frobeniusschur(a) * Fsymbol(dual(a), a, b, b, leftone(a), c))
    else
        reshape(
            (sqrtdim(a) * sqrtdim(b) * invsqrtdim(c)) *
                conj(frobeniusschur(a) * Fsymbol(dual(a), a, b, b, leftone(a), c)),
            (Nsymbol(a, b, c), Nsymbol(dual(a), c, b))
        )
    end
end

"""
    Bsymbol(a::I, b::I, c::I) where {I<:Sector}

Return the value of ``B^{ab}_c`` which appears in transforming a splitting vertex
into a fusion vertex using the transformation
```
a -<-μ-<- c                                                    a -<-ν-<- c
     ∨          -> √(dim(c)/dim(a)) * Bsymbol(a,b,c)[μ,ν]           ∧
     b                                                            dual(b)
```
If `FusionStyle(I)` is `UniqueFusion()` or `SimpleFusion()`, the B-symbol is a
number. Otherwise it is a square matrix with row and column size
`Nsymbol(a, b, c) == Nsymbol(c, dual(b), a)`.
"""
function Bsymbol(a::I, b::I, c::I) where {I <: Sector}
    return if FusionStyle(I) isa UniqueFusion || FusionStyle(I) isa SimpleFusion
        (sqrtdim(a) * sqrtdim(b) * invsqrtdim(c)) * Fsymbol(a, b, dual(b), a, c, rightone(a))
    else
        reshape(
            (sqrtdim(a) * sqrtdim(b) * invsqrtdim(c)) * Fsymbol(a, b, dual(b), a, c, rightone(a)),
            (Nsymbol(a, b, c), Nsymbol(c, dual(b), a))
        )
    end
end

# Braiding:
#-------------------------------------------------
# trait to describe type to denote how the elementary spaces in a tensor product space
# interact under permutations or actions of the braid group
"""
    BraidingStyle(::Sector) -> ::BraidingStyle
    BraidingStyle(I::Type{<:Sector}) -> ::BraidingStyle

Return the type of braiding and twist behavior of sectors of type `I`, which can be either
*   `Bosonic()`: symmetric braiding with trivial twist (i.e. identity)
*   `Fermionic()`: symmetric braiding with non-trivial twist (squares to identity)
*   `Anyonic()`: general ``R_(a,b)^c`` phase or matrix (depending on `SimpleFusion` or
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

Base.:&(b::B, ::B) where {B <: BraidingStyle} = b
Base.:&(B1::BraidingStyle, B2::BraidingStyle) = B2 & B1
Base.:&(::Bosonic, ::Fermionic) = Fermionic()
Base.:&(::Bosonic, ::Anyonic) = Anyonic()
Base.:&(::Fermionic, ::Anyonic) = Anyonic()
Base.:&(::Bosonic, ::NoBraiding) = NoBraiding()
Base.:&(::Fermionic, ::NoBraiding) = NoBraiding()
Base.:&(::Anyonic, ::NoBraiding) = NoBraiding()

"""
    Rsymbol(a::I, b::I, c::I) where {I<:Sector}

Returns the R-symbol ``R^{ab}_c`` that maps between ``c → a ⊗ b`` and ``c → b ⊗ a`` as in
```
a -<-μ-<- c                                 b -<-ν-<- c
     ∨          -> Rsymbol(a,b,c)[μ,ν]           v
     b                                           a
```
If `FusionStyle(I)` is `UniqueFusion()` or `SimpleFusion()`, the R-symbol is a
number. Otherwise it is a square matrix with row and column size
`Nsymbol(a,b,c) == Nsymbol(b,a,c)`.
"""
function Rsymbol end

# properties that can be determined in terms of the R symbol

"""
    twist(a::Sector)

Return the twist of a sector `a`
"""
twist(a::Sector) = sum(dim(b) / dim(a) * tr(Rsymbol(a, a, b)) for b in a ⊗ a)

# Triangle equation
#-------------------------------------------------------------------------------
# requirement that certain F-moves involving unit objects are trivial
function triangle_equation(a::I, b::I; kwargs...) where {I <: Sector}
    for c in ⊗(a, b)
        F1 = Fsymbol(leftone(a), a, b, c, a, c)
        F2 = Fsymbol(a, rightone(a), b, c, a, b)
        F3 = Fsymbol(a, b, rightone(b), c, c, b)

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
#-------------------------------------------------------------------------------
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

# Time reversed sector
struct TimeReversed{I <: Sector} <: Sector
    a::I
    function TimeReversed{I}(a::I) where {I <: Sector}
        if BraidingStyle(I) isa NoBraiding
            throw(ArgumentError("TimeReversed is not defined for sectors $I with no braiding"))
        end
        return new{I}(a)
    end
end
FusionStyle(::Type{TimeReversed{I}}) where {I <: Sector} = FusionStyle(I)
BraidingStyle(::Type{TimeReversed{I}}) where {I <: Sector} = BraidingStyle(I)
function Nsymbol(
        a::TimeReversed{I}, b::TimeReversed{I}, c::TimeReversed{I}
    ) where {I <: Sector}
    return Nsymbol(a.a, b.a, c.a)
end

function Fsymbol(
        a::TimeReversed{I}, b::TimeReversed{I}, c::TimeReversed{I},
        d::TimeReversed{I}, e::TimeReversed{I}, f::TimeReversed{I}
    ) where {I <: Sector}
    return Fsymbol(a.a, b.a, c.a, d.a, e.a, f.a)
end
function Rsymbol(
        a::TimeReversed{I}, b::TimeReversed{I}, c::TimeReversed{I}
    ) where {I <: Sector}
    return adjoint(Rsymbol(a.a, b.a, c.a))
end

Base.one(::Type{TimeReversed{I}}) where {I <: Sector} = TimeReversed{I}(one(I))
Base.conj(c::TimeReversed{I}) where {I <: Sector} = TimeReversed{I}(conj(c.a))
function ⊗(c1::TimeReversed{I}, c2::TimeReversed{I}) where {I <: Sector}
    return Iterators.map(TimeReversed{I}, c1.a ⊗ c2.a)
end
function Base.IteratorSize(::Type{SectorValues{TimeReversed{I}}}) where {I <: Sector}
    return Base.IteratorSize(values(I))
end
function Base.length(::SectorValues{TimeReversed{I}}) where {I <: Sector}
    return length(values(I))
end
function Base.getindex(::SectorValues{TimeReversed{I}}, i::Int) where {I <: Sector}
    return TimeReversed{I}(getindex(values(I), i))
end
function Base.iterate(::SectorValues{TimeReversed{I}}, state...) where {I <: Sector}
    next = iterate(values(I), state...)
    if isnothing(next)
        return nothing
    else
        obj, nextstate = next
        return TimeReversed{I}(obj), nextstate
    end
end
function findindex(
        ::SectorValues{TimeReversed{I}}, a::TimeReversed{I}
    ) where {I <: Sector}
    return findindex(values(I), a.a)
end

function Base.isless(c1::TimeReversed{I}, c2::TimeReversed{I}) where {I <: Sector}
    return isless(c1.a, c2.a)
end
