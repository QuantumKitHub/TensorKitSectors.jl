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
and optionally, if `FusionStyle(I) isa GenericFusion`
*   `vertex_ind2label(i::Int, a::I, b::I, c::I)` -> a custom label for the `i`th copy of
    `c` appearing in `a ⊗ b`

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
*   `Base.IteratorSize(::Type{SectorValues{I}})`: `HasLenght()`, `SizeUnkown()`
    or `IsInfinite()` depending on whether the number of values of type `I` is finite
    (and sufficiently small) or infinite; for a large number of values, `SizeUnknown()` is
    recommend because this will trigger the use of `GenericGradedSpace`.
If `IteratorSize(I) == HasLength()`, also the following must be implemented:
*   `Base.length(::SectorValues{I})`: the number of different values
*   `Base.getindex(::SectorValues{I}, i::Int)`: a mapping between an index `i` and an
    instance of `I`
*   `findindex(::SectorValues{I}, c::I)`: reverse mapping between a value `c::I` and an
    index `i::Integer ∈ 1:length(values(I))`
"""
struct SectorValues{I<:Sector} end
Base.IteratorEltype(::Type{<:SectorValues}) = HasEltype()
Base.eltype(::Type{SectorValues{I}}) where {I<:Sector} = I
Base.values(::Type{I}) where {I<:Sector} = SectorValues{I}()

"""
    one(::Sector) -> Sector
    one(::Type{<:Sector}) -> Sector

Return the unit element within this type of sector.
"""
Base.one(a::Sector) = one(typeof(a))

"""
    dual(a::Sector) -> Sector

Return the conjugate label `conj(a)`.
"""
dual(a::Sector) = conj(a)

"""
    sectorscalartype(I::Type{<:Sector}) -> Type

Return the scalar type of the topological data (Fsymbol, Rsymbol) of the sector `I`.
"""
function sectorscalartype(::Type{I}) where {I<:Sector}
    if BraidingStyle(I) isa NoBraiding
        return eltype(Core.Compiler.return_type(Fsymbol, NTuple{6,I}))
    else
        Feltype = eltype(Core.Compiler.return_type(Fsymbol, NTuple{6,I}))
        Reltype = eltype(Core.Compiler.return_type(Rsymbol, NTuple{3,I}))
        return Base.promote_op(*, Feltype, Reltype)
    end
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
@inline function ⊗(a::I, b::I, c::I, rest::Vararg{I}) where {I<:Sector}
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
const MultiplicityFreeFusion = Union{UniqueFusion,SimpleFusion}

# combine fusion properties of tensor products of sectors
Base.:&(f::F, ::F) where {F<:FusionStyle} = f
Base.:&(f₁::FusionStyle, f₂::FusionStyle) = f₂ & f₁

Base.:&(::SimpleFusion, ::UniqueFusion) = SimpleFusion()
Base.:&(::GenericFusion, ::UniqueFusion) = GenericFusion()
Base.:&(::GenericFusion, ::SimpleFusion) = GenericFusion()

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

# If a I::Sector with `fusion(I) == GenericFusion` fusion wants to have custom vertex
# labels, a specialized method for `vertindex2label` should be added
"""
    vertex_ind2label(k::Int, a::I, b::I, c::I) where {I<:Sector}

Convert the index `k` of the fusion vertex (a,b)->c into a label. For
`FusionStyle(I) == UniqueFusion()` or `FusionStyle(I) == MultipleFusion()`, where every
fusion output occurs only once and `k == 1`, the default is to suppress vertex labels by
setting them equal to `nothing`. For `FusionStyle(I) == GenericFusion()`, the default is to
just use `k`, unless a specialized method is provided.
"""
function vertex_ind2label(k::Int, a::I, b::I, c::I) where {I<:Sector}
    return _ind2label(FusionStyle(I), k::Int, a::I, b::I, c::I)
end
_ind2label(::UniqueFusion, k, a, b, c) = nothing
_ind2label(::SimpleFusion, k, a, b, c) = nothing
_ind2label(::GenericFusion, k, a, b, c) = k

"""
    vertex_labeltype(I::Type{<:Sector}) -> Type

Return the type of labels for the fusion vertices of sectors of type `I`.
"""
vertex_labeltype(I::Type{<:Sector}) = typeof(vertex_ind2label(1, one(I), one(I), one(I)))

# properties that can be determined in terms of the F symbol
# TODO: find mechanism for returning these numbers with custom type T<:AbstractFloat
"""
    dim(a::Sector)

Return the (quantum) dimension of the sector `a`.
"""
function dim(a::Sector)
    if FusionStyle(a) isa UniqueFusion
        1
    elseif FusionStyle(a) isa SimpleFusion
        abs(1 / Fsymbol(a, conj(a), a, a, one(a), one(a)))
    else
        abs(1 / Fsymbol(a, conj(a), a, a, one(a), one(a))[1])
    end
end
sqrtdim(a::Sector) = (FusionStyle(a) isa UniqueFusion) ? 1 : sqrt(dim(a))
invsqrtdim(a::Sector) = (FusionStyle(a) isa UniqueFusion) ? 1 : inv(sqrt(dim(a)))

"""
    frobeniusschur(a::Sector)

Return the Frobenius-Schur indicator of a sector `a`.
"""
function frobeniusschur(a::Sector)
    if FusionStyle(a) isa UniqueFusion || FusionStyle(a) isa SimpleFusion
        sign(Fsymbol(a, conj(a), a, a, one(a), one(a)))
    else
        sign(Fsymbol(a, conj(a), a, a, one(a), one(a))[1])
    end
end

# Not necessary
function Asymbol(a::I, b::I, c::I) where {I<:Sector}
    if FusionStyle(I) isa UniqueFusion || FusionStyle(I) isa SimpleFusion
        (sqrtdim(a) * sqrtdim(b) * invsqrtdim(c)) *
        conj(frobeniusschur(a) * Fsymbol(dual(a), a, b, b, one(a), c))
    else
        reshape((sqrtdim(a) * sqrtdim(b) * invsqrtdim(c)) *
                conj(frobeniusschur(a) * Fsymbol(dual(a), a, b, b, one(a), c)),
                (Nsymbol(a, b, c), Nsymbol(dual(a), c, b)))
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
function Bsymbol(a::I, b::I, c::I) where {I<:Sector}
    if FusionStyle(I) isa UniqueFusion || FusionStyle(I) isa SimpleFusion
        (sqrtdim(a) * sqrtdim(b) * invsqrtdim(c)) * Fsymbol(a, b, dual(b), a, c, one(a))
    else
        reshape((sqrtdim(a) * sqrtdim(b) * invsqrtdim(c)) *
                Fsymbol(a, b, dual(b), a, c, one(a)),
                (Nsymbol(a, b, c), Nsymbol(c, dual(b), a)))
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

Base.:&(b::B, ::B) where {B<:BraidingStyle} = b
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

# Pentagon and Hexagon equations
#-------------------------------------------------------------------------------
# Consistency equations for F- and R-symbols
function pentagon_equation(a::I, b::I, c::I, d::I; kwargs...) where {I<:Sector}
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
                        @tensor p2[λ, μ, ν, κ, ρ, σ] += Fsymbol(a, b, c, g, f, j)[κ, λ, α,
                                                                                  β] *
                                                        Fsymbol(a, j, d, e, g, i)[β, μ, τ,
                                                                                  σ] *
                                                        Fsymbol(b, c, d, i, j, h)[α, τ, ν,
                                                                                  ρ]
                    end
                end
                isapprox(p1, p2; kwargs...) || return false
            end
        end
    end
    return true
end

function hexagon_equation(a::I, b::I, c::I; kwargs...) where {I<:Sector}
    BraidingStyle(I) isa NoBraiding &&
        throw(ArgumentError("Hexagon equation only defined for sectors with braiding"))
    for e in ⊗(c, a), f in ⊗(c, b)
        for d in intersect(⊗(e, b), ⊗(a, f))
            if FusionStyle(I) isa MultiplicityFreeFusion
                p1 = Rsymbol(a, c, e) * Fsymbol(a, c, b, d, e, f) * Rsymbol(b, c, f)
                p2 = zero(p1)
                for g in ⊗(a, b)
                    p2 += Fsymbol(c, a, b, d, e, g) * Rsymbol(g, c, d) *
                          Fsymbol(a, b, c, d, g, f)
                end
            else
                @tensor p1[α, β, μ, ν] := Rsymbol(a, c, e)[α, λ] *
                                          Fsymbol(a, c, b, d, e, f)[λ, β, γ, ν] *
                                          Rsymbol(b, c, f)[γ, μ]
                p2 = zero(p1)
                for g in ⊗(a, b)
                    @tensor p2[α, β, μ, ν] += Fsymbol(c, a, b, d, e, g)[α, β, δ,
                                                                        σ] *
                                              Rsymbol(g, c, d)[σ, ψ] *
                                              Fsymbol(a, b, c, d, g, f)[δ, ψ, μ, ν]
                end
            end
            isapprox(p1, p2; kwargs...) || return false
        end
    end
    return true
end

# SectorSet:
#-------------------------------------------------------------------------------
# Custum generator to represent sets of sectors with type inference
struct SectorSet{I<:Sector,F,S}
    f::F
    set::S
end
SectorSet{I}(::Type{F}, set::S) where {I<:Sector,F,S} = SectorSet{I,Type{F},S}(F, set)
SectorSet{I}(f::F, set::S) where {I<:Sector,F,S} = SectorSet{I,F,S}(f, set)
SectorSet{I}(set) where {I<:Sector} = SectorSet{I}(identity, set)

Base.IteratorEltype(::Type{<:SectorSet}) = HasEltype()
Base.IteratorSize(::Type{SectorSet{I,F,S}}) where {I<:Sector,F,S} = Base.IteratorSize(S)

Base.eltype(::SectorSet{I}) where {I<:Sector} = I
Base.length(s::SectorSet) = length(s.set)
Base.size(s::SectorSet) = size(s.set)

function Base.iterate(s::SectorSet{I}, args...) where {I<:Sector}
    next = iterate(s.set, args...)
    next === nothing && return nothing
    val, state = next
    return convert(I, s.f(val)), state
end