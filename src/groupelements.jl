# Sectors corresponding to irreducible representations of compact groups
#------------------------------------------------------------------------------#
"""
    abstract type AbstractGroupElement{G<:Group,ω} <: Sector end

Abstract supertype for sectors which corresponds to group elements of a group `G`, with
the associator given by a 3-cocycle `ω`.

Actual concrete implementations of those irreps can be obtained as `Element[G]`, or via their
actual name, which generically takes the form `(asciiG)Element`, i.e. the ASCII spelling of
the group name followed by `Irrep`.

All group elements have [`FusionStyle`](@ref) equal to `UniqueFusion()`.
For the fusion structure, a specific `SomeGroupElement<:AbstractElment{SomeGroup,ω}`
should only implement the following methods
```julia
Base.:*(c1::GroupElement, c2::GroupElement) -> GroupElement
Base.one(::Type{GroupElement}) -> GroupElement
Base.inv(c::GroupElement) -> GroupElement
TensorKitSectors.cocycle(c1::GroupElement, c2::GroupElement, c3::GroupElement) -> Number
```
and the methods `conj`, `⊗`, `Nsymbol`, `Fsymbol`, `dim`, `Asymbol`, `Bsymbol` and
`frobeniusschur` will be automatically defined.
"""
abstract type AbstractGroupElement{G <: Group, ω} <: Sector end # irreps have integer quantum dimensions
FusionStyle(::Type{<:AbstractGroupElement}) = UniqueFusion()

⊗(a::I, b::I) where {I <: AbstractGroupElement} = (a * b,)
Base.one(a::AbstractGroupElement) = one(typeof(a))
Base.conj(a::AbstractGroupElement) = inv(a)
Nsymbol(a::I, b::I, c::I) where {I <: AbstractGroupElement} = c == a * b
function Fsymbol(a::I, b::I, c::I, d::I, e::I, f::I) where {I <: AbstractGroupElement}
    ω = cocycle(a, b, c)
    if e == a * b && d == e * c && f == b * c
        return ω
    else
        return zero(ω)
    end
end
dim(c::AbstractGroupElement) = 1
frobeniusschur(c::AbstractGroupElement) = cocycle(c, inv(c), c)
function Asymbol(a::I, b::I, c::I) where {I <: AbstractGroupElement}
    A = frobeniusschur(dual(a)) * cocycle(inv(a), a, b)
    return c == a * b ? A : zero(A)
end
function Bsymbol(a::I, b::I, c::I) where {I <: AbstractGroupElement}
    B = cocycle(a, b, inv(b))
    return c == a * b ? B : zero(B)
end

struct ElementTable end
"""
    const Element

A constant of a singleton type used as `Element[G,ω]` with `G<:Group` a type of group, to
construct or obtain a concrete subtype of `AbstractElement{G,ω}` that implements the data
structure used to represent elements of the group `G` with associator given by the
3-cocycle `ω`.
"""
const Element = ElementTable()

type_repr(::Type{<:AbstractGroupElement{G, ω}}) where {G <: Group, ω} = "Element[" * type_repr(G) * ", " * sprint(show, ω) * "]"
function Base.show(io::IO, c::AbstractGroupElement)
    I = typeof(c)
    return if get(io, :typeinfo, nothing) !== I
        print(io, type_repr(I), "(")
        for k in 1:fieldcount(I)
            k > 1 && print(io, ", ")
            print(io, getfield(c, k))
        end
        print(io, ")")
    else
        fieldcount(I) > 1 && print(io, "(")
        for k in 1:fieldcount(I)
            k > 1 && print(io, ", ")
            print(io, getfield(c, k))
        end
        fieldcount(I) > 1 && print(io, ")")
    end
end


# ZNElement: elements of Z_N are labelled by integers mod N; do we ever want N > 64?
"""
    struct ZNElement{N, ω} <: AbstractGroupElement{ℤ{N}, ω}
    ZNElement{N, ω}(n::Integer)
    Element[ℤ{N}, ω](n::Integer)

Represents an element of the group ``ℤ_N`` for some value of `N<64`. (We need 2*(N-1) <= 127 in
order for a ⊗ b to work correctly.) For `N` equals `2`, `3` or `4`, `ℤ{N}` can be replaced
by `ℤ₂`, `ℤ₃`, `ℤ₄`. An arbitrary `Integer` `n` can be provided to the constructor, but only
the value `mod(n, N)` is relevant. The cocycle `ω` should also be specified as an integer 
` 0 <= p < N`, and then leads to the 3-cocycle being given by 

```julia
cocycle(a, b, c) = cis(2* π* p * a.n * (b.n + c.n - mod(b.n + c.n, N)) / N))
```

## Fields
- `n::Int8`: the integer label of the element, modulo `N`.
"""
struct ZNElement{N, p} <: AbstractGroupElement{ℤ{N}, p}
    n::Int8
    function ZNElement{N, p}(n::Integer) where {N, p}
        @assert N < 64
        @assert 0 <= p < N
        return new{N, p}(mod(n, N))
    end
end
ZNElement{N}(n::Integer) where {N} = ZNElement{N, 0}(n)
Base.getindex(::ElementTable, ::Type{ℤ{N}}, p::Int) where {N} = ZNElement{N, mod(p, N)}
Base.convert(::Type{ZNElement{N, p}}, n::Real) where {N, p} = ZNElement{N, p}(n)
const Z2Element{p} = ZNElement{2, p}
const Z3Element{p} = ZNElement{3, p}
const Z4Element{p} = ZNElement{4, p}

Base.one(::Type{ZNElement{N, p}}) where {N, p} = ZNElement{N, p}(0)
Base.inv(c::ZNElement{N, p}) where {N, p} = ZNElement{N, p}(-c.n)
Base.:*(c1::ZNElement{N, p}, c2::ZNElement{N, p}) where {N, p} =
    ZNElement{N, p}(mod(c1.n + c2.n, N))

function cocycle(a::ZNElement{N, p}, b::ZNElement{N, p}, c::ZNElement{N, p}) where {N, p}
    return cis(2 * π * p * a.n * (b.n + c.n - mod(b.n + c.n, N)) / N)
end

Base.IteratorSize(::Type{SectorValues{ZNElement{N, p}}}) where {N, p} = HasLength()
Base.length(::SectorValues{ZNElement{N, p}}) where {N, p} = N
function Base.iterate(::SectorValues{ZNElement{N, p}}, i = 0) where {N, p}
    return i == N ? nothing : (ZNElement{N, p}(i), i + 1)
end
function Base.getindex(::SectorValues{ZNElement{N, p}}, i::Int) where {N, p}
    return 1 <= i <= N ? ZNElement{N, p}(i - 1) : throw(BoundsError(values(ZNElement{N, p}), i))
end
findindex(::SectorValues{ZNElement{N, p}}, c::ZNElement{N, p}) where {N, p} = c.n + 1

Base.hash(c::ZNElement, h::UInt) = hash(c.n, h)
Base.isless(c1::ZNElement{N, p}, c2::ZNElement{N, p}) where {N, p} = isless(c1.n, c2.n)

# TODO: is this true?
# const AbelianGroupElement{G, ω} = AbelianGroupElement{G, ω} where {G <: AbelianGroup}
# BraidingStyle(::Type{<:AbelianGroupElement}) = HasBraiding()

BraidingStyle(::Type{<:AbstractGroupElement}) = NoBraiding()
