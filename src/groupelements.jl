# Sectors corresponding to the elements of a group, possibly with a nontrivial
# associator given by a 3-cocycle.
#------------------------------------------------------------------------------#
"""
    abstract type AbstractGroupElement{G<:Group} <: Sector end

Abstract supertype for sectors which corresponds to group elements of a group `G`.

Actual concrete implementations of those irreps can be obtained as `Element[G]`, or via their
actual name, which generically takes the form `(asciiG)Element`, i.e. the ASCII spelling of
the group name followed by `Element`.

All group elements have [`FusionStyle`](@ref) equal to `UniqueFusion()`.
Furthermore, the [`BraidingStyle`](@ref) is set to `NoBraiding()`, although this can be
overridden by a concrete implementation of `AbstractGroupElement`.

For the fusion structure, a specific `SomeGroupElement<:AbstractGroupElement{SomeGroup}`
should only implement the following methods
```julia
Base.:*(c1::GroupElement, c2::GroupElement) -> GroupElement
Base.one(::Type{GroupElement}) -> GroupElement
Base.inv(c::GroupElement) -> GroupElement
# and optionally
TensorKitSectors.cocycle(c1::GroupElement, c2::GroupElement, c3::GroupElement) -> Number
```
The methods `conj`, `dual`, `⊗`, `Nsymbol`, `Fsymbol`, `dim`, `Asymbol`, `Bsymbol` and
`frobeniusschur` will then be automatically defined. If no `cocycle` method is defined,
the cocycle will be assumed to be trivial, i.e. equal to `1`.

"""
abstract type AbstractGroupElement{G <: Group} <: Sector end
FusionStyle(::Type{<:AbstractGroupElement}) = UniqueFusion()
BraidingStyle(::Type{<:AbstractGroupElement}) = NoBraiding()

cocycle(a::I, b::I, c::I) where {I <: AbstractGroupElement} = 1
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
    const GroupElement

A constant of a singleton type used as `Element[G]` or `Element[G,ω]` with `G<:Group`
a type of group, to construct or obtain a concrete subtype of `AbstractElement{G}`
that implements the data structure used to represent elements of the group `G`, possibly
with a second argument `ω` that specifies the associated 3-cocycle.
"""
const GroupElement = ElementTable()

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
    struct ZNElement{N, p} <: AbstractGroupElement{ℤ{N}}
    ZNElement{N, p}(n::Integer)
    GroupElement[ℤ{N}, p](n::Integer)

Represents an element of the group ``ℤ_N`` for some value of `N<64`. (We need 2*(N-1) <= 127 in
order for a ⊗ b to work correctly.) For `N` equals `2`, `3` or `4`, `ℤ{N}` can be replaced
by `ℤ₂`, `ℤ₃`, `ℤ₄`. An arbitrary `Integer` `n` can be provided to the constructor, but only
the value `mod(n, N)` is relevant. The second type parameter `p` should also be specified as
an integer ` 0 <= p < N` and specifies the 3-cocycle, which is then being given by 

```julia
cocycle(a, b, c) = cispi(2 * p * a.n * (b.n + c.n - mod(b.n + c.n, N)) / N))
```

If `p` is not specified, it defaults to `0`, i.e. the trivial cocycle.

## Fields
- `n::Int8`: the integer label of the element, modulo `N`.
"""
struct ZNElement{N, p} <: AbstractGroupElement{ℤ{N}}
    n::Int8
    function ZNElement{N, p}(n::Integer) where {N, p}
        @assert N < 64
        @assert 0 <= p < N
        return new{N, p}(mod(n, N))
    end
end
ZNElement{N}(n::Integer) where {N} = ZNElement{N, 0}(n)
Base.getindex(::ElementTable, ::Type{ℤ{N}}, p::Int = 0) where {N} = ZNElement{N, mod(p, N)}
type_repr(::Type{ZNElement{N, p}}) where {N, p} = "GroupElement[ℤ{$N}, $p]"

Base.convert(T::Type{<:ZNElement}, n::Real) = T(n)
const Z2Element{p} = ZNElement{2, p}
const Z3Element{p} = ZNElement{3, p}
const Z4Element{p} = ZNElement{4, p}

Base.one(::Type{Z}) where {Z <: ZNElement} = Z(0)
Base.inv(c::ZNElement) = typeof(c)(-c.n)
Base.:*(c1::ZNElement{N, p}, c2::ZNElement{N, p}) where {N, p} =
    ZNElement{N, p}(mod(c1.n + c2.n, N))

function cocycle(a::ZNElement{N, p}, b::ZNElement{N, p}, c::ZNElement{N, p}) where {N, p}
    return cispi(2 * p * a.n * (b.n + c.n - mod(b.n + c.n, N)) / N)
end

Base.IteratorSize(::Type{SectorValues{<:ZNElement}}) = HasLength()
Base.length(::SectorValues{<:ZNElement{N}}) where {N} = N
function Base.iterate(::SectorValues{ZNElement{N, p}}, i = 0) where {N, p}
    return i == N ? nothing : (ZNElement{N, p}(i), i + 1)
end
function Base.getindex(::SectorValues{ZNElement{N, p}}, i::Int) where {N, p}
    return 1 <= i <= N ? ZNElement{N, p}(i - 1) : throw(BoundsError(values(ZNElement{N, p}), i))
end
findindex(::SectorValues{ZNElement{N, p}}, c::ZNElement{N, p}) where {N, p} = c.n + 1

Base.hash(c::ZNElement, h::UInt) = hash(c.n, h)
Base.isless(c1::ZNElement{N, p}, c2::ZNElement{N, p}) where {N, p} = isless(c1.n, c2.n)

# Experimental
BraidingStyle(::Type{ZNElement{2, 0}}) = Bosonic()
BraidingStyle(::Type{ZNElement{2, 1}}) = Fermionic()
BraidingStyle(::Type{ZNElement{N, 0}}) where {N} = Bosonic()
BraidingStyle(::Type{ZNElement{N, p}}) where {N, p} = Anyonic()
function Rsymbol(a::ZNElement{N, p}, b::ZNElement{N, p}, c::ZNElement{N, p}) where {N, p}
    if p == 0
        R = 1
    elseif N == 2 && p == 1
        R = ifelse(a.n == b.n == 1, -1.0, 1.0)
    else
        R = cispi(2 * p * a.n * b.n / N)
    end
    return ifelse(c == a * b, R, zero(R))
end
