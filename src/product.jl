# Deligne tensor product of different sectors: ⊠
#------------------------------------------------------------------------------#
const SectorTuple = Tuple{Vararg{Sector}}

"""
    ProductSector{T<:SectorTuple}

Represents the Deligne tensor product of sectors. The type parameter `T` is a tuple of the
component sectors. The recommended way to construct a `ProductSector` is using the
[`deligneproduct`](@ref) (`⊠`) operator on the components.
"""
struct ProductSector{T <: SectorTuple} <: Sector
    sectors::T
end

Base.getindex(s::ProductSector, i::Int) = getindex(s.sectors, i)
Base.iterate(s::ProductSector, args...) = iterate(s.sectors, args...)
Base.indexed_iterate(s::ProductSector, args...) = Base.indexed_iterate(s.sectors, args...)

_sectors(::Type{Tuple{}}) = ()
Base.@pure function _sectors(::Type{T}) where {T <: SectorTuple}
    return (Base.tuple_type_head(T), _sectors(Base.tuple_type_tail(T))...)
end

function Base.IteratorSize(::Type{SectorValues{ProductSector{T}}}) where {T <: SectorTuple}
    return Base.IteratorSize(Base.Iterators.product(map(values, _sectors(T))...))
end
function Base.size(::SectorValues{ProductSector{T}}) where {T <: SectorTuple}
    return map(s -> length(values(s)), _sectors(T))
end
Base.length(P::SectorValues{<:ProductSector}) = *(size(P)...)

function _length(iter::SectorValues{I}) where {I <: Sector}
    return Base.IteratorSize(iter) === Base.IsInfinite() ? typemax(Int) : length(iter)
end
function _size(::SectorValues{ProductSector{T}}) where {T <: SectorTuple}
    return map(s -> _length(values(s)), _sectors(T))
end
function Base.getindex(P::SectorValues{ProductSector{T}}, i::Int) where {T <: SectorTuple}
    I = manhattan_to_multidimensional_index(i, _size(P))
    return ProductSector{T}(getindex.(values.(_sectors(T)), I))
end
function findindex(
        P::SectorValues{ProductSector{T}},
        c::ProductSector{T}
    ) where {T <: SectorTuple}
    return to_manhattan_index(findindex.(values.(_sectors(T)), c.sectors), _size(P))
end

function Base.iterate(P::SectorValues{ProductSector{T}}, i = 1) where {T <: SectorTuple}
    Base.IteratorSize(P) != Base.IsInfinite() && i > length(P) && return nothing
    return getindex(P, i), i + 1
end

ProductSector{T}(args...) where {T <: SectorTuple} = ProductSector{T}(args)
function Base.convert(::Type{ProductSector{T}}, t::Tuple) where {T <: SectorTuple}
    return ProductSector{T}(convert(T, t))
end

Base.one(::Type{ProductSector{T}}) where {I <: Sector, T <: Tuple{I}} = ProductSector((one(I),))
function Base.one(::Type{ProductSector{T}}) where {I <: Sector, T <: Tuple{I, Vararg{Sector}}}
    return one(I) ⊠ one(ProductSector{Base.tuple_type_tail(T)})
end

Base.conj(p::ProductSector) = ProductSector(map(conj, p.sectors))
function ⊗(p1::P, p2::P) where {P <: ProductSector}
    if FusionStyle(P) isa UniqueFusion
        (P(first(product(map(⊗, p1.sectors, p2.sectors)...))),)
    else
        return SectorSet{P}(product(map(⊗, p1.sectors, p2.sectors)...))
    end
end

function Nsymbol(a::P, b::P, c::P) where {P <: ProductSector}
    return prod(map(Nsymbol, a.sectors, b.sectors, c.sectors))
end

_firstsector(x::ProductSector) = x.sectors[1]
_tailsector(x::ProductSector) = ProductSector(Base.tail(x.sectors))

function Fsymbol(a::P, b::P, c::P, d::P, e::P, f::P) where {P <: ProductSector}
    heads = map(_firstsector, (a, b, c, d, e, f))
    tails = map(_tailsector, (a, b, c, d, e, f))
    F₁ = Fsymbol(heads...)
    F₂ = Fsymbol(tails...)
    if F₁ isa Number && F₂ isa Number
        return F₁ * F₂
    elseif F₁ isa Number
        a₁, b₁, c₁, d₁, e₁, f₁ = heads
        sz₁ = (
            Nsymbol(a₁, b₁, e₁), Nsymbol(e₁, c₁, d₁), Nsymbol(b₁, c₁, f₁), Nsymbol(a₁, f₁, d₁),
        )
        F₁′ = fill(F₁, sz₁)
        return _kron(F₁′, F₂)
    elseif F₂ isa Number
        a₂, b₂, c₂, d₂, e₂, f₂ = tails
        sz₂ = (
            Nsymbol(a₂, b₂, e₂), Nsymbol(e₂, c₂, d₂), Nsymbol(b₂, c₂, f₂), Nsymbol(a₂, f₂, d₂),
        )
        F₂′ = fill(F₂, sz₂)
        return _kron(F₁, F₂′)
    else
        return _kron(F₁, F₂)
    end
end
function Fsymbol(
        a::P, b::P, c::P, d::P, e::P, f::P
    ) where {P <: ProductSector{<:Tuple{Sector}}}
    return Fsymbol(map(_firstsector, (a, b, c, d, e, f))...)
end

function Rsymbol(a::P, b::P, c::P) where {P <: ProductSector}
    heads = map(_firstsector, (a, b, c))
    tails = map(_tailsector, (a, b, c))
    R₁ = Rsymbol(heads...)
    R₂ = Rsymbol(tails...)
    if R₁ isa Number && R₂ isa Number
        R₁ * R₂
    elseif R₁ isa Number
        a₁, b₁, c₁ = heads
        sz₁ = (Nsymbol(a₁, b₁, c₁), Nsymbol(b₁, a₁, c₁)) # 0 x 0 or 1 x 1
        R₁′ = fill(R₁, sz₁)
        return _kron(R₁′, R₂)
    elseif R₂ isa Number
        a₂, b₂, c₂ = tails
        sz₂ = (Nsymbol(a₂, b₂, c₂), Nsymbol(b₂, a₂, c₂)) # 0 x 0 or 1 x 1
        R₂′ = fill(R₂, sz₂)
        return _kron(R₁, R₂′)
    else
        return _kron(R₁, R₂)
    end
end
function Rsymbol(a::P, b::P, c::P) where {P <: ProductSector{<:Tuple{Sector}}}
    return Rsymbol(map(_firstsector, (a, b, c))...)
end

function Bsymbol(a::P, b::P, c::P) where {P <: ProductSector}
    heads = map(_firstsector, (a, b, c))
    tails = map(_tailsector, (a, b, c))
    B₁ = Bsymbol(heads...)
    B₂ = Bsymbol(tails...)
    if B₁ isa Number && B₂ isa Number
        B₁ * B₂
    elseif B₁ isa Number
        a₁, b₁, c₁ = heads
        sz₁ = (Nsymbol(a₁, b₁, c₁), Nsymbol(c₁, dual(b₁), a₁)) # 0 x 0 or 1 x 1
        B₁′ = fill(B₁, sz₁)
        return _kron(B₁′, B₂)
    elseif B₂ isa Number
        a₂, b₂, c₂ = tails
        sz₂ = (Nsymbol(a₂, b₂, c₂), Nsymbol(c₂, dual(b₂), a₂)) # 0 x 0 or 1 x 1
        B₂′ = fill(B₂, sz₂)
        return _kron(B₁, B₂′)
    else
        return _kron(B₁, B₂)
    end
end
function Bsymbol(a::P, b::P, c::P) where {P <: ProductSector{<:Tuple{Sector}}}
    return Bsymbol(map(_firstsector, (a, b, c))...)
end

function Asymbol(a::P, b::P, c::P) where {P <: ProductSector}
    heads = map(_firstsector, (a, b, c))
    tails = map(_tailsector, (a, b, c))
    A₁ = Asymbol(heads...)
    A₂ = Asymbol(tails...)
    if A₁ isa Number && A₂ isa Number
        A₁ * A₂
    elseif A₁ isa Number
        a₁, b₁, c₁ = heads
        sz₁ = (Nsymbol(a₁, b₁, c₁), Nsymbol(dual(a₁), c₁, b₁)) # 0 x 0 or 1 x 1
        A₁′ = fill(A₁, sz₁)
        return _kron(A₁′, A₂)
    elseif A₂ isa Number
        a₂, b₂, c₂ = tails
        sz₂ = (Nsymbol(a₂, b₂, c₂), Nsymbol(dual(a₂), c₂, b₂)) # 0 x 0 or 1 x 1
        A₂′ = fill(A₂, sz₂)
        return _kron(A₁, A₂′)
    else
        return _kron(A₁, A₂)
    end
end
function Asymbol(a::P, b::P, c::P) where {P <: ProductSector{<:Tuple{Sector}}}
    return Asymbol(map(_firstsector, (a, b, c))...)
end

frobeniusschur(p::ProductSector) = prod(map(frobeniusschur, p.sectors))

function fusiontensor(a::P, b::P, c::P) where {P <: ProductSector}
    return _kron(
        fusiontensor(map(_firstsector, (a, b, c))...),
        fusiontensor(map(_tailsector, (a, b, c))...)
    )
end

function fusiontensor(a::P, b::P, c::P) where {P <: ProductSector{<:Tuple{Sector}}}
    return fusiontensor(map(_firstsector, (a, b, c))...)
end

function FusionStyle(::Type{<:ProductSector{T}}) where {T <: SectorTuple}
    return Base.:&(map(FusionStyle, _sectors(T))...)
end
function BraidingStyle(::Type{<:ProductSector{T}}) where {T <: SectorTuple}
    return Base.:&(map(BraidingStyle, _sectors(T))...)
end
Base.isreal(::Type{<:ProductSector{T}}) where {T <: SectorTuple} = _isreal(T)
_isreal(::Type{Tuple{}}) = true
function _isreal(T::Type{<:SectorTuple})
    return isreal(Base.tuple_type_head(T)) && _isreal(Base.tuple_type_tail(T))
end

fermionparity(P::ProductSector) = mapreduce(fermionparity, xor, P.sectors)

dim(p::ProductSector) = *(dim.(p.sectors)...)

Base.isequal(p1::ProductSector, p2::ProductSector) = isequal(p1.sectors, p2.sectors)
Base.hash(p::ProductSector, h::UInt) = hash(p.sectors, h)
function Base.isless(p1::ProductSector{T}, p2::ProductSector{T}) where {T <: SectorTuple}
    I1 = findindex.(values.(_sectors(T)), p1.sectors)
    I2 = findindex.(values.(_sectors(T)), p2.sectors)
    d1 = sum(I1) - length(I1)
    d2 = sum(I2) - length(I2)
    d1 < d2 && return true
    d1 > d2 && return false
    return isless(I1, I2)
end

# Default construction from tensor product of sectors
#-----------------------------------------------------
⊠(s1, s2, s3, s4...) = ⊠(⊠(s1, s2), s3, s4...)
const deligneproduct = ⊠

"""
    ⊠(s₁::Sector, s₂::Sector)
    deligneproduct(s₁::Sector, s₂::Sector)

Given two sectors `s₁` and `s₂`, which label an isomorphism class of simple objects in a
fusion category ``C₁`` and ``C₂``, `s1 ⊠ s2` (obtained as `\\boxtimes+TAB`) labels the
isomorphism class of simple objects in the Deligne tensor product category ``C₁ ⊠ C₂``.

The Deligne tensor product also works in the type domain and for spaces and tensors. For
group representations, we have `Irrep[G₁] ⊠ Irrep[G₂] == Irrep[G₁ × G₂]`.
"""
⊠(s1::Sector, s2::Sector) = ProductSector((s1, s2))
⊠(s1::Trivial, s2::Trivial) = s1
⊠(s1::Sector, s2::Trivial) = s1
⊠(s1::Trivial, s2::Sector) = s2
⊠(p1::ProductSector, s2::Trivial) = p1
⊠(p1::ProductSector, s2::Sector) = ProductSector(tuple(p1.sectors..., s2))
⊠(s1::Trivial, p2::ProductSector) = p2
⊠(s1::Sector, p2::ProductSector) = ProductSector(tuple(s1, p2.sectors...))
⊠(p1::ProductSector, p2::ProductSector) = ProductSector(tuple(p1.sectors..., p2.sectors...))

⊠(I1::Type{Trivial}, I2::Type{Trivial}) = Trivial
⊠(I1::Type{Trivial}, I2::Type{<:ProductSector}) = I2
⊠(I1::Type{Trivial}, I2::Type{<:Sector}) = I2

⊠(I1::Type{<:ProductSector}, I2::Type{Trivial}) = I1
@static if VERSION >= v"1.8"
    Base.@assume_effects :foldable function ⊠(
            I1::Type{<:ProductSector}, I2::Type{<:ProductSector}
        )
        T1 = I1.parameters[1]
        T2 = I2.parameters[1]
        return ProductSector{Tuple{T1.parameters..., T2.parameters...}}
    end
else
    Base.@pure function ⊠(I1::Type{<:ProductSector}, I2::Type{<:ProductSector})
        T1 = I1.parameters[1]
        T2 = I2.parameters[1]
        return ProductSector{Tuple{T1.parameters..., T2.parameters...}}
    end
end
⊠(I1::Type{<:ProductSector}, I2::Type{<:Sector}) = I1 ⊠ ProductSector{Tuple{I2}}

⊠(I1::Type{<:Sector}, I2::Type{Trivial}) = I1
⊠(I1::Type{<:Sector}, I2::Type{<:ProductSector}) = ProductSector{Tuple{I1}} ⊠ I2
⊠(I1::Type{<:Sector}, I2::Type{<:Sector}) = ProductSector{Tuple{I1, I2}}

function Base.show(io::IO, P::ProductSector)
    sectors = P.sectors
    compact = get(io, :typeinfo, nothing) === typeof(P)
    sep = compact ? ", " : " ⊠ "
    print(io, "(")
    for i in 1:length(sectors)
        i == 1 || print(io, sep)
        io2 = compact ? IOContext(io, :typeinfo => typeof(sectors[i])) : io
        print(io2, sectors[i])
    end
    return print(io, ")")
end

function type_repr(P::Type{<:ProductSector})
    sectors = P.parameters[1].parameters
    if length(sectors) == 1
        s = "ProductSector{Tuple{" * type_repr(sectors[1]) * "}}"
    else
        s = "("
        for i in 1:length(sectors)
            if i != 1
                s *= " ⊠ "
            end
            s *= type_repr(sectors[i])
        end
        s *= ")"
    end
    return s
end

#==============================================================================
TODO: the following would implement pretty-printing of product sectors, i.e.
`ProductSector{Tuple{Irrep[G]}}` would be printed as `Irrep[G]`, and
`ProductSector{Tuple{Irrep[G]}}(x)` would be printed as `Irrep[G](x)`.
However, defining show for a type is considered type piracy/treason, and can lead
to unexpected behavior. While we can avoid this by only defining show for the
instances, this would lead to the following behavior:

```julia-repl
julia> [Irrep[ℤ₂ × U₁](0, 0)]
1-element Vector{TensorKit.ProductSector{Tuple{Z2Irrep, U1Irrep}}}:
 (0, 0)
```

See Julia issues #29988, #29428, #22363, #28983.

Base.show(io::IO, P::Type{<:ProductSector}) = print(io, type_repr(P))
==============================================================================#
function Base.show(io::IO, P::ProductSector{T}) where {T <: Tuple{Vararg{AbstractIrrep}}}
    sectors = P.sectors
    get(io, :typeinfo, nothing) === typeof(P) || print(io, type_repr(typeof(P)))
    print(io, "(")
    for i in 1:length(sectors)
        i == 1 || print(io, ", ")
        print(IOContext(io, :typeinfo => typeof(sectors[i])), sectors[i])
    end
    return print(io, ")")
end

function type_repr(::Type{ProductSector{T}}) where {T <: Tuple{Vararg{AbstractIrrep}}}
    sectors = T.parameters
    s = "Irrep["
    for i in 1:length(sectors)
        if i != 1
            s *= " × "
        end
        s *= type_repr(supertype(sectors[i]).parameters[1])
    end
    s *= "]"
    return s
end

function Base.getindex(::IrrepTable, ::Type{ProductGroup{Gs}}) where {Gs <: GroupTuple}
    G1 = tuple_type_head(Gs)
    Grem = tuple_type_tail(Gs)
    return ProductSector{Tuple{Irrep[G1]}} ⊠ Irrep[ProductGroup{tuple_type_tail(Gs)}]
end
function Base.getindex(::IrrepTable, ::Type{ProductGroup{Tuple{G}}}) where {G <: Group}
    return ProductSector{Tuple{Irrep[G]}}
end
