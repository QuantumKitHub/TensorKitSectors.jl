# Deligne tensor product of different sectors: ⊠
#------------------------------------------------------------------------------#
const SectorTuple = Tuple{Sector, Vararg{Sector}}
const SectorNamedTuple = NamedTuple{<:Any, <:SectorTuple}
const AnySectorTuple = Union{SectorTuple, SectorNamedTuple}

"""
    ProductSector{T <: Union{SectorTuple, SectorNamedTuple}

Represents the Deligne tensor product of sectors. The type parameter `T` is a (named) tuple
of the component sectors. The recommended way to construct a `ProductSector` is using the
[`deligneproduct`](@ref) (`⊠`) operator on the components.

It is also possible to use names for the different factors in the product, either by passing
a `NamedTuple` of `Sector`s to the `ProductSector` constructor, or by making use of the
[`@NamedProductSector`](@ref) macro.
"""
struct ProductSector{T <: AnySectorTuple} <: Sector
    sectors::T
    ProductSector{T}(x::T) where {T} = new{T}(x)
end

ProductSector(x::AnySectorTuple) = ProductSector{typeof(x)}(x)
ProductSector(x, y, z...) = ProductSector((x, y, z...))
ProductSector(x::Sector) = ProductSector((x,))

ProductSector{T}(x) where {T} = ProductSector{T}(convert(T, x))
ProductSector{T}(x, y...) where {T} = ProductSector{T}((x, y...))
ProductSector{T}(x::Sector) where {T} = ProductSector{T}((x,))

# necessary in Julia 1.10:
ProductSector{NamedTuple{K, V}}(x::V) where {K, V} = ProductSector{NamedTuple{K, V}}(NamedTuple{K}(x))
ProductSector{NamedTuple{K, V}}(x::Tuple) where {K, V} = ProductSector{NamedTuple{K, V}}(NamedTuple{K}(V(x)))

const TupleProductSector{T <: SectorTuple} = ProductSector{T}
const NamedProductSector{T <: SectorTuple} = ProductSector{<:NamedTuple{<:Any, T}}

TupleProductSector(s::NamedProductSector) = ProductSector(values(Tuple(s)))

Base.Tuple(a::ProductSector) = Tuple(a.sectors)

Base.getindex(s::ProductSector, i::Int) = getindex(Tuple(s), i)
Base.length(s::ProductSector) = length(Tuple(s))
Base.iterate(s::ProductSector, args...) = iterate(Tuple(s), args...)
Base.indexed_iterate(s::ProductSector, args...) = Base.indexed_iterate(s.sectors, args...)

Base.@constprop :aggressive function Base.getproperty(s::ProductSector{<:NamedTuple}, f::Symbol)
    sectors = getfield(s, :sectors)
    if f === :sectors
        return sectors
    else
        return getproperty(sectors, f)
    end
end
Base.propertynames(s::NamedProductSector) = (:sectors, propertynames(fieldtype(typeof(s), :sectors))...)

_sectors(::Type{NamedTuple{K, V}}) where {K, V} = _sectors(V)
_sectors(::Type{Tuple{}}) = ()
Base.@pure function _sectors(::Type{T}) where {T <: SectorTuple}
    return (Base.tuple_type_head(T), _sectors(Base.tuple_type_tail(T))...)
end

function Base.IteratorSize(::Type{SectorValues{ProductSector{T}}}) where {T}
    return Base.IteratorSize(Base.Iterators.product(map(values, _sectors(T))...))
end
function Base.size(::SectorValues{ProductSector{T}}) where {T}
    return map(s -> length(values(s)), _sectors(T))
end
Base.length(P::SectorValues{<:ProductSector}) = *(size(P)...)

function _length(iter::SectorValues{I}) where {I <: Sector}
    return Base.IteratorSize(iter) === Base.IsInfinite() ? typemax(Int) : length(iter)
end
function _size(::SectorValues{ProductSector{T}}) where {T}
    return map(s -> _length(values(s)), _sectors(T))
end
function Base.getindex(P::SectorValues{ProductSector{T}}, i::Int) where {T}
    I = manhattan_to_multidimensional_index(i, _size(P))
    return ProductSector{T}(getindex.(values.(_sectors(T)), I))
end
function findindex(P::SectorValues{ProductSector{T}}, c::ProductSector{T}) where {T}
    return to_manhattan_index(findindex.(values.(_sectors(T)), Tuple(c.sectors)), _size(P))
end

function Base.iterate(P::SectorValues{<:ProductSector}, i = 1)
    Base.IteratorSize(P) != Base.IsInfinite() && i > length(P) && return nothing
    return getindex(P, i), i + 1
end

function Base.convert(::Type{ProductSector{T}}, t::Tuple) where {T <: SectorTuple}
    return ProductSector{T}(convert(T, t))
end
function Base.convert(P::Type{<:NamedProductSector{T}}, t::Tuple) where {T <: SectorTuple}
    return P(convert(T, t))
end

function unit(::Type{T}) where {T <: ProductSector}
    UnitStyle(T) === GenericUnit() && throw_genericunit_error(T)
    return only(allunits(T))
end
function allunits(::Type{ProductSector{T}}) where {T}
    iterators = map(allunits, _sectors(T))
    return SectorSet{ProductSector{T}}(Base.Iterators.product(iterators...))
end

dual(p::ProductSector) = ProductSector(map(dual, p.sectors))
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

_firstsector(x::ProductSector) = Base.first(Tuple(x))
_tailsector(x::ProductSector) = ProductSector(Base.tail(Tuple(x)))

function Fsymbol(a::P, b::P, c::P, d::P, e::P, f::P) where {P <: ProductSector}
    heads = map(_firstsector, (a, b, c, d, e, f))
    F₁ = Fsymbol(heads...)
    length(a) == 1 && return F₁
    tails = map(_tailsector, (a, b, c, d, e, f))
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

function Rsymbol(a::P, b::P, c::P) where {P <: ProductSector}
    heads = map(_firstsector, (a, b, c))
    R₁ = Rsymbol(heads...)
    length(a) == 1 && return R₁
    tails = map(_tailsector, (a, b, c))
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

function Bsymbol(a::P, b::P, c::P) where {P <: ProductSector}
    heads = map(_firstsector, (a, b, c))
    B₁ = Bsymbol(heads...)
    length(a) == 1 && return B₁
    tails = map(_tailsector, (a, b, c))
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

function Asymbol(a::P, b::P, c::P) where {P <: ProductSector}
    heads = map(_firstsector, (a, b, c))
    A₁ = Asymbol(heads...)
    length(a) == 1 && return A₁
    tails = map(_tailsector, (a, b, c))
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

frobenius_schur_phase(p::ProductSector) = prod(frobenius_schur_phase, Tuple(p))
frobenius_schur_indicator(p::ProductSector) = prod(frobenius_schur_indicator, Tuple(p))

function fusiontensor(a::P, b::P, c::P) where {P <: ProductSector}
    heads = map(_firstsector, (a, b, c))
    C₁ = fusiontensor(heads...)
    length(a) == 1 && return C₁
    tails = map(_tailsector, (a, b, c))
    C₂ = fusiontensor(tails...)
    return _kron(C₁, C₂)
end

function FusionStyle(::Type{<:ProductSector{T}}) where {T}
    return mapreduce(FusionStyle, &, _sectors(T))
end
function UnitStyle(::Type{<:ProductSector{T}}) where {T}
    return mapreduce(UnitStyle, &, _sectors(T))
end
function BraidingStyle(::Type{<:ProductSector{T}}) where {T}
    return mapreduce(BraidingStyle, &, _sectors(T))
end
function Base.isreal(::Type{<:ProductSector{T}}) where {T}
    return mapreduce(isreal, &, _sectors(T))
end

fermionparity(P::ProductSector) = mapreduce(fermionparity, xor, P.sectors)

dim(p::ProductSector) = prod(dim, p.sectors)

Base.isequal(p1::ProductSector, p2::ProductSector) = isequal(p1.sectors, p2.sectors)
Base.hash(p::ProductSector, h::UInt) = hash(p.sectors, h)
function Base.isless(p1::ProductSector{T}, p2::ProductSector{T}) where {T}
    I1 = findindex.(values.(_sectors(T)), Tuple(p1.sectors))
    I2 = findindex.(values.(_sectors(T)), Tuple(p2.sectors))
    d1 = sum(I1) - length(I1)
    d2 = sum(I2) - length(I2)
    d1 < d2 && return true
    d1 > d2 && return false
    return isless(I1, I2)
end

# Default construction from tensor product of sectors
#-----------------------------------------------------

"""
    deligneproduct(s₁::Sector, s₂::Sector)
    deligneproduct(S₁::Type{<:Sector}, S₂::Type{<:Sector})
    s₁ ⊠ s₂

Given two sectors `s₁` and `s₂`, which label an isomorphism class of simple objects in a
fusion category ``C₁`` and ``C₂``, `s1 ⊠ s2` (obtained as `\\boxtimes+TAB`) labels the
isomorphism class of simple objects in the Deligne tensor product category ``C₁ ⊠ C₂``.

The Deligne tensor product also works in the type domain and for spaces and tensors. For
group representations, we have `Irrep[G₁] ⊠ Irrep[G₂] == Irrep[G₁ × G₂]`.
"""
function ⊠ end
const deligneproduct = ⊠

⊠(s1, s2, s3, s4...) = ⊠(⊠(s1, s2), s3, s4...)

⊠(s1::Sector, s2::Sector) = ProductSector((s1, s2))
⊠(s1::Trivial, s2::Trivial) = s1
⊠(s1::Sector, s2::Trivial) = s1
⊠(s1::Trivial, s2::Sector) = s2
⊠(p1::ProductSector, s2::Trivial) = p1
⊠(p1::ProductSector, s2::Sector) = ProductSector(tuple(p1.sectors..., s2))
⊠(s1::Trivial, p2::ProductSector) = p2
⊠(s1::Sector, p2::ProductSector) = ProductSector(tuple(s1, p2.sectors...))
⊠(p1::ProductSector, p2::ProductSector) = ProductSector(tuple(p1.sectors..., p2.sectors...))
function ⊠(p1::NamedProductSector, p2::NamedProductSector)
    K = (keys(p1.sectors)..., keys(p2.sectors)...)
    V = (values(p1.sectors)..., values(p2.sectors)...)
    return ProductSector(NamedTuple{K}(V))
end

Base.@assume_effects :foldable function ⊠(::Type{I₁}, ::Type{I₂}) where {I₁ <: Sector, I₂ <: Sector}
    I₁ === Trivial && return I₂
    I₂ === Trivial && return I₁
    I₁ <: ProductSector || return ⊠(ProductSector{Tuple{I₁}}, I₂)
    I₂ <: ProductSector || return ⊠(I₁, ProductSector{Tuple{I₂}})
    return _deligneproduct_impl(I₁, I₂)
end
Base.@assume_effects :foldable function _deligneproduct_impl(
        ::Type{ProductSector{T₁}}, ::Type{ProductSector{T₂}}
    ) where {T₁, T₂}
    return ProductSector{Tuple{_sectors(T₁)..., _sectors(T₂)...}}
end
Base.@assume_effects :foldable function _deligneproduct_impl(
        ::Type{ProductSector{NamedTuple{K₁, V₁}}}, ::Type{ProductSector{NamedTuple{K₂, V₂}}}
    ) where {K₁, V₁, K₂, V₂}
    return ProductSector{NamedTuple{(K₁..., K₂...), Tuple{_sectors(V₁)..., _sectors(V₂)...}}}
end

function Base.show(io::IO, P::TupleProductSector)
    sectors = Tuple(P)
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
function Base.show(io::IO, P::NamedProductSector)
    if get(io, :typeinfo, nothing) !== typeof(P)
        print(io, type_repr(typeof(P)))
    end
    print(io, "(")
    first = true
    for sector in P.sectors
        first || print(io, ", ")
        ioc = IOContext(io, :typeinfo => typeof(sector))
        show(ioc, sector)
        first = false
    end
    print(io, ")")
    return nothing
end

function type_repr(P::Type{<:TupleProductSector})
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

"""
    @NamedProductSector{key1::Type1, key2::Type2, ...}
    @NamedProductSector begin key1::Type1; key2::Type2; ...; end

This macro gives a more convenient syntax for declaring `ProductSector` types which have
names associated to the components. This is also used when printing `ProductSector` types
to e.g. the REPL.
"""
macro NamedProductSector(ex)
    return esc(:(ProductSector{@NamedTuple($ex)}))
end
function type_repr(P::Type{<:NamedProductSector})
    names = P.parameters[1].parameters[1]
    sectors = P.parameters[1].parameters[2].parameters
    iob = IOBuffer()
    print(iob, "@NamedProductSector{")
    first = true
    for (name, sector) in zip(names, sectors)
        first || print(iob, ", ")
        print(iob, name)
        print(iob, "::")
        print(iob, type_repr(sector))
        first = false
    end
    print(iob, "}")
    return String(take!(iob))
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
    sectors = Tuple(P)
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
