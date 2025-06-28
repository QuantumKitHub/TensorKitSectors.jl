# Rep Z2 ⊕ Rep Z2 ≅ Ising is worked out here
# 𝒞 = 𝒟 = RepZ2 ≅ {1, ψ}, while ℳ = Vec ≅ {σ}
# this is mainly meant for testing within TensorKit without relying on MultiTensorKit
#------------------------------------------------------------------------------#

# 𝒞   ℳ 
# ℳᵒᵖ 𝒟
struct IsingBimod <: Sector
    row::Int
    col::Int
    label::Int
    function IsingBimod(row::Int, col::Int, label::Int)
        1 <= row <= 2 && 1 <= col <= 2 ||
            throw(DomainError(lazy"Invalid subcategory ($row, $col)"))
        if label < 0 || label > (row == col)
            throw(ArgumentError(lazy"Invalid label $label for IsingBimod subcategory ($row, $col)"))
        end
        return new(row, col, label)
    end
end

const all_isingbimod_objects = (IsingBimod(1, 1, 0), IsingBimod(1, 1, 1),
                                IsingBimod(2, 1, 0), IsingBimod(1, 2, 0),
                                IsingBimod(2, 2, 0), IsingBimod(2, 2, 1))

Base.IteratorSize(::Type{<:SectorValues{IsingBimod}}) = Base.SizeUnknown()
Base.iterate(::SectorValues{IsingBimod}, i=1) = iterate(all_isingbimod_objects, i)
Base.length(::SectorValues{IsingBimod}) = length(all_isingbimod_objects)

⊗(a::IsingBimod, b::IsingBimod) = IsingBimodIterator(a, b)

struct IsingBimodIterator
    a::IsingBimod
    b::IsingBimod
end

Base.IteratorSize(::Type{IsingBimodIterator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{IsingBimodIterator}) = Base.HasEltype()
Base.eltype(::Type{IsingBimodIterator}) = IsingBimod

function Base.iterate(iter::IsingBimodIterator, state=0)
    a, b = iter.a, iter.b
    a.col == b.row || return nothing

    _state = (a.row == b.col == a.col) ? mod(a.label + b.label, 2) : state
    return state < (1 + (a.row == b.col && a.row != a.col)) ?
           (IsingBimod(a.row, b.col, _state), state + 1) : nothing
end

function Base.convert(::Type{IsingAnyon}, a::IsingBimod) # identify RepZ2 ⊕ RepZ2 ≅ Ising
    (a.row != a.col) && return IsingAnyon(:σ)
    return IsingAnyon(a.label == 0 ? :I : :ψ)
end

FusionStyle(::Type{IsingBimod}) = SimpleFusion() # no multiplicities
BraidingStyle(::Type{IsingBimod}) = NoBraiding() # because of module categories

function Nsymbol(a::IsingBimod, b::IsingBimod, c::IsingBimod)
    if (a.row != c.row) || (a.col != b.row) || (b.col != c.col)
        throw(ArgumentError("invalid fusion channel"))
    end
    return Nsymbol(convert(IsingAnyon, a), convert(IsingAnyon, b), convert(IsingAnyon, c))
end

function Fsymbol(a::I, b::I, c::I, d::I, e::I, f::I) where {I<:IsingBimod}
    Nsymbol(a, b, e) && Nsymbol(e, c, d) &&
    Nsymbol(b, c, f) && Nsymbol(a, f, d) || return 0.0
    return Fsymbol(convert(IsingAnyon, a), convert(IsingAnyon, b), convert(IsingAnyon, c),
                   convert(IsingAnyon, d), convert(IsingAnyon, e), convert(IsingAnyon, f))
end

# ℳ ↔ ℳop when conjugating elements within these
Base.conj(a::IsingBimod) = IsingBimod(a.col, a.row, a.label)

rightone(a::IsingBimod) = IsingBimod(a.col, a.col, 0)
leftone(a::IsingBimod) = IsingBimod(a.row, a.row, 0)

function Base.one(a::IsingBimod)
    a.row == a.col ||
        throw(DomainError("unit of module category ($(a.row), $(a.col)) doesn't exist"))
    return IsingBimod(a.row, a.col, 0)
end

Base.one(::Type{IsingBimod}) = throw(ArgumentError("one of Type IsingBimod doesn't exist"))

function Base.isless(a::IsingBimod, b::IsingBimod)
    return isless((a.col, a.row, a.label), (b.col, b.row, b.label))
end

function Base.hash(a::IsingBimod, h::UInt)
    return hash(a.label, hash(a.row, hash(a.col, h)))
end
