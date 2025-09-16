# Rep Z2 ⊕ Rep Z2 ≅ Ising is worked out here
# 𝒞 = 𝒟 = RepZ2 ≅ {1, ψ}, while ℳ = Vec ≅ {σ}
# this is mainly meant for testing within TensorKit without relying on MultiTensorKit
#------------------------------------------------------------------------------#

# 𝒞   ℳ
# ℳᵒᵖ 𝒟
struct IsingBimodule <: Sector
    row::Int
    col::Int
    label::Int
    function IsingBimodule(row::Int, col::Int, label::Int)
        1 <= row <= 2 && 1 <= col <= 2 ||
            throw(DomainError(lazy"Invalid subcategory ($row, $col)"))
        0 <= label <= (row == col) ||
            throw(ArgumentError(lazy"Invalid label $label for IsingBimodule subcategory ($row, $col)"))
        return new(row, col, label)
    end
end

const all_isingbimod_objects = (
    IsingBimodule(1, 1, 0), IsingBimodule(1, 1, 1), IsingBimodule(2, 1, 0),
    IsingBimodule(1, 2, 0), IsingBimodule(2, 2, 0), IsingBimodule(2, 2, 1),
)

Base.IteratorSize(::Type{SectorValues{IsingBimodule}}) = Base.SizeUnknown()
Base.iterate(::SectorValues{IsingBimodule}, i = 1) = iterate(all_isingbimod_objects, i)
Base.length(::SectorValues{IsingBimodule}) = length(all_isingbimod_objects)

⊗(a::IsingBimodule, b::IsingBimodule) = IsingBimoduleIterator(a, b)

struct IsingBimoduleIterator
    a::IsingBimodule
    b::IsingBimodule
end

Base.IteratorSize(::Type{IsingBimoduleIterator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{IsingBimoduleIterator}) = Base.HasEltype()
Base.eltype(::Type{IsingBimoduleIterator}) = IsingBimodule

function Base.iterate(iter::IsingBimoduleIterator, state = 0)
    a, b = iter.a, iter.b
    a.col == b.row || return nothing

    _state = (a.row == b.col == a.col) ? mod(a.label + b.label, 2) : state
    return state < (1 + (a.row == b.col && a.row != a.col)) ?
        (IsingBimodule(a.row, b.col, _state), state + 1) : nothing
end

Base.convert(::Type{IsingBimodule}, labels::NTuple{3, Int}) = IsingBimodule(labels...)

function Base.convert(::Type{IsingAnyon}, a::IsingBimodule) # identify RepZ2 ⊕ RepZ2 ≅ Ising
    (a.row != a.col) && return IsingAnyon(:σ)
    return IsingAnyon(a.label == 0 ? :I : :ψ)
end

FusionStyle(::Type{IsingBimodule}) = SimpleFusion() # no multiplicities
BraidingStyle(::Type{IsingBimodule}) = NoBraiding() # because of module categories

function Nsymbol(a::IsingBimodule, b::IsingBimodule, c::IsingBimodule)
    if (a.row != c.row) || (a.col != b.row) || (b.col != c.col)
        throw(ArgumentError("invalid fusion channel"))
    end
    return Nsymbol(convert(IsingAnyon, a), convert(IsingAnyon, b), convert(IsingAnyon, c))
end

function Fsymbol(a::I, b::I, c::I, d::I, e::I, f::I) where {I <: IsingBimodule}
    Nsymbol(a, b, e) && Nsymbol(e, c, d) &&
        Nsymbol(b, c, f) && Nsymbol(a, f, d) || return 0.0
    return Fsymbol(
        convert(IsingAnyon, a), convert(IsingAnyon, b), convert(IsingAnyon, c),
        convert(IsingAnyon, d), convert(IsingAnyon, e), convert(IsingAnyon, f)
    )
end

# ℳ ↔ ℳop when conjugating elements within these
Base.conj(a::IsingBimodule) = IsingBimodule(a.col, a.row, a.label)

rightunit(a::IsingBimodule) = IsingBimodule(a.col, a.col, 0)
leftunit(a::IsingBimodule) = IsingBimodule(a.row, a.row, 0)

function unit(a::IsingBimodule)
    a.row == a.col ||
        throw(DomainError(a, "unit of module category ($(a.row), $(a.col)) doesn't exist"))
    return IsingBimodule(a.row, a.col, 0)
end

Base.isone(a::IsingBimodule) = leftunit(a) == a == rightunit(a)

function unit(::Type{IsingBimodule})
    throw(ArgumentError("one of Type IsingBimodule doesn't exist"))
end

allunits(::Type{IsingBimodule}) = (IsingBimodule(1, 1, 0), IsingBimodule(2, 2, 0))

function Base.isless(a::IsingBimodule, b::IsingBimodule)
    return isless((a.col, a.row, a.label), (b.col, b.row, b.label))
end

function Base.hash(a::IsingBimodule, h::UInt)
    return hash(a.label, hash(a.row, hash(a.col, h)))
end

function Base.show(io::IO, a::IsingBimodule)
    if get(io, :typeinfo, nothing) === IsingBimodule
        print(io, (a.row, a.col, a.label))
    else
        print(io, "IsingBimodule", (a.row, a.col, a.label))
    end
    return nothing
end
