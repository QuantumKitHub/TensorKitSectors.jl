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
        row <= 2 && col <= 2 || throw(DomainError(lazy"Invalid subcategory ($row, $col)"))
        if label < 0 || label > 1 || (label == 1 && (row != col))
            throw(ArgumentError(lazy"Invalid label $label for IsingBimod subcategory ($row, $col)"))
        end
        return new(row, col, label)
    end
end

isC(x::IsingBimod) = (x.row == x.col == 1)
isM(x::IsingBimod) = (x.row == 1 && x.col == 2)
isMop(x::IsingBimod) = (x.row == 2 && x.col == 1)
isD(x::IsingBimod) = (x.row == 2 && x.col == 2)

function ismodulecategory(a::IsingBimod)
    return (isM(a) || isMop(a))
end

const all_isingbimod_objects = (IsingBimod(1, 1, 0), IsingBimod(1, 1, 1),
                                IsingBimod(2, 1, 0),
                                IsingBimod(1, 2, 0), IsingBimod(2, 2, 0),
                                IsingBimod(2, 2, 1))

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
    if isC(iter.a) && isC(iter.b) # 𝒞 × 𝒞 -> 𝒞
        return state == 0 ? (IsingBimod(1, 1, mod(iter.a.label + iter.b.label, 2)), 1) :
               nothing
    end

    if isD(iter.a) && isD(iter.b) # 𝒟 × 𝒟 -> 𝒟
        return state == 0 ? (IsingBimod(2, 2, mod(iter.a.label + iter.b.label, 2)), 1) :
               nothing
    end

    if isM(iter.a) && isMop(iter.b) # ℳ × ℳop -> 𝒞
        return state < 2 ? (IsingBimod(1, 1, state), state + 1) : nothing
    end

    if isMop(iter.a) && isM(iter.b) # ℳop × ℳ -> 𝒟
        return state < 2 ? (IsingBimod(2, 2, state), state + 1) : nothing
    end

    if isC(iter.a) && isM(iter.b) # 𝒞 × ℳ -> ℳ
        return state == 0 ? (iter.b, 1) : nothing
    end

    if isMop(iter.a) && isC(iter.b) # ℳop x 𝒞 -> ℳop
        return state == 0 ? (iter.a, 1) : nothing
    end

    if isM(iter.a) && isD(iter.b) # ℳ x 𝒟 -> ℳ
        return state == 0 ? (iter.a, 1) : nothing
    end

    if isD(iter.a) && isMop(iter.b) # 𝒟 x ℳop -> ℳop
        return state == 0 ? (iter.b, 1) : nothing
    end

    return nothing
end

function Base.convert(::Type{IsingAnyon}, a::IsingBimod) # identify RepZ2 ⊕ RepZ2 ≅ Ising
    ismodulecategory(a) && return IsingAnyon(:σ)
    return IsingAnyon(a.label == 0 ? :I : :ψ)
end

FusionStyle(::Type{IsingBimod}) = SimpleFusion() # no multiplicities
BraidingStyle(::Type{IsingBimod}) = NoBraiding() # because of module categories

function Nsymbol(a::IsingBimod, b::IsingBimod, c::IsingBimod)
    # if a and b can fuse, then so can dual(a) and c, and c and dual(b)
    # only needs to be explicitly checked when CatTypes differ or when there's a module category involved
    if (a.row, a.col) != (b.row, b.col) || (a.row, a.col) != (c.row, c.col) ||
       (b.row, b.col) != (c.row, c.col) ||
       any(ismodulecategory, (a, b, c))
        c ∈ a ⊗ b && dual(b) ∈ dual(c) ⊗ a && dual(a) ∈ b ⊗ dual(c) ||
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

function Base.conj(a::IsingBimod) # ℳ ↔ ℳop when conjugating elements within these
    if isC(a) || isD(a) # self-conjugate within RepZ2
        return a
    elseif isM(a)
        return IsingBimod(2, 1, a.label)
    else
        return IsingBimod(1, 2, a.label)
    end
end

function rightone(a::IsingBimod)
    if isC(a) || isD(a)
        return IsingBimod(a.row, a.col, 0)
    elseif isM(a) # ℳ as right-𝒟 module category
        return IsingBimod(2, 2, 0)
    else
        return IsingBimod(1, 1, 0)
    end
end

function leftone(a::IsingBimod)
    if isC(a) || isD(a)
        return IsingBimod(a.row, a.col, 0)
    elseif isM(a) # ℳ as left-𝒞 module category
        return IsingBimod(1, 1, 0)
    else
        return IsingBimod(2, 2, 0)
    end
end

function Base.one(a::IsingBimod)
    if isC(a) || isD(a)
        return IsingBimod(a.row, a.col, 0)
    else
        throw(DomainError("unit of module category ($(a.row), $(a.col)) doesn't exist"))
    end
end

Base.one(::Type{IsingBimod}) = throw(ArgumentError("one of Type IsingBimod doesn't exist"))

function Base.show(io::IO, ::MIME"text/plain", a::IsingBimod)
    if isC(a)
        print(io, "𝒞[$(a.label)]")
    elseif isM(a)
        print(io, "ℳ[$(a.label)]")
    elseif isMop(a)
        print(io, "ℳᵒᵖ[$(a.label)]")
    elseif isD(a)
        print(io, "𝒟[$(a.label)]")
    end
end

function Base.isless(a::IsingBimod, b::IsingBimod)
    vals = SectorValues{IsingBimod}()
    return isless(findindex(vals, a), findindex(vals, b)) # order 𝒞 < ℳᵒᵖ < ℳ < 𝒟, and then within each cat according to label
end

function Base.hash(a::IsingBimod, h::UInt)
    return hash(a.label, hash(a.row, hash(a.col, h)))
end
