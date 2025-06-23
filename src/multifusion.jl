# Rep Z2 ⊕ Rep Z2 ≅ Ising is worked out here
# 𝒞 = 𝒟 = RepZ2 ≅ {1, ψ}, while ℳ = Vec ≅ {σ}
# this is mainly meant for testing within TensorKit without relying on MultiTensorKit
#------------------------------------------------------------------------------#

@enum CatType 𝒞 = 1 ℳ = 3 ℳᵒᵖ = 2 𝒟 = 4
# possible TODO: get rid of CatType and use Int instead -> prevents the need to export 𝒞 etc
function Base.getindex(a::CatType, label::Int)
    return IsingBimod(a, label)
end

struct IsingBimod <: Sector
    type::CatType
    label::Int
    function IsingBimod(type::CatType, label::Int)
        if label < 0 || label > 1 || (label == 1 && (type == ℳ || type == ℳᵒᵖ))
            throw(ArgumentError(lazy"Invalid $type label for Ising bimodule: $(label)"))
        end
        return new(type, label)
    end
end

isC(x::IsingBimod) = x.type == 𝒞
isM(x::IsingBimod) = x.type == ℳ
isMop(x::IsingBimod) = x.type == ℳᵒᵖ
isD(x::IsingBimod) = x.type == 𝒟

function ismodulecategory(a::IsingBimod)
    return a.type in (ℳ, ℳᵒᵖ)
end

const all_isingbimod_objects = (IsingBimod(𝒞, 0), IsingBimod(𝒞, 1), IsingBimod(ℳᵒᵖ, 0),
                                IsingBimod(ℳ, 0), IsingBimod(𝒟, 0), IsingBimod(𝒟, 1))

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
        return state == 0 ? (IsingBimod(𝒞, mod(iter.a.label + iter.b.label, 2)), 1) :
               nothing
    end

    if isD(iter.a) && isD(iter.b) # 𝒟 × 𝒟 -> 𝒟
        return state == 0 ? (IsingBimod(𝒟, mod(iter.a.label + iter.b.label, 2)), 1) :
               nothing
    end

    if isM(iter.a) && isMop(iter.b) # ℳ × ℳop -> 𝒞
        return state < 2 ? (IsingBimod(𝒞, state), state + 1) : nothing
    end

    if isMop(iter.a) && isM(iter.b) # ℳop × ℳ -> 𝒟
        return state < 2 ? (IsingBimod(𝒟, state), state + 1) : nothing
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
    if a.type != b.type || a.type != c.type || b.type != c.type ||
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
        return IsingBimod(ℳᵒᵖ, a.label)
    else
        return IsingBimod(ℳ, a.label)
    end
end

function rightone(a::IsingBimod)
    if isC(a) || isD(a)
        return IsingBimod(a.type, 0)
    elseif isM(a) # ℳ as right-𝒟 module category
        return IsingBimod(𝒟, 0)
    else
        return IsingBimod(𝒞, 0)
    end
end

function leftone(a::IsingBimod)
    if isC(a) || isD(a)
        return IsingBimod(a.type, 0)
    elseif isM(a) # ℳ as left-𝒞 module category
        return IsingBimod(𝒞, 0)
    else
        return IsingBimod(𝒟, 0)
    end
end

function Base.one(a::IsingBimod)
    if isC(a) || isD(a)
        return IsingBimod(a.type, 0)
    else
        throw(DomainError("unit of module category $(a.type) doesn't exist"))
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
    return isless((a.type, a.label), (b.type, b.label)) # order 𝒞 < ℳᵒᵖ < ℳ < 𝒟, and then within each cat according to label
end

function Base.hash(a::IsingBimod, h::UInt)
    return hash(a.label, hash(a.type, h))
end

dim(a::IsingBimod) = dim(convert(IsingAnyon, a))
