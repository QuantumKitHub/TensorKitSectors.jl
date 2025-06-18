# Rep Z2 ⊕ Rep Z2 ≅ Ising is worked out here
# 𝒞 = 𝒟 = RepZ2 ≅ {1, ψ}, while ℳ = Vec ≅ {σ}
# this is mainly meant for testing within TensorKit without relying on MultiTensorKit

abstract type Bimodule <: Sector end
"""
    CatType
    
𝒞   ℳ
ℳᵒᵖ 𝒟
"""
@enum CatType 𝒞 = 1 ℳ = 3 ℳᵒᵖ = 2 𝒟 = 4

struct IsingBimod <: Bimodule
    type::CatType
    label::Int
    function IsingBimod(type::CatType, label::Int)
        if type == 𝒞
            0 ≤ label ≤ 1 ||
                throw(ArgumentError("Invalid 𝒞 label for Ising bimodule: $(label)"))
        elseif type == ℳ
            label == 0 ||
                throw(ArgumentError("Invalid ℳ label for Ising bimodule: $(label)"))
        elseif type == ℳᵒᵖ
            label == 0 ||
                throw(ArgumentError("Invalid ℳᵒᵖ label for Ising bimodule: $(label)"))
        elseif type == 𝒟
            0 ≤ label ≤ 1 ||
                throw(ArgumentError("Invalid 𝒟 label for Ising bimodule: $(label)"))
        end
        return new(type, label)
    end
end

isC(x::IsingBimod) = x.type == 𝒞
isM(x::IsingBimod) = x.type == ℳ
isMop(x::IsingBimod) = x.type == ℳᵒᵖ
isD(x::IsingBimod) = x.type == 𝒟

function isModule(a::IsingBimod)
    return a.type in (ℳ, ℳᵒᵖ)
end

⊗(a::IsingBimod, b::IsingBimod) = IsingBimodIterator(a, b)

struct IsingBimodIterator
    a::IsingBimod
    b::IsingBimod
end

Base.IteratorSize(::Type{IsingBimodIterator}) = Base.HasLength()
Base.IteratorEltype(::Type{IsingBimodIterator}) = Base.HasEltype()
Base.eltype(::Type{IsingBimodIterator}) = IsingBimod

function Base.length(iter::IsingBimodIterator)
    if (isC(iter.a) && isC(iter.b)) || (isD(iter.a) && isD(iter.b))
        return 1
    elseif (isM(iter.a) && isMop(iter.b)) || (isMop(iter.a) && isM(iter.b)) # σ x σ = 1 + ψ
        return 2
    elseif (isM(iter.a) && isD(iter.b)) || (isC(iter.a) && isM(iter.b)) ||
           (isD(iter.a) && isMop(iter.b)) || (isMop(iter.a) && isC(iter.b))
        return 1
    else
        return 0 # some a⊗b can't exist, e.g. 𝒞 ⊗ 𝒟
    end
end

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
    isModule(a) && return IsingAnyon(:σ)
    return IsingAnyon(a.label == 0 ? :I : :ψ)
end

function Nsymbol(a::IsingBimod, b::IsingBimod, c::IsingBimod)
    return c ∈ a ⊗ b ?
           Nsymbol(convert(IsingAnyon, a), convert(IsingAnyon, b), convert(IsingAnyon, c)) :
           false
end

function Fsymbol(a::I, b::I, c::I, d::I, e::I, f::I) where {I<:IsingBimod}
    return Fsymbol(convert(IsingAnyon, a), convert(IsingAnyon, b), convert(IsingAnyon, c),
                   convert(IsingAnyon, d), convert(IsingAnyon, e), convert(IsingAnyon, f))
end

function Rsymbol(a::IsingBimod, b::IsingBimod, c::IsingBimod)
    a.type == b.type == c.type ||
        throw(ArgumentError("can't braid between different categories"))
    ℳ ∉ map(i -> i.type, (a, b, c)) && ℳᵒᵖ ∉ map(i -> i.type, (a, b, c)) ||
        throw(ArgumentError("can't braid with module category"))
    return Rsymbol(convert(IsingAnyon, a), convert(IsingAnyon, b), convert(IsingAnyon, c))
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

FusionStyle(::Type{IsingBimod}) = SimpleFusion() # no multiplicities
BraidingStyle(::Type{IsingBimod}) = Anyonic() # see R symbols

Base.isreal(::Type{IsingBimod}) = false

function Base.show(io::IO, a::IsingBimod)
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

vertex_labeltype(::Type{IsingBimod}) = Nothing
VectorInterface.scalartype(::Type{IsingBimod}) = ComplexF64

dim(a::IsingBimod) = dim(convert(IsingAnyon, a))