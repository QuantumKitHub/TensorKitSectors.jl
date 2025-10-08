# ZNIrrep: irreps of Z_N are labelled by integers mod N; do we ever want N > 64?
"""
    struct ZNIrrep{N} <: AbstractIrrep{ℤ{N}}
    ZNIrrep{N}(n::Integer)
    Irrep[ℤ{N}](n::Integer)

Represents irreps of the group ``ℤ_N`` for some value of `N<64`. (We need 2*(N-1) <= 127 in
order for a ⊗ b to work correctly.) For `N` equals `2`, `3` or `4`, `ℤ{N}` can be replaced
by `ℤ₂`, `ℤ₃`, `ℤ₄`. An arbitrary `Integer` `n` can be provided to the constructor, but only
the value `mod(n, N)` is relevant.

## Fields
- `n::Int8`: the integer label of the irrep, modulo `N`.
"""
struct ZNIrrep{N, T <: Unsigned} <: AbstractIrrep{ℤ{N}}
    n::T
    function ZNIrrep{N}(n::Integer) where {N}
        T = _integer_type(N)
        return new{N, T}(T(mod(n, N)))
    end
    function ZNIrrep{N, T}(n::Integer) where {N, T <: Unsigned}
        N ≤ typemax(T) ||
            throw(TypeError(:ZNIrrep, ZNIrrep{N, T}, ZNIrrep{N, _integer_type(N)}))
        return new{N, T}(mod(n, N))
    end
end

modulus(c::ZNIrrep) = modulus(typeof(c))
modulus(::Type{<:ZNIrrep{N}}) where {N} = N

charge(c::ZNIrrep) = Int(c.N)

Base.@assume_effects :foldable function _integer_type(N::Integer)
    N <= 0 && throw(DomainError(N, "N should be positive"))
    for T in (UInt8, UInt16, UInt32, UInt64)
        N ≤ typemax(T) && return T
    end
    throw(DomainError(N, "N is too large"))
end

Base.getindex(::IrrepTable, ::Type{ℤ{N}}) where {N} = ZNIrrep{N, _integer_type(N)}
Base.convert(Z::Type{<:ZNIrrep}, n::Real) = Z(n)
const Z2Irrep = ZNIrrep{2, UInt8}
const Z3Irrep = ZNIrrep{3, UInt8}
const Z4Irrep = ZNIrrep{4, UInt8}

unit(::Type{ZNIrrep{N, T}}) where {N, T} = ZNIrrep{N, T}(zero(T))
dual(c::ZNIrrep{N, T}) where {N, T} = ZNIrrep{N, T}(N - c.n)
⊗(c1::ZNIrrep{N, T}, c2::ZNIrrep{N, T}) where {N, T} = (ZNIrrep{N, T}(c1.n + c2.n),)

Base.IteratorSize(::Type{SectorValues{ZNIrrep{N, T}}}) where {N, T} = HasLength()
Base.length(::SectorValues{ZNIrrep{N, T}}) where {N, T} = N
function Base.iterate(::SectorValues{ZNIrrep{N, T}}, i = 0) where {N, T}
    return i == N ? nothing : (ZNIrrep{N, T}(i), i + 1)
end
function Base.getindex(::SectorValues{ZNIrrep{N, T}}, i::Int) where {N, T}
    return 1 <= i <= N ? ZNIrrep{N, T}(i - 1) : throw(BoundsError(values(ZNIrrep{N}), i))
end
findindex(::SectorValues{ZNIrrep{N, T}}, c::ZNIrrep{N, T}) where {N, T} = c.n + 1

Base.hash(c::ZNIrrep{N, T}, h::UInt) where {N, T} = hash(c.n, h)
Base.isless(c1::ZNIrrep{N, T}, c2::ZNIrrep{N, T}) where {N, T} = isless(c1.n, c2.n)
