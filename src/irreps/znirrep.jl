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
struct ZNIrrep{N} <: AbstractIrrep{ℤ{N}}
    n::Int8
    function ZNIrrep{N}(n::Integer) where {N}
        @assert N < 64
        return new{N}(mod(n, N))
    end
end
Base.getindex(::IrrepTable, ::Type{ℤ{N}}) where {N} = ZNIrrep{N}
Base.convert(Z::Type{<:ZNIrrep}, n::Real) = Z(n)
const Z2Irrep = ZNIrrep{2}
const Z3Irrep = ZNIrrep{3}
const Z4Irrep = ZNIrrep{4}

Base.one(::Type{ZNIrrep{N}}) where {N} = ZNIrrep{N}(0)
Base.conj(c::ZNIrrep{N}) where {N} = ZNIrrep{N}(-c.n)
⊗(c1::ZNIrrep{N}, c2::ZNIrrep{N}) where {N} = (ZNIrrep{N}(c1.n + c2.n),)

Base.IteratorSize(::Type{SectorValues{ZNIrrep{N}}}) where {N} = HasLength()
Base.length(::SectorValues{ZNIrrep{N}}) where {N} = N
function Base.iterate(::SectorValues{ZNIrrep{N}}, i = 0) where {N}
    return i == N ? nothing : (ZNIrrep{N}(i), i + 1)
end
function Base.getindex(::SectorValues{ZNIrrep{N}}, i::Int) where {N}
    return 1 <= i <= N ? ZNIrrep{N}(i - 1) : throw(BoundsError(values(ZNIrrep{N}), i))
end
findindex(::SectorValues{ZNIrrep{N}}, c::ZNIrrep{N}) where {N} = c.n + 1

Base.hash(c::ZNIrrep{N}, h::UInt) where {N} = hash(c.n, h)
Base.isless(c1::ZNIrrep{N}, c2::ZNIrrep{N}) where {N} = isless(c1.n, c2.n)
