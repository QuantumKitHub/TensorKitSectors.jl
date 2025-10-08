# ZNIrrep: irreps of Z_N are labelled by integers mod N; do we ever want N > 64?
"""
    struct ZNIrrep{N, T <: Unsigned} <: AbstractIrrep{ℤ{N}}
    ZNIrrep{N}(n::Integer)
    Irrep[ℤ{N}](n::Integer)

Represents irreps of the group ``ℤ_N`` for some value of `N`.
For `N` equals `2`, `3` or `4`, `ℤ{N}` can be replaced by [`ℤ₂`](@ref), [`ℤ₃`](@ref), and [`ℤ₄`](@ref).
An arbitrary `Integer` `n` can be provided to the constructor, but only the value `mod(n, N)` is relevant.

The type of the stored integer `T` can either be explicitly provided, or will automatically be determined
to be the smallest unsigned integer type that fits all possible irreps for the given `N`.

See also [`charge`](@ref)` and [`modulus`](@ref) to extract the relevant data.

## Fields
- `n::T`: the integer label of the irrep, modulo `N`.
"""
struct ZNIrrep{N, T <: Unsigned} <: AbstractIrrep{ℤ{N}}
    n::T
    function ZNIrrep{N}(n::Integer) where {N}
        T = _integer_type(N)
        return new{N, T}(T(mod(n, N)))
    end
    function ZNIrrep{N, T}(n::Integer) where {N, T <: Unsigned}
        N ≤ typemax(T) + 1 ||
            throw(TypeError(:ZNIrrep, ZNIrrep{N, T}, ZNIrrep{N, _integer_type(N)}))
        return new{N, T}(mod(n, N))
    end
end

"""
    modulus(c::ZNIrrep{N}) -> N
    modulus(::Type{<:ZNIrrep{N}}) -> N

The order of the cyclic group, or the modulus of the charge labels.
"""
modulus(c::ZNIrrep) = modulus(typeof(c))
modulus(::Type{<:ZNIrrep{N}}) where {N} = N

"""
    charge(c::ZNIrrep) -> Int

The charge label of the irrep `c`.
"""
charge(c::ZNIrrep) = Int(c.n)

Base.@assume_effects :foldable function _integer_type(N::Integer)
    N <= 0 && throw(DomainError(N, "N should be positive"))
    for T in (UInt8, UInt16, UInt32, UInt64)
        # T needs to fit a
        N ≤ (typemax(T) + 1) && return T
    end
    throw(DomainError(N, "N is too large"))
end

Base.getindex(::IrrepTable, ::Type{ℤ{N}}) where {N} = ZNIrrep{N, _integer_type(N)}
Base.convert(Z::Type{<:ZNIrrep}, n::Real) = Z(n)
const Z2Irrep = ZNIrrep{2, UInt8}
const Z3Irrep = ZNIrrep{3, UInt8}
const Z4Irrep = ZNIrrep{4, UInt8}

unit(::Type{ZNIrrep{N, T}}) where {N, T} = ZNIrrep{N, T}(zero(T))
# be careful with `-` for unsigned integers!
dual(c::ZNIrrep{N, T}) where {N, T} = ZNIrrep{N, T}(N - c.n)
⊗(c1::ZNIrrep{N, T}, c2::ZNIrrep{N, T}) where {N, T} = (ZNIrrep{N, T}(modular_add(c1.n, c2.n, Val(N))),)

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

# compute x + y mod N, requires 0 <= x < N and 0 <= y < N (unchecked!)
function modular_add(x::T, y::T, ::Val{N}) where {T <: Unsigned, N}
    Tmax = typemax(T) + 1
    0 < N ≤ Tmax || throw(DomainError(N, "N is too large"))
    N ≤ (typemax(T) + 1) ÷ 2 && return x + y
    r, flag = Base.add_with_overflow(x, y)
    return ifelse(flag, (Tmax - N) + r, r)
end
