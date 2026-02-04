# ZNIrrep: irreps of Z_N are labelled by integers mod N; stored as UInt or UInt8
const SMALL_ZN_CUTOFF = (typemax(UInt8) + 1) ÷ 2

"""
    struct ZNIrrep{N} <: AbstractIrrep{ℤ{N}}
    ZNIrrep{N}(n::Integer)
    Irrep[ℤ{N}](n::Integer)

Represents irreps of the group ``ℤ_N`` for some value of `N`.
For `N` equals `2`, `3` or `4`, `ℤ{N}` can be replaced by `ℤ₂`, `ℤ₃`, and `ℤ₄`.
An arbitrary `Integer` `n` can be provided to the constructor, but only the value `mod(n, N)` is relevant.

The type of the stored integer (`UInt8`) requires `N ≤ $SMALL_ZN_CUTOFF`.
Larger values of `N` should use the [`LargeZNIrrep`](@ref) instead.
The constructor `Irrep[ℤ{N}]` should be preferred, as it will automatically select the most efficient storage type for a given value of `N`.

See also [`charge`](@ref) and [`modulus`](@ref) to extract the relevant data.

## Fields
- `n::UInt8`: the integer label of the irrep, modulo `N`.
"""
struct ZNIrrep{N} <: AbstractIrrep{ℤ{N}}
    n::UInt8
    function ZNIrrep{N}(n::Integer) where {N}
        N ≤ SMALL_ZN_CUTOFF || throw(DomainError(N, "N exceeds the maximal value, use `LargeZNIrrep` instead"))
        return new{N}(UInt8(mod(n, N)))
    end
end

"""
    struct LargeZNIrrep{N} <: AbstractIrrep{ℤ{N}}
    LargeZNIrrep{N}(n::Integer)
    Irrep[ℤ{N}](n::Integer)

Represents irreps of the group ``ℤ_N`` for some value of `N`, which is typically larger than $SMALL_ZN_CUTOFF.
For smaller values of `N`, the [`ZNIrrep`](@ref) sector type should be used instead.
An arbitrary `Integer` `n` can be provided to the constructor, but only the value `mod(n, N)` is relevant.

The constructor `Irrep[ℤ{N}]` should be preferred, as it will automatically select the most efficient storage type for a given value of `N`.

See also [`charge`](@ref) and [`modulus`](@ref) to extract the relevant data.

## Fields
- `n::UInt`: the integer label of the irrep, modulo `N`.
"""
struct LargeZNIrrep{N} <: AbstractIrrep{ℤ{N}}
    n::UInt
    function LargeZNIrrep{N}(n::Integer) where {N}
        N ≤ (typemax(UInt) ÷ 2) || throw(DomainError(N, "N exceeds the maximal value"))
        return new{N}(UInt(mod(n, N)))
    end

end

const AnyZNIrrep{N} = Union{ZNIrrep{N}, LargeZNIrrep{N}}
const Z2Irrep = ZNIrrep{2}
const Z3Irrep = ZNIrrep{3}
const Z4Irrep = ZNIrrep{4}

"""
    modulus(c::ZNIrrep{N}) -> N
    modulus(::Type{<:ZNIrrep{N}}) -> N

The order of the cyclic group, or the modulus of the charge labels.
"""
modulus(c::AnyZNIrrep) = modulus(typeof(c))
modulus(::Type{<:AnyZNIrrep{N}}) where {N} = N

"""
    charge(c::ZNIrrep) -> Int

The charge label of the irrep `c`.
"""
charge(c::AnyZNIrrep) = Int(c.n)

Base.getindex(::IrrepTable, ::Type{ℤ{N}}) where {N} = N ≤ SMALL_ZN_CUTOFF ? ZNIrrep{N} : LargeZNIrrep{N}
Base.convert(Z::Type{<:AnyZNIrrep}, n::Real) = Z(n)

unit(::Type{ZNIrrep{N}}) where {N} = ZNIrrep{N}(zero(UInt8))
unit(::Type{LargeZNIrrep{N}}) where {N} = LargeZNIrrep{N}(zero(UInt))
# be careful with `-` for unsigned integers!
dual(c::AnyZNIrrep{N}) where {N} = typeof(c)(N - c.n)
⊗(c1::I, c2::I) where {I <: AnyZNIrrep} = (I(c1.n + c2.n),)

Base.IteratorSize(::Type{SectorValues{<:ZNIrrep}}) = HasLength()
# for larger values it doesn't make sense to store the sectors as a tuple
Base.IteratorSize(::Type{SectorValues{<:LargeZNIrrep}}) = SizeUnknown()

Base.length(::SectorValues{I}) where {I <: AnyZNIrrep} = modulus(I)
Base.iterate(::SectorValues{I}, i = 0) where {I <: AnyZNIrrep} = i == modulus(I) ? nothing : (I(i), i + 1)
function Base.getindex(::SectorValues{I}, i::Int) where {I <: AnyZNIrrep}
    return 1 <= i <= modulus(I) ? I(i - 1) : throw(BoundsError(values(I), i))
end
findindex(::SectorValues{I}, c::I) where {I <: AnyZNIrrep} = charge(c) + 1

Base.hash(c::AnyZNIrrep, h::UInt) = hash(c.n, h)
Base.isless(c1::I, c2::I) where {I <: AnyZNIrrep} = isless(c1.n, c2.n)

# ensure the printing uses `Int`.
function Base.show(io::IO, c::AnyZNIrrep)
    I = typeof(c)
    print_type = get(io, :typeinfo, nothing) !== I
    print_type && print(io, type_repr(I), '(')
    print(io, charge(c))
    print_type && print(io, ')')
    return nothing
end

# TimeReversed print
_tr_repr(a::AnyZNIrrep) = charge(a)
