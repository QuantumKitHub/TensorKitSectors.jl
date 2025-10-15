# U1Irrep: irreps of U1 are labelled by integers
"""
    struct U1Irrep <: AbstractIrrep{U₁}
    U1Irrep(charge::Real)
    Irrep[U₁](charge::Real)

Represents irreps of the group ``U₁``. The irrep is labelled by a charge, which should be
an integer for a linear representation. However, it is often useful to allow half integers
to represent irreps of ``U₁`` subgroups of ``SU₂``, such as the ``S^z`` of spin-1/2 system.
Hence, the charge is stored as a `HalfInt` from the package HalfIntegers.jl, but can be
entered as arbitrary `Real`. The sequence of the charges is: 0, 1/2, -1/2, 1, -1, ...

## Fields
- `charge::HalfInt`: the label of the irrep, which can be any half integer.
"""
struct U1Irrep <: AbstractIrrep{U₁}
    charge::HalfInt
end
Base.getindex(::IrrepTable, ::Type{U₁}) = U1Irrep
Base.convert(::Type{U1Irrep}, c::Real) = U1Irrep(c)

"""
    charge(c::U1Irrep) -> HalfInt

The charge label of the irrep `c`.
"""
charge(c::U1Irrep) = Int(c.charge)

unit(::Type{U1Irrep}) = U1Irrep(0)
dual(c::U1Irrep) = U1Irrep(-charge(c))
⊗(c1::U1Irrep, c2::U1Irrep) = (U1Irrep(charge(c1) + charge(c2)),)

Base.IteratorSize(::Type{SectorValues{U1Irrep}}) = IsInfinite()
function Base.iterate(::SectorValues{U1Irrep}, i::Int = 0)
    return i <= 0 ? (U1Irrep(half(i)), (-i + 1)) : (U1Irrep(half(i)), -i)
end
function Base.getindex(::SectorValues{U1Irrep}, i::Int)
    i < 1 && throw(BoundsError(values(U1Irrep), i))
    return U1Irrep(iseven(i) ? half(i >> 1) : -half(i >> 1))
end
function findindex(::SectorValues{U1Irrep}, c::U1Irrep)
    return (n = twice(charge(c)); 2 * abs(n) + (n <= 0))
end

Base.hash(c::U1Irrep, h::UInt) = hash(c.charge, h)
@inline function Base.isless(c1::U1Irrep, c2::U1Irrep)
    return isless(abs(charge(c1)), abs(charge(c2))) || zero(HalfInt) < charge(c1) == -charge(c2)
end
