# SU2Irrep: irreps of SU2 are labelled by half integers j
struct SU2IrrepException <: Exception end
function Base.show(io::IO, ::SU2IrrepException)
    return print(
        io,
        "Irreps of (bosonic or fermionic) `SU₂` should be labelled by non-negative half integers, i.e. elements of `Rational{Int}` with denominator 1 or 2"
    )
end

"""
    struct SU2Irrep <: AbstractIrrep{SU₂}
    SU2Irrep(j::Real)
    Irrep[SU₂](j::Real)

Represents irreps of the group ``SU₂``. The irrep is labelled by a half integer `j` which
can be entered as an abitrary `Real`, but is stored as a `HalfInt` from the HalfIntegers.jl
package.

## Fields
- `j::HalfInt`: the label of the irrep, which can be any non-negative half integer.
"""
struct SU2Irrep <: AbstractIrrep{SU₂}
    j::HalfInt
    function SU2Irrep(j)
        j >= zero(j) || error("Not a valid SU₂ irrep")
        return new(j)
    end
end
Base.getindex(::IrrepTable, ::Type{SU₂}) = SU2Irrep
Base.convert(::Type{SU2Irrep}, j::Real) = SU2Irrep(j)

const _su2one = SU2Irrep(zero(HalfInt))
Base.one(::Type{SU2Irrep}) = _su2one
allones(::Type{SU2Irrep}) = (_su2one,)
Base.conj(s::SU2Irrep) = s
⊗(s1::SU2Irrep, s2::SU2Irrep) = SectorSet{SU2Irrep}(abs(s1.j - s2.j):(s1.j + s2.j))

Base.IteratorSize(::Type{SectorValues{SU2Irrep}}) = IsInfinite()
Base.iterate(::SectorValues{SU2Irrep}, i::Int = 0) = (SU2Irrep(half(i)), i + 1)
function Base.getindex(::SectorValues{SU2Irrep}, i::Int)
    return 1 <= i ? SU2Irrep(half(i - 1)) : throw(BoundsError(values(SU2Irrep), i))
end
findindex(::SectorValues{SU2Irrep}, s::SU2Irrep) = twice(s.j) + 1

dim(s::SU2Irrep) = twice(s.j) + 1

FusionStyle(::Type{SU2Irrep}) = SimpleFusion()
sectorscalartype(::Type{SU2Irrep}) = Float64
Base.isreal(::Type{SU2Irrep}) = true

Nsymbol(sa::SU2Irrep, sb::SU2Irrep, sc::SU2Irrep) = WignerSymbols.δ(sa.j, sb.j, sc.j)

function Fsymbol(
        s1::SU2Irrep, s2::SU2Irrep, s3::SU2Irrep,
        s4::SU2Irrep, s5::SU2Irrep, s6::SU2Irrep
    )
    if all(==(_su2one), (s1, s2, s3, s4, s5, s6))
        return 1.0
    else
        return sqrtdim(s5) * sqrtdim(s6) *
            WignerSymbols.racahW(
            sectorscalartype(SU2Irrep), s1.j, s2.j, s4.j, s3.j,
            s5.j, s6.j
        )
    end
end

function Rsymbol(sa::SU2Irrep, sb::SU2Irrep, sc::SU2Irrep)
    Nsymbol(sa, sb, sc) || return zero(sectorscalartype(SU2Irrep))
    return iseven(convert(Int, sa.j + sb.j - sc.j)) ? one(sectorscalartype(SU2Irrep)) :
        -one(sectorscalartype(SU2Irrep))
end

function fusiontensor(a::SU2Irrep, b::SU2Irrep, c::SU2Irrep)
    C = Array{Float64}(undef, dim(a), dim(b), dim(c), 1)
    ja, jb, jc = a.j, b.j, c.j

    for kc in 1:dim(c), kb in 1:dim(b), ka in 1:dim(a)
        C[ka, kb, kc, 1] = WignerSymbols.clebschgordan(
            ja, ja + 1 - ka, jb, jb + 1 - kb, jc,
            jc + 1 - kc
        )
    end
    return C
end

Base.hash(s::SU2Irrep, h::UInt) = hash(s.j, h)
Base.isless(s1::SU2Irrep, s2::SU2Irrep) = isless(s1.j, s2.j)
