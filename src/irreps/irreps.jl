# Sectors corresponding to irreducible representations of compact groups
#------------------------------------------------------------------------------#
"""
    abstract type AbstractIrrep{G <: Group} <: Sector

Abstract supertype for sectors which corresponds to irreps (irreducible representations) of
a group `G`. As we assume unitary representations, these would be finite groups or compact
Lie groups. Note that this could also include projective rather than linear representations.

Actual concrete implementations of those irreps can be obtained as `Irrep[G]`, or via their
actual name, which generically takes the form `(asciiG)Irrep`, i.e. the ASCII spelling of
the group name followed by `Irrep`.

All irreps have [`BraidingStyle`](@ref) equal to `Bosonic()` and thus trivial twists.
"""
abstract type AbstractIrrep{G <: Group} <: Sector end # irreps have integer quantum dimensions
BraidingStyle(::Type{<:AbstractIrrep}) = Bosonic()

struct IrrepTable end
"""
    const Irrep

A constant of a singleton type used as `Irrep[G]` with `G <: Group` a type of group, to
construct or obtain a concrete subtype of `AbstractIrrep{G}` that implements the data
structure used to represent irreducible representations of the group `G`.
"""
const Irrep = IrrepTable()

type_repr(::Type{IR}) where {G <: Group, IR <: AbstractIrrep{G}} = "Irrep[" * type_repr(G) * "]"
function Base.show(io::IO, c::AbstractIrrep)
    I = typeof(c)
    return if get(io, :typeinfo, nothing) !== I
        print(io, type_repr(I), "(")
        for k in 1:fieldcount(I)
            k > 1 && print(io, ", ")
            print(io, getfield(c, k))
        end
        print(io, ")")
    else
        fieldcount(I) > 1 && print(io, "(")
        for k in 1:fieldcount(I)
            k > 1 && print(io, ", ")
            print(io, getfield(c, k))
        end
        fieldcount(I) > 1 && print(io, ")")
    end
end

const AbelianIrrep{G} = AbstractIrrep{G} where {G <: AbelianGroup}
FusionStyle(::Type{<:AbelianIrrep}) = UniqueFusion()
Base.isreal(::Type{<:AbelianIrrep}) = true

Nsymbol(a::I, b::I, c::I) where {I <: AbelianIrrep} = c == first(a ⊗ b)
function Fsymbol(a::I, b::I, c::I, d::I, e::I, f::I) where {I <: AbelianIrrep}
    return Int(Nsymbol(a, b, e) * Nsymbol(e, c, d) * Nsymbol(b, c, f) * Nsymbol(a, f, d))
end
frobenius_schur_phase(a::AbelianIrrep) = 1
Asymbol(a::I, b::I, c::I) where {I <: AbelianIrrep} = Int(Nsymbol(a, b, c))
Bsymbol(a::I, b::I, c::I) where {I <: AbelianIrrep} = Int(Nsymbol(a, b, c))
Rsymbol(a::I, b::I, c::I) where {I <: AbelianIrrep} = Int(Nsymbol(a, b, c))

function fusiontensor(a::I, b::I, c::I) where {I <: AbelianIrrep}
    return fill(Int(Nsymbol(a, b, c)), (1, 1, 1, 1))
end

include("znirrep.jl")
include("u1irrep.jl")
include("dnirrep.jl")
include("cu1irrep.jl")
include("su2irrep.jl")
include("a4irrep.jl")
