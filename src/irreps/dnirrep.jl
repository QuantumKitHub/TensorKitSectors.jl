"""
    struct DNIrrep{N} <: AbstractIrrep{D{N}}
    DNIrrep{N}(n::Integer, isodd::Bool=false)
    Irrep[D{N}](n::Integer, isodd::Bool=false)

Represents irreps of the dihedral group ``D_N = Z_N ⋊ C`` (``Z_N`` and charge conjugation or reflection).

## Properties

- `j::Int`: the value of the ``Z_N`` charge.
- `isodd::Bool`: the representation of charge conjugation.

Combined these take the values ``+0, -0, 1, ..., (N - 1) / 2`` for odd ``N``, and
``+0, -0, 1, ..., N / 2 - 1, +(N/2), -(N/2)`` for even ``N``, where the ``+`` (``-``)
refer to the even (odd) one-dimensional irreps, while the others are two-dimensional.
"""
struct DNIrrep{N} <: AbstractIrrep{D{N}}
    # store isodd in right bit, use rest for storing the angle
    data::UInt8
    DNIrrep{N}(data::UInt8) where {N} = new{N}(data)
end

function DNIrrep{N}(j::Integer, isodd::Bool = false) where {N}
    @assert 0 < N < 32 "DNIrrep requires 0 < N < 32 to function properly"
    0 <= j <= (N >> 1) ||
        throw(DomainError(j, "DNIrrep only has irreps 0 <= j <= (N >> 1)"))
    !isodd || j == 0 || (iseven(N) && (j == (N >> 1))) ||
        throw(DomainError(j, "DNIrrep only has odd irreps when `j == 0` or `iseven(N) && j == N / 2`"))
    return DNIrrep{N}((j % UInt8) << 1 | isodd)
end
