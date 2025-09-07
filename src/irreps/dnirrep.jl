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

function Base.getproperty(a::DNIrrep{N}, sym::Symbol) where {N}
    if sym === :j
        return getfield(a, :data) >> 1
    elseif sym === :isodd
        return Bool(getfield(a, :data) & true)
    elseif sym === :data
        return getfield(a, :data)
    else
        error("Unknown property $sym")
    end
end

Base.propertynames(x::DNIrrep) = (:j, :isodd, :data)

Base.convert(::Type{DNIrrep{N}}, (j, n)::Tuple{Integer, Bool}) where {N} = DNIrrep{N}(j, n)

function Base.show(io::IO, a::DNIrrep)
    if get(io, :typeinfo, nothing) !== typeof(a)
        print(io, type_repr(typeof(a)))
    end
    print(io, "(", a.j, ", ", a.isodd, ")")
    return nothing
end
function Base.show(io::IO, ::MIME"text/plain", a::DNIrrep{N}) where {N}
    print_type = get(io, :typeinfo, nothing) !== typeof(a)
    print_type && print(io, type_repr(typeof(a)), "(")
    print(io, '"')
    if a.isodd
        print(io, '-')
    elseif a.j == 0 || (iseven(N) && a.j == N >> 1)
        print(io, '+')
    end
    print(io, a.j, '"')
    print_type && print(io, ")")
    return nothing
end

