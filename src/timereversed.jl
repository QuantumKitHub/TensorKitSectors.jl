# Time reversed sector

"""
    struct TimeReversed{I <: Sector}
    TimeReversed(a::Sector)
    timereversed(a::Sector)

Represents the time-reversed version of the sector `a`, i.e. the
sector with the same fusion rules and `F`-symbols, but with the
inverse braiding.
Time reversal acts trivially on sectors with symmetric braiding.
In such cases, `timereversed(a)` simply returns `a`.
"""
struct TimeReversed{I <: Sector} <: Sector
    a::I
    function TimeReversed{I}(a) where {I <: Sector}
        BraidingStyle(I) isa NoBraiding &&
            throw(ArgumentError("TimeReversed is not defined for sectors $I with no braiding"))
        return new{I}(a)
    end
end

TimeReversed(a::Sector) = TimeReversed{typeof(a)}(a)

"""
    timereversed(a::I) where {I <: Sector}

Return the time-reversed version of the sector `a`. 
If `a` is already a time-reversed sector, return the original sector.
"""
function timereversed(a::Sector)
    return BraidingStyle(typeof(a)) isa SymmetricBraiding ? a : TimeReversed(a)
end
timereversed(a::TimeReversed) = a.a

FusionStyle(::Type{TimeReversed{I}}) where {I <: Sector} = FusionStyle(I)
BraidingStyle(::Type{TimeReversed{I}}) where {I <: Sector} = BraidingStyle(I)
function Nsymbol(
        a::TimeReversed{I}, b::TimeReversed{I}, c::TimeReversed{I}
    ) where {I <: Sector}
    return Nsymbol(a.a, b.a, c.a)
end

function Fsymbol(
        a::TimeReversed{I}, b::TimeReversed{I}, c::TimeReversed{I},
        d::TimeReversed{I}, e::TimeReversed{I}, f::TimeReversed{I}
    ) where {I <: Sector}
    return Fsymbol(a.a, b.a, c.a, d.a, e.a, f.a)
end
function Rsymbol(
        a::TimeReversed{I}, b::TimeReversed{I}, c::TimeReversed{I}
    ) where {I <: Sector}
    return adjoint(Rsymbol(a.a, b.a, c.a))
end
dim(c::TimeReversed{I}) where {I <: Sector} = dim(c.a)

unit(::Type{TimeReversed{I}}) where {I <: Sector} = TimeReversed{I}(unit(I))
allunits(::Type{TimeReversed{I}}) where {I <: Sector} = SectorSet{TimeReversed{I}}(TimeReversed{I}, allunits(I))
dual(c::TimeReversed{I}) where {I <: Sector} = TimeReversed{I}(dual(c.a))
function ⊗(c1::TimeReversed{I}, c2::TimeReversed{I}) where {I <: Sector}
    return SectorSet{TimeReversed{I}}(TimeReversed{I}, c1.a ⊗ c2.a)
end
function Base.IteratorSize(::Type{SectorValues{TimeReversed{I}}}) where {I <: Sector}
    return Base.IteratorSize(values(I))
end
function Base.length(::SectorValues{TimeReversed{I}}) where {I <: Sector}
    return length(values(I))
end
function Base.getindex(::SectorValues{TimeReversed{I}}, i::Int) where {I <: Sector}
    return TimeReversed{I}(getindex(values(I), i))
end
function Base.iterate(::SectorValues{TimeReversed{I}}, state...) where {I <: Sector}
    next = iterate(values(I), state...)
    if isnothing(next)
        return nothing
    else
        obj, nextstate = next
        return TimeReversed{I}(obj), nextstate
    end
end
function findindex(
        ::SectorValues{TimeReversed{I}}, a::TimeReversed{I}
    ) where {I <: Sector}
    return findindex(values(I), a.a)
end

function Base.isless(c1::TimeReversed{I}, c2::TimeReversed{I}) where {I <: Sector}
    return isless(c1.a, c2.a)
end

fusionscalartype(::Type{TimeReversed{I}}) where {I} = fusionscalartype(I)
braidingscalartype(::Type{TimeReversed{I}}) where {I} = braidingscalartype(I)
sectorscalartype(::Type{TimeReversed{I}}) where {I} = sectorscalartype(I)
dimscalartype(::Type{TimeReversed{I}}) where {I} = dimscalartype(I)

type_repr(::Type{TimeReversed{I}}) where {I <: Sector} = "TimeReversed{$(type_repr(I))}"
function Base.show(io::IO, c::TimeReversed{I}) where {I <: Sector}
    print_info = get(io, :typeinfo, nothing) !== typeof(c)
    print_info && print(io, type_repr(typeof(c)), "(")
    show(IOContext(io, :typeinfo => I), c.a)
    print_info && print(io, ")")
    return nothing
end
