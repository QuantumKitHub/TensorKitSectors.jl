# Time reversed sector

"""
    struct TimeReversed{I <: Sector}
    TimeReversed(a::I) where {I <: Sector}

Represents the time-reversed version of the sector `I`, i.e. the
sector with the same fusion rules and `F`-symbols, but with the
inverse braiding.
"""
struct TimeReversed{I <: Sector} <: Sector
    a::I
    function TimeReversed{I}(a::I) where {I <: Sector}
        if BraidingStyle(I) isa NoBraiding
            throw(ArgumentError("TimeReversed is not defined for sectors $I with no braiding"))
        end
        return new{I}(a)
    end
end

TimeReversed(a::I) where {I <: Sector} = TimeReversed{I}(a)
TimeReversed(a::TimeReversed{<:Sector}) = a.a

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
function fusiontensor(
        a::TimeReversed{I}, b::TimeReversed{I}, c::TimeReversed{I}
    ) where {I <: Sector}
    return fusiontensor(a.a, b.a, c.a) # conj here? all fusiontensors available are real
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
function Base.show(io::IO, a::TimeReversed{I}) where {I <: Sector}
    if get(io, :typeinfo, nothing) !== typeof(a)
        print(io, type_repr(typeof(a)))
    end
    # print(io, "(", _tr_repr(a.a), ")") # FIXME: this adds extra brackets to sectorset printing as well
    print(io, "(", a.a, ")")
    return nothing
end
