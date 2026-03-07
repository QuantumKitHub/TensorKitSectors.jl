# Named sector: similar to ProductSector but with named tuple storage
#------------------------------------------------------------------------------#
const NamedSectorTuple = NamedTuple{<:Any, <:SectorTuple}

"""
    struct NamedSector{NT <: NamedSectorTuple}
    NamedSector((; name₁=s₁, name₂=s₂, ...))

Represents the Deligne tensor product of sectors stored with named components.
Similar to [`ProductSector`](@ref), but components are accessible by name via
`s.name`.

# Examples
```julia
julia> s = NamedSector((; charge=U1Irrep(1), spin=SU2Irrep(1//2)))
(charge=1, spin=1/2)

julia> s.charge
1

julia> s.spin
1/2
```
"""
struct NamedSector{NT <: NamedSectorTuple} <: Sector
    sectors::NT
end

# Construction
NamedSector(; kwargs...) = NamedSector(values(kwargs))
NamedSector{T}(args...) where {T <: NamedSectorTuple} = NamedSector{T}(args)
NamedSector{NT}(args::Vararg{Sector}) where {NT <: NamedSectorTuple} = NamedSector{NT}(NT(args))
function Base.convert(::Type{NamedSector{NT}}, nt::NamedTuple) where {NT <: NamedSectorTuple}
    return NamedSector{NT}(convert(NT, nt))
end

Base.NamedTuple(a::NamedSector) = a.sectors
Base.Tuple(a::NamedSector) = Tuple(a.sectors)

function Base.getproperty(s::NamedSector, name::Symbol)
    name === :sectors && return getfield(s, :sectors)
    return getfield(s, :sectors)[name]
end
Base.propertynames(s::NamedSector) = keys(s.sectors)

Base.getindex(s::NamedSector, i::Int) = s.sectors[i]
Base.getindex(s::NamedSector, name::Symbol) = s.sectors[name]
Base.length(s::NamedSector) = length(s.sectors)
Base.iterate(s::NamedSector, args...) = iterate(s.sectors, args...)
Base.indexed_iterate(s::NamedSector, args...) = Base.indexed_iterate(s.sectors, args...)
Base.keys(s::NamedSector) = keys(s.sectors)

# Type-level helpers
_sectornames(::Type{NamedSector{NT}}) where {NT} = _sectornames(NT)
_sectornames(::Type{NamedTuple{names, T}}) where {names, T} = names
_sectortupletype(::Type{NamedSector{NT}}) where {NT} = _sectortupletype(NT)
_sectortupletype(::Type{NamedTuple{names, T}}) where {names, T} = T

function _to_productsector(p::NamedSector{NT}) where {NT <: NamedSectorTuple}
    return ProductSector{_sectortupletype(NT)}(Tuple(p.sectors))
end

# SectorValues iteration
function Base.IteratorSize(::Type{SectorValues{NamedSector{NT}}}) where {NT <: NamedSectorTuple}
    T = _sectortupletype(NT)
    return Base.IteratorSize(Base.Iterators.product(map(values, _sectors(T))...))
end
function Base.size(::SectorValues{NamedSector{NT}}) where {NT <: NamedSectorTuple}
    T = _sectortupletype(NT)
    return map(s -> length(values(s)), _sectors(T))
end
Base.length(P::SectorValues{<:NamedSector}) = *(size(P)...)

function _size(::SectorValues{NamedSector{NT}}) where {NT <: NamedSectorTuple}
    T = _sectortupletype(NT)
    return map(s -> _length(values(s)), _sectors(T))
end
function Base.getindex(P::SectorValues{NamedSector{NT}}, i::Int) where {NT <: NamedSectorTuple}
    names = _sectornames(NT)
    T = _sectortupletype(NT)
    I = manhattan_to_multidimensional_index(i, _size(P))
    sectors_i = getindex.(values.(_sectors(T)), I)
    return NamedSector(NamedTuple{names}(sectors_i))
end
function findindex(
        P::SectorValues{NamedSector{NT}},
        c::NamedSector{NT}
    ) where {NT <: NamedSectorTuple}
    T = _sectortupletype(NT)
    return to_manhattan_index(findindex.(values.(_sectors(T)), Tuple(c.sectors)), _size(P))
end
function Base.iterate(P::SectorValues{NamedSector{NT}}, i = 1) where {NT <: NamedSectorTuple}
    Base.IteratorSize(P) != Base.IsInfinite() && i > length(P) && return nothing
    return getindex(P, i), i + 1
end

# Unit
function unit(::Type{T}) where {T <: NamedSector}
    UnitStyle(T) === GenericUnit() && throw_genericunit_error(T)
    return only(allunits(T))
end
function allunits(::Type{NamedSector{NT}}) where {NT <: NamedSectorTuple}
    names = _sectornames(NT)
    T = _sectortupletype(NT)
    iterators = map(allunits, _sectors(T))
    f = t -> NamedSector(NamedTuple{names}(t))
    return SectorSet{NamedSector{NT}}(f, Base.Iterators.product(iterators...))
end
leftunit(a::P) where {P <: NamedSector} = P(map(leftunit, a.sectors))
rightunit(a::P) where {P <: NamedSector} = P(map(rightunit, a.sectors))

# Sector operations — delegate to ProductSector where recursion is needed
dual(p::NamedSector) = NamedSector(map(dual, p.sectors))

function ⊗(p1::P, p2::P) where {P <: NamedSector}
    t1, t2 = Tuple(p1.sectors), Tuple(p2.sectors)
    names = _sectornames(P)
    if FusionStyle(P) isa UniqueFusion
        t = first(product(map(⊗, t1, t2)...))
        return (P(NamedTuple{names}(t)),)
    else
        return SectorSet{P}(t -> P(NamedTuple{names}(t)), product(map(⊗, t1, t2)...))
    end
end

function Nsymbol(a::P, b::P, c::P) where {P <: NamedSector}
    return prod(map(Nsymbol, Tuple(a.sectors), Tuple(b.sectors), Tuple(c.sectors)))
end

function Fsymbol(a::P, b::P, c::P, d::P, e::P, f::P) where {P <: NamedSector}
    return Fsymbol(map(_to_productsector, (a, b, c, d, e, f))...)
end
function Rsymbol(a::P, b::P, c::P) where {P <: NamedSector}
    return Rsymbol(map(_to_productsector, (a, b, c))...)
end
function Bsymbol(a::P, b::P, c::P) where {P <: NamedSector}
    return Bsymbol(map(_to_productsector, (a, b, c))...)
end
function Asymbol(a::P, b::P, c::P) where {P <: NamedSector}
    return Asymbol(map(_to_productsector, (a, b, c))...)
end
function fusiontensor(a::P, b::P, c::P) where {P <: NamedSector}
    return fusiontensor(map(_to_productsector, (a, b, c))...)
end

# Style traits
function FusionStyle(::Type{NamedSector{NT}}) where {NT <: NamedSectorTuple}
    return mapreduce(FusionStyle, &, _sectors(_sectortupletype(NT)))
end
function fusionscalartype(::Type{NamedSector{NT}}) where {NT <: NamedSectorTuple}
    return typeof(prod(zero ∘ fusionscalartype, _sectors(_sectortupletype(NT))))
end
function UnitStyle(::Type{NamedSector{NT}}) where {NT <: NamedSectorTuple}
    return mapreduce(UnitStyle, &, _sectors(_sectortupletype(NT)))
end
function BraidingStyle(::Type{NamedSector{NT}}) where {NT <: NamedSectorTuple}
    return mapreduce(BraidingStyle, &, _sectors(_sectortupletype(NT)))
end
function braidingscalartype(::Type{NamedSector{NT}}) where {NT <: NamedSectorTuple}
    return typeof(prod(zero ∘ braidingscalartype, _sectors(_sectortupletype(NT))))
end
function sectorscalartype(::Type{NamedSector{NT}}) where {NT <: NamedSectorTuple}
    return if BraidingStyle(NamedSector{NT}) == NoBraiding()
        typeof(prod(zero ∘ fusionscalartype, _sectors(_sectortupletype(NT))))
    else
        typeof(prod(zero ∘ sectorscalartype, _sectors(_sectortupletype(NT))))
    end
end
function dimscalartype(::Type{NamedSector{NT}}) where {NT <: NamedSectorTuple}
    return typeof(prod(zero ∘ dimscalartype, _sectors(_sectortupletype(NT))))
end

fermionparity(p::NamedSector) = mapreduce(fermionparity, xor, Tuple(p.sectors))
frobenius_schur_phase(p::NamedSector) = prod(frobenius_schur_phase, Tuple(p.sectors))
frobenius_schur_indicator(p::NamedSector) = prod(frobenius_schur_indicator, Tuple(p.sectors))

dim(p::NamedSector) = *(dim.(Tuple(p.sectors))...)

# Equality, hashing, ordering
Base.isequal(p1::NamedSector, p2::NamedSector) = isequal(p1.sectors, p2.sectors)
Base.hash(p::NamedSector, h::UInt) = hash(p.sectors, h)
function Base.isless(p1::NamedSector{NT}, p2::NamedSector{NT}) where {NT <: NamedSectorTuple}
    T = _sectortupletype(NT)
    I1 = findindex.(values.(_sectors(T)), Tuple(p1.sectors))
    I2 = findindex.(values.(_sectors(T)), Tuple(p2.sectors))
    d1 = sum(I1) - length(I1)
    d2 = sum(I2) - length(I2)
    d1 < d2 && return true
    d1 > d2 && return false
    return isless(I1, I2)
end

# Display
function Base.show(io::IO, P::NamedSector)
    sectors = P.sectors
    if get(io, :typeinfo, nothing) === typeof(P)
        # compact mode: just show values
        print(io, "(")
        for (i, name) in enumerate(keys(sectors))
            i == 1 || print(io, ", ")
            print(io, name, "=")
            print(IOContext(io, :typeinfo => typeof(sectors[name])), sectors[name])
        end
        print(io, ")")
    else
        # full mode: emit valid Julia that reconstructs the instance
        print(io, "NamedSector(;")
        for (i, name) in enumerate(keys(sectors))
            i == 1 || print(io, ",")
            print(io, " ", name, "=")
            print(io, sectors[name])
        end
        print(io, ")")
    end
    return nothing
end

function type_repr(::Type{NamedSector{NT}}) where {NT <: NamedSectorTuple}
    names = _sectornames(NT)
    sectors = _sectors(_sectortupletype(NT))
    s = "@NamedSector{"
    for (i, (name, I)) in enumerate(zip(names, sectors))
        i == 1 || (s *= ", ")
        s *= string(name) * "::" * type_repr(I)
    end
    s *= "}"
    return s
end

"""
    @NamedSector{name₁::T₁, name₂::T₂, ...}

Convenience macro for constructing a `NamedSector` type with named components,
analogous to [`@NamedTuple`](@ref).

# Examples
```julia
julia> @NamedSector{charge::U1Irrep, spin::SU2Irrep}
NamedSector{@NamedTuple{charge::U1Irrep, spin::SU2Irrep}}

julia> (@NamedSector{charge::U1Irrep, spin::SU2Irrep})(U1Irrep(1), SU2Irrep(1//2))
(charge=1, spin=1/2)
```
"""
macro NamedSector(ex)
    Meta.isexpr(ex, :braces) || Meta.isexpr(ex, :block) ||
        throw(ArgumentError("@NamedSector expects {...} or begin...end"))
    decls = filter(e -> !(e isa LineNumberNode), ex.args)
    all(e -> e isa Symbol || Meta.isexpr(e, :(::)), decls) ||
        throw(ArgumentError("@NamedSector must contain a sequence of name or name::T expressions"))
    vars = [QuoteNode(e isa Symbol ? e : e.args[1]) for e in decls]
    types = [esc(e isa Symbol ? :Any : e.args[2]) for e in decls]
    return :(NamedSector{NamedTuple{($(vars...),), Tuple{$(types...)}}})
end
