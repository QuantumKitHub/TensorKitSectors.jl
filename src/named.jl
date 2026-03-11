# Named sector: similar to ProductSector but with named tuple storage
#------------------------------------------------------------------------------#
const NamedSectorTuple = NamedTuple{<:Any, <:SectorTuple}

"""
    struct NamedSector{NT <: NamedSectorTuple}
    NamedSector((; name₁=s₁, name₂=s₂, ...))

Represents the Deligne tensor product of sectors stored with named components.
Similar to [`ProductSector`](@ref), but components are accessible by name via `s.name`.

# Examples

```julia
julia> s = NamedSector(; charge=U1Irrep(1), spin=SU2Irrep(1//2))
(charge=1, spin=1/2)

julia> s.charge
1

julia> s.spin
1/2
```

!!! warning
    The sectors are stored internally in the `s.sectors` field, so the named components should not contain the `sectors` as a name.
"""
struct NamedSector{NT <: NamedSectorTuple} <: Sector
    sectors::NT
end

# Construction
NamedSector(; kwargs...) = NamedSector(values(kwargs))
NamedSector{NT}(arg1, args...) where {NT} = NamedSector{NT}((arg1, args...))
NamedSector{NT}(args::Tuple) where {NT <: NamedSectorTuple} =
    NamedSector{NT}(NT(convert(_sectortupletype(NamedSector{NT}), args)))

function Base.convert(::Type{NamedSector{NT}}, nt::Union{Tuple, NamedTuple}) where {NT}
    return NamedSector{NT}(convert(NT, nt))
end

Base.NamedTuple(a::NamedSector) = a.sectors
Base.Tuple(a::NamedSector) = Tuple(a.sectors)

ProductSector(a::NamedSector) = ProductSector(Tuple(a))
NamedSector{NT}(a::ProductSector) where {NT <: NamedSectorTuple} = NamedSector{NT}(Tuple(a))

function Base.getproperty(s::NamedSector, name::Symbol)
    name === :sectors && return getfield(s, :sectors)
    return getfield(s, :sectors)[name]
end
Base.propertynames(s::NamedSector) = (:sectors, keys(s.sectors)...)

Base.getindex(s::NamedSector, i::Int) = s.sectors[i]
Base.getindex(s::NamedSector, name::Symbol) = s.sectors[name]
Base.length(s::NamedSector) = length(s.sectors)
Base.iterate(s::NamedSector, args...) = iterate(s.sectors, args...)
Base.indexed_iterate(s::NamedSector, args...) = Base.indexed_iterate(s.sectors, args...)
Base.keys(s::NamedSector) = keys(s.sectors)

# Type-level helpers
_sectornames(::Type{NamedSector{NamedTuple{Names, T}}}) where {Names, T} = Names
_sectortupletype(::Type{NamedSector{NamedTuple{Names, T}}}) where {Names, T} = T
_sectortupletype(::Type) = error("should never be reached") # keeps JET happy
_sectors(::Type{I}) where {I <: NamedSector} = fieldtypes(_sectortupletype(I))
_productsectortype(::Type{I}) where {I <: NamedSector} = ProductSector{_sectortupletype(I)}

# SectorValues iteration
Base.IteratorSize(::Type{SectorValues{I}}) where {I <: NamedSector} =
    Base.IteratorSize(values(_productsectortype(I)))
Base.size(::SectorValues{I}) where {I <: NamedSector} =
    size(values(_productsectortype(I)))
Base.length(::SectorValues{I}) where {I <: NamedSector} =
    length(values(_productsectortype(I)))

Base.getindex(::SectorValues{I}, i::Int) where {I <: NamedSector} =
    I(getindex(values(_productsectortype(I)), i))
function findindex(P::SectorValues{I}, c::I) where {I <: NamedSector}
    return findindex(values(_productsectortype(I)), ProductSector(c))
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
allunits(::Type{I}) where {I <: NamedSector} = SectorSet{I}(I, allunits(_productsectortype(I)))
leftunit(a::P) where {P <: NamedSector} = P(leftunit(ProductSector(a)))
rightunit(a::P) where {P <: NamedSector} = P(rightunit(ProductSector(a)))

# Sector operations
dual(p::P) where {P <: NamedSector} = P(dual(ProductSector(p)))

⊗(p1::P, p2::P) where {P <: NamedSector} = SectorSet{P}(P, ProductSector(p1) ⊗ ProductSector(p2))

function Nsymbol(a::P, b::P, c::P) where {P <: NamedSector}
    return Nsymbol(map(ProductSector, (a, b, c))...)
end
function Fsymbol(a::P, b::P, c::P, d::P, e::P, f::P) where {P <: NamedSector}
    return Fsymbol(map(ProductSector, (a, b, c, d, e, f))...)
end
function Rsymbol(a::P, b::P, c::P) where {P <: NamedSector}
    return Rsymbol(map(ProductSector, (a, b, c))...)
end
function Bsymbol(a::P, b::P, c::P) where {P <: NamedSector}
    return Bsymbol(map(ProductSector, (a, b, c))...)
end
function Asymbol(a::P, b::P, c::P) where {P <: NamedSector}
    return Asymbol(map(ProductSector, (a, b, c))...)
end
function fusiontensor(a::P, b::P, c::P) where {P <: NamedSector}
    return fusiontensor(map(ProductSector, (a, b, c))...)
end

# Style traits
FusionStyle(::Type{I}) where {I <: NamedSector} = FusionStyle(_productsectortype(I))
fusionscalartype(::Type{I}) where {I <: NamedSector} = fusionscalartype(_productsectortype(I))
UnitStyle(::Type{I}) where {I <: NamedSector} = UnitStyle(_productsectortype(I))
BraidingStyle(::Type{I}) where {I <: NamedSector} = BraidingStyle(_productsectortype(I))
braidingscalartype(::Type{I}) where {I <: NamedSector} = braidingscalartype(_productsectortype(I))
sectorscalartype(::Type{I}) where {I <: NamedSector} = sectorscalartype(_productsectortype(I))
dimscalartype(::Type{I}) where {I <: NamedSector} = dimscalartype(_productsectortype(I))

fermionparity(p::NamedSector) = fermionparity(ProductSector(p))
frobenius_schur_phase(p::NamedSector) = frobenius_schur_phase(ProductSector(p))
frobenius_schur_indicator(p::NamedSector) = frobenius_schur_indicator(ProductSector(p))

dim(p::NamedSector) = dim(ProductSector(p))

# Equality, hashing, ordering
Base.isequal(p1::NamedSector, p2::NamedSector) = isequal(p1.sectors, p2.sectors)
Base.hash(p::NamedSector, h::UInt) = hash(p.sectors, h)
Base.isless(a::I, b::I) where {I <: NamedSector} = isless(ProductSector(a), ProductSector(b))

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

function type_repr(::Type{P}) where {P <: NamedSector}
    names = _sectornames(P)
    sectors = _sectors(P)
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
analogous to `@NamedTuple`.

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
