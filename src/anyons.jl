# Sectors corresponding to simplest fusion categories that are not group representations.
#---------------------------------------------------------------------------------------#

# PlanarTrivial
"""
    struct PlanarTrivial <: Sector
    PlanarTrivial()

Represents a trivial anyon sector, i.e. a trivial sector without braiding. This is mostly
useful for testing.
"""
struct PlanarTrivial <: Sector end

type_repr(::Type{PlanarTrivial}) = "PlanarTrivial"

Base.IteratorSize(::Type{SectorValues{PlanarTrivial}}) = HasLength()
Base.length(::SectorValues{PlanarTrivial}) = 1
Base.iterate(::SectorValues{PlanarTrivial}, i = 0) = i == 0 ? (PlanarTrivial(), 1) : nothing
function Base.getindex(::SectorValues{PlanarTrivial}, i::Int)
    return i == 1 ? PlanarTrivial() :
        throw(BoundsError(values(PlanarTrivial), i))
end
findindex(::SectorValues{PlanarTrivial}, c::PlanarTrivial) = 1
Base.isless(::PlanarTrivial, ::PlanarTrivial) = false

unit(::Type{PlanarTrivial}) = PlanarTrivial()
dual(::PlanarTrivial) = PlanarTrivial()

FusionStyle(::Type{PlanarTrivial}) = UniqueFusion()
BraidingStyle(::Type{PlanarTrivial}) = NoBraiding()
fusionscalartype(::Type{PlanarTrivial}) = Int
sectorscalartype(::Type{PlanarTrivial}) = Int

Nsymbol(::Vararg{PlanarTrivial, 3}) = 1
Fsymbol(::Vararg{PlanarTrivial, 6}) = 1

⊗(::PlanarTrivial, ::PlanarTrivial) = (PlanarTrivial(),)

# FibonacciAnyons
"""
    struct FibonacciAnyon <: Sector
    FibonacciAnyon(s::Symbol)

Represents the anyons of the Fibonacci modular fusion category. It can take two values,
corresponding to the trivial sector `FibonacciAnyon(:I)` and the non-trivial sector
`FibonacciAnyon(:τ)` with fusion rules ``τ ⊗ τ = 1 ⊕ τ``.

## Fields
- `isunit::Bool`: indicates whether the sector corresponds to the trivial anyon `:I`
  (`true`), or the non-trivial anyon `:τ` (`false`).
"""
struct FibonacciAnyon <: Sector
    isunit::Bool
    function FibonacciAnyon(s::Symbol)
        s in (:I, :τ, :tau) || throw(ArgumentError("Unknown FibonacciAnyon $s."))
        return new(s === :I)
    end
end

type_repr(::Type{FibonacciAnyon}) = "FibonacciAnyon"

Base.IteratorSize(::Type{SectorValues{FibonacciAnyon}}) = HasLength()
Base.length(::SectorValues{FibonacciAnyon}) = 2
function Base.iterate(::SectorValues{FibonacciAnyon}, i = 0)
    return i == 0 ? (FibonacciAnyon(:I), 1) : (i == 1 ? (FibonacciAnyon(:τ), 2) : nothing)
end
function Base.getindex(S::SectorValues{FibonacciAnyon}, i::Int)
    if i == 1
        return FibonacciAnyon(:I)
    elseif i == 2
        return FibonacciAnyon(:τ)
    else
        throw(BoundsError(S, i))
    end
end
findindex(::SectorValues{FibonacciAnyon}, s::FibonacciAnyon) = 2 - s.isunit

Base.convert(::Type{FibonacciAnyon}, s::Symbol) = FibonacciAnyon(s)
unit(::Type{FibonacciAnyon}) = FibonacciAnyon(:I)
dual(s::FibonacciAnyon) = s

const _goldenratio = Float64(MathConstants.golden)
dim(a::FibonacciAnyon) = isunit(a) ? one(_goldenratio) : _goldenratio

FusionStyle(::Type{FibonacciAnyon}) = SimpleFusion()
fusionscalartype(::Type{FibonacciAnyon}) = Float64
BraidingStyle(::Type{FibonacciAnyon}) = Anyonic()
braidingscalartype(::Type{FibonacciAnyon}) = ComplexF64
sectorscalartype(::Type{FibonacciAnyon}) = ComplexF64

const FibonacciAnyonProdIterator = SectorProductIterator{FibonacciAnyon}

Base.IteratorSize(::Type{FibonacciAnyonProdIterator}) = Base.HasLength()
Base.length(iter::FibonacciAnyonProdIterator) = (isunit(iter.a) || isunit(iter.b)) ? 1 : 2

function Base.iterate(iter::FibonacciAnyonProdIterator, state = 1)
    I = FibonacciAnyon(:I)
    τ = FibonacciAnyon(:τ)
    if state == 1 # first iteration
        iter.a == I && return (iter.b, 2)
        iter.b == I && return (iter.a, 2)
        return (I, 2)
    elseif state == 2
        (iter.a == iter.b == τ) && return (τ, 3)
        return nothing
    else
        return nothing
    end
end

function Nsymbol(a::FibonacciAnyon, b::FibonacciAnyon, c::FibonacciAnyon)
    return isunit(a) + isunit(b) + isunit(c) != 2
end # zero if one tau and two ones

function Fsymbol(
        a::FibonacciAnyon, b::FibonacciAnyon, c::FibonacciAnyon,
        d::FibonacciAnyon, e::FibonacciAnyon, f::FibonacciAnyon
    )
    Nsymbol(a, b, e) || return zero(_goldenratio)
    Nsymbol(e, c, d) || return zero(_goldenratio)
    Nsymbol(b, c, f) || return zero(_goldenratio)
    Nsymbol(a, f, d) || return zero(_goldenratio)

    I = FibonacciAnyon(:I)
    τ = FibonacciAnyon(:τ)
    if a == b == c == d == τ
        if e == f == I
            return +1 / _goldenratio
        elseif e == f == τ
            return -1 / _goldenratio
        else
            return +1 / sqrt(_goldenratio)
        end
    else
        return one(_goldenratio)
    end
end

function Rsymbol(a::FibonacciAnyon, b::FibonacciAnyon, c::FibonacciAnyon)
    Nsymbol(a, b, c) || return 0 * cispi(0 / 1)
    if isunit(a) || isunit(b)
        return cispi(0 / 1)
    else
        return isunit(c) ? cispi(4 / 5) : cispi(-3 / 5)
    end
end

function Base.show(io::IO, a::FibonacciAnyon)
    s = isunit(a) ? ":I" : ":τ"
    return get(io, :typeinfo, nothing) === FibonacciAnyon ?
        print(io, s) : print(io, "FibonacciAnyon(", s, ")")
end

Base.hash(a::FibonacciAnyon, h::UInt) = hash(a.isunit, h)
Base.isless(a::FibonacciAnyon, b::FibonacciAnyon) = isless(!a.isunit, !b.isunit)

# IsingAnyons
"""
    struct IsingAnyon <: Sector
    IsingAnyon(s::Symbol)

Represents the anyons of the Ising modular fusion category. It can take three values,
corresponding to the trivial sector `IsingAnyon(:I)` and the non-trivial sectors
`IsingAnyon(:σ)` and `IsingAnyon(:ψ)`, with fusion rules ``ψ ⊗ ψ = 1``, ``σ ⊗ ψ = σ``, and
``σ ⊗ σ = 1 ⊕ ψ``.

## Fields
- `label::UInt8`: internal encoding of the anyon, with `0`, `1`, `2` corresponding to `:I`,
  `:σ`, and `:ψ` respectively. The symbol label can be obtained as `a.s`.
"""
struct IsingAnyon <: Sector
    label::UInt8
    function IsingAnyon(s::Symbol)
        label = if s === :I
            0x00
        elseif s === :σ || s === :sigma
            0x01
        elseif s === :ψ || s === :psi
            0x02
        else
            throw(ArgumentError("Unknown IsingAnyon $s."))
        end
        return new(label)
    end
    function IsingAnyon(label::UInt8)
        label <= 0x02 || throw(ArgumentError("Unknown IsingAnyon label $label."))
        return new(label)
    end
end

const _isinganyon_symbols = (:I, :σ, :ψ)

function Base.getproperty(a::IsingAnyon, f::Symbol)
    f === :s && return _isinganyon_symbols[getfield(a, :label) + 1]
    return getfield(a, f)
end
Base.propertynames(::IsingAnyon, private::Bool = false) =
    private ? (:s, :label) : (:s,)

type_repr(::Type{IsingAnyon}) = "IsingAnyon"

# with the encoding I = 0, σ = 1, ψ = 2, σ is the only odd label, and the abelian
# subcategory {I, ψ} ≅ Z2 fuses as `xor` of the labels.
_isσ(a::IsingAnyon) = isodd(a.label)

Base.IteratorSize(::Type{SectorValues{IsingAnyon}}) = HasLength()
Base.length(::SectorValues{IsingAnyon}) = 3
function Base.iterate(::SectorValues{IsingAnyon}, label::UInt8 = 0x00)
    return label > 0x02 ? nothing : (IsingAnyon(label), label + 0x01)
end
function Base.getindex(vals::SectorValues{IsingAnyon}, i::Int)
    checkbounds(Bool, 1:length(vals), i) || throw(BoundsError(vals, i))
    return IsingAnyon(i % UInt8 - 0x01)
end

findindex(::SectorValues{IsingAnyon}, a::IsingAnyon) = Int(a.label) + 1

Base.convert(::Type{IsingAnyon}, s::Symbol) = IsingAnyon(s)
unit(::Type{IsingAnyon}) = IsingAnyon(:I)
dual(s::IsingAnyon) = s

dim(a::IsingAnyon) = _isσ(a) ? sqrt(2) : 1.0

FusionStyle(::Type{IsingAnyon}) = SimpleFusion()
fusionscalartype(::Type{IsingAnyon}) = Float64
BraidingStyle(::Type{IsingAnyon}) = Anyonic()
braidingscalartype(::Type{IsingAnyon}) = ComplexF64
sectorscalartype(::Type{IsingAnyon}) = ComplexF64

const IsingAnyonProdIterator = SectorProductIterator{IsingAnyon}

Base.IteratorSize(::Type{IsingAnyonProdIterator}) = Base.HasLength()
Base.length(iter::IsingAnyonProdIterator) = 1 + (_isσ(iter.a) & _isσ(iter.b))

function Base.iterate(iter::IsingAnyonProdIterator, state = 1)
    a, b = iter.a, iter.b
    if state == 1
        _isσ(a) & _isσ(b) && return (IsingAnyon(:I), 2)       # σ ⊗ σ ∋ 1
        _isσ(a) | _isσ(b) && return (IsingAnyon(:σ), 2)       # σ ⊗ (1|ψ) = σ
        return (IsingAnyon(xor(a.label, b.label)), 2)        # Z2 fusion of {1, ψ}
    elseif state == 2 && _isσ(a) & _isσ(b)
        return (IsingAnyon(:ψ), 3)                           # σ ⊗ σ ∋ ψ
    else
        return nothing
    end
end

function Nsymbol(a::IsingAnyon, b::IsingAnyon, c::IsingAnyon)
    # either exactly two σ's, or no σ's and the Z2 charges cancel
    nσ = _isσ(a) + _isσ(b) + _isσ(c)
    return (nσ == 2) | ((nσ == 0) & iszero(xor(a.label, b.label, c.label)))
end

function Fsymbol(
        a::IsingAnyon, b::IsingAnyon, c::IsingAnyon,
        d::IsingAnyon, e::IsingAnyon, f::IsingAnyon
    )
    Nsymbol(a, b, e) || return 0.0
    Nsymbol(e, c, d) || return 0.0
    Nsymbol(b, c, f) || return 0.0
    Nsymbol(a, f, d) || return 0.0
    σ, ψ = IsingAnyon(:σ), IsingAnyon(:ψ)
    if a == b == c == d == σ
        if e == f == ψ
            return -1.0 / sqrt(2.0)
        else
            return 1.0 / sqrt(2.0)
        end
    end
    if e == f == σ
        if a == c == σ && b == d == ψ
            return -1.0
        elseif a == c == ψ && b == d == σ
            return -1.0
        end
    end
    return 1.0
end

function Rsymbol(a::IsingAnyon, b::IsingAnyon, c::IsingAnyon)
    Nsymbol(a, b, c) || return complex(0.0)
    I, σ, ψ = IsingAnyon(:I), IsingAnyon(:σ), IsingAnyon(:ψ)
    if c == I
        if b == a == σ
            return cispi(-1 / 8)
        elseif b == a == ψ
            return complex(-1.0)
        end
    elseif c == σ && (a == σ && b == ψ || a == ψ && b == σ)
        return -1.0im
    elseif c == ψ && a == b == σ
        return cispi(3 / 8)
    end
    return complex(1.0)
end

function Base.show(io::IO, a::IsingAnyon)
    if get(io, :typeinfo, nothing) === IsingAnyon
        return print(io, ":$(a.s)")
    else
        return print(io, "IsingAnyon(:$(a.s))")
    end
end

Base.hash(a::IsingAnyon, h::UInt) = hash(a.label, h)
Base.isless(a::IsingAnyon, b::IsingAnyon) = isless(a.label, b.label)
