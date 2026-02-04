# U₁ ⋊ C (U₁ and charge conjugation)
"""
    struct CU1Irrep <: AbstractIrrep{CU₁}
    CU1Irrep(j, s = ifelse(j>zero(j), 2, 0))
    Irrep[CU₁](j, s = ifelse(j>zero(j), 2, 0))

Represents irreps of the group ``U₁ ⋊ C`` (``U₁`` and charge conjugation or reflection),
which is also known as just `O₂`. 

## Fields
- `j::HalfInt`: the value of the ``U₁`` charge.
- `s::Int`: the representation of charge conjugation.

They can take values:
*   if `j == 0`, `s = 0` (trivial charge conjugation) or
    `s = 1` (non-trivial charge conjugation)
*   if `j > 0`, `s = 2` (two-dimensional representation)
"""
struct CU1Irrep <: AbstractIrrep{CU₁}
    j::HalfInt # value of the U1 charge
    s::Int # rep of charge conjugation:
    # if j == 0, s = 0 (trivial) or s = 1 (non-trivial),
    # else s = 2 (two-dimensional representation)
    # Let constructor take the actual half integer value j
    function CU1Irrep(j::Real, s::Integer = ifelse(j > zero(j), 2, 0))
        return if ((j > zero(j) && s == 2) || (j == zero(j) && (s == 0 || s == 1)))
            new(j, s)
        else
            error("Not a valid CU₁ irrep")
        end
    end
end
Base.getindex(::IrrepTable, ::Type{CU₁}) = CU1Irrep
Base.convert(::Type{CU1Irrep}, (j, s)::Tuple{Real, Integer}) = CU1Irrep(j, s)

Base.IteratorSize(::Type{SectorValues{CU1Irrep}}) = IsInfinite()
function Base.iterate(::SectorValues{CU1Irrep}, state::Tuple{Int, Int} = (0, 0))
    j, s = state
    if iszero(j) && s == 0
        return CU1Irrep(j, s), (j, 1)
    elseif iszero(j) && s == 1
        return CU1Irrep(j, s), (j + 1, 2)
    else
        return CU1Irrep(half(j), s), (j + 1, 2)
    end
end
function Base.getindex(::SectorValues{CU1Irrep}, i::Int)
    i < 1 && throw(BoundsError(values(CU1Irrep), i))
    if i == 1
        return CU1Irrep(0, 0)
    elseif i == 2
        return CU1Irrep(0, 1)
    else
        return CU1Irrep(half(i - 2), 2)
    end
end
findindex(::SectorValues{CU1Irrep}, c::CU1Irrep) = twice(c.j) + iszero(c.j) + c.s

Base.hash(c::CU1Irrep, h::UInt) = hash(c.s, hash(c.j, h))
function Base.isless(c1::CU1Irrep, c2::CU1Irrep)
    return isless(c1.j, c2.j) || (c1.j == c2.j == zero(HalfInt) && c1.s < c2.s)
end

_tr_repr(a::CU1Irrep) = (a.j, a.s)

# CU1Irrep(j::Real, s::Int = ifelse(j>0, 2, 0)) = CU1Irrep(convert(HalfInteger, j), s)

Base.convert(::Type{CU1Irrep}, j::Real) = CU1Irrep(j)
Base.convert(::Type{CU1Irrep}, js::Tuple{Real, Int}) = CU1Irrep(js...)

unit(::Type{CU1Irrep}) = CU1Irrep(zero(HalfInt), 0)
dual(c::CU1Irrep) = c

const CU1IrrepProdIterator = SectorProductIterator{CU1Irrep}

Base.IteratorSize(::Type{CU1IrrepProdIterator}) = Base.HasLength()
function Base.length(p::CU1IrrepProdIterator)
    if p.a.j == zero(HalfInt) || p.b.j == zero(HalfInt)
        return 1
    elseif p.a == p.b
        return 3
    else
        return 2
    end
end

function Base.iterate(p::CU1IrrepProdIterator, s::Int = 1)
    if s == 1
        if p.a.j == p.b.j == zero(HalfInt)
            return CU1Irrep(zero(HalfInt), xor(p.a.s, p.b.s)), 4
        elseif p.a.j == zero(HalfInt)
            return p.b, 4
        elseif p.b.j == zero(HalfInt)
            return p.a, 4
        elseif p.a == p.b # != zero
            return unit(CU1Irrep), 2
        else
            return CU1Irrep(abs(p.a.j - p.b.j)), 3
        end
    elseif s == 2
        return CU1Irrep(zero(HalfInt), 1), 3
    elseif s == 3
        CU1Irrep(p.a.j + p.b.j), 4
    else
        return nothing
    end
end

dim(c::CU1Irrep) = ifelse(c.j == zero(HalfInt), 1, 2)

FusionStyle(::Type{CU1Irrep}) = SimpleFusion()
fusionscalartype(::Type{CU1Irrep}) = Float64
braidingscalartype(::Type{CU1Irrep}) = Float64
sectorscalartype(::Type{CU1Irrep}) = Float64
Base.isreal(::Type{CU1Irrep}) = true

function Nsymbol(a::CU1Irrep, b::CU1Irrep, c::CU1Irrep)
    return ifelse(
        c.s == 0, (a.j == b.j) & ((a.s == b.s == 2) | (a.s == b.s)),
        ifelse(
            c.s == 1, (a.j == b.j) & ((a.s == b.s == 2) | (a.s != b.s)),
            (c.j == a.j + b.j) | (c.j == abs(a.j - b.j))
        )
    )
end

function Fsymbol(
        a::CU1Irrep, b::CU1Irrep, c::CU1Irrep,
        d::CU1Irrep, e::CU1Irrep, f::CU1Irrep
    )
    Nabe = convert(Int, Nsymbol(a, b, e))
    Necd = convert(Int, Nsymbol(e, c, d))
    Nbcf = convert(Int, Nsymbol(b, c, f))
    Nafd = convert(Int, Nsymbol(a, f, d))

    T = fusionscalartype(CU1Irrep)
    Nabe * Necd * Nbcf * Nafd == 0 && return zero(T)

    op = CU1Irrep(0, 0)
    om = CU1Irrep(0, 1)

    if a == op || b == op || c == op
        return one(T)
    end
    if (a == b == om) || (a == c == om) || (b == c == om)
        return one(T)
    end
    if a == om
        if d.j == zero(HalfInt)
            return one(T)
        else
            return (d.j == c.j - b.j) ? -one(T) : one(T)
        end
    end
    if b == om
        return (d.j == abs(a.j - c.j)) ? -one(T) : one(T)
    end
    if c == om
        return (d.j == a.j - b.j) ? -one(T) : one(T)
    end
    # from here on, a, b, c are neither 0+ or 0-
    s = T(sqrt(2) / 2)
    if a == b == c
        if d == a
            if e.j == 0
                if f.j == 0
                    return f.s == 1 ? T(-1 // 2) : T(1 // 2)
                else
                    return e.s == 1 ? -s : s
                end
            else
                return f.j == 0 ? s : zero(T)
            end
        else
            return one(T)
        end
    end
    if a == b # != c
        if d == c
            if f.j == b.j + c.j
                return e.s == 1 ? -s : s
            else
                return s
            end
        else
            return one(T)
        end
    end
    if b == c
        if d == a
            if e.j == a.j + b.j
                return s
            else
                return f.s == 1 ? -s : s
            end
        else
            return one(T)
        end
    end
    if a == c
        if d == b
            if e.j == f.j
                return zero(T)
            else
                return one(T)
            end
        else
            return d.s == 1 ? -one(T) : one(T)
        end
    end
    if d == om
        return b.j == a.j + c.j ? -one(T) : one(T)
    end
    return one(T)
end

function Rsymbol(a::CU1Irrep, b::CU1Irrep, c::CU1Irrep)
    R = convert(braidingscalartype(CU1Irrep), Nsymbol(a, b, c))
    return c.s == 1 && a.j > 0 ? -R : R
end

function fusiontensor(a::CU1Irrep, b::CU1Irrep, c::CU1Irrep)
    C = fill(zero(sectorscalartype(CU1Irrep)), dim(a), dim(b), dim(c), 1)
    !Nsymbol(a, b, c) && return C
    if c.j == 0
        if a.j == b.j == 0
            C[1, 1, 1, 1] = 1.0
        else
            if c.s == 0
                C[1, 2, 1, 1] = 1.0 / sqrt(2)
                C[2, 1, 1, 1] = 1.0 / sqrt(2)
            else
                C[1, 2, 1, 1] = 1.0 / sqrt(2)
                C[2, 1, 1, 1] = -1.0 / sqrt(2)
            end
        end
    elseif a.j == 0
        C[1, 1, 1, 1] = 1.0
        C[1, 2, 2, 1] = a.s == 1 ? -1.0 : 1.0
    elseif b.j == 0
        C[1, 1, 1, 1] = 1.0
        C[2, 1, 2, 1] = b.s == 1 ? -1.0 : 1.0
    elseif c.j == a.j + b.j
        C[1, 1, 1, 1] = 1.0
        C[2, 2, 2, 1] = 1.0
    elseif c.j == a.j - b.j
        C[1, 2, 1, 1] = 1.0
        C[2, 1, 2, 1] = 1.0
    elseif c.j == b.j - a.j
        C[2, 1, 1, 1] = 1.0
        C[1, 2, 2, 1] = 1.0
    end
    return C
end
frobenius_schur_phase(::CU1Irrep) = 1
