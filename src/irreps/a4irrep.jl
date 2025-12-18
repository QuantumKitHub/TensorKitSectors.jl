"""
    struct A4Irrep <: AbstractIrrep{A₄}
    A4Irrep(n::Integer)
    Irrep[A₄](n::Integer)

Represents irreps of the alternating group ``A₄``.

## Fields

- `n::Int8`: the label of the irrep, corresponding to ``1``, ``1′``, ``1″`` and ``3``.
"""
struct A4Irrep <: AbstractIrrep{A₄}
    n::Int8
    function A4Irrep(n::Integer)
        0 ≤ n < 4 || throw(ArgumentError("A4Irrep only has irrep labels in `0:3`"))
        return new(n)
    end
end

FusionStyle(::Type{A4Irrep}) = GenericFusion()
FusionDataStyle(::Type{A4Irrep}) = NonTrivialFusionData()
sectorscalartype(::Type{A4Irrep}) = Float64

unit(::Type{A4Irrep}) = A4Irrep(0)
dual(a::A4Irrep) = a.n == 3 ? a : A4Irrep((3 - a.n) % 3)

Base.hash(a::A4Irrep, h::UInt) = hash(a.n, h)
Base.convert(::Type{A4Irrep}, n::Integer) = A4Irrep(n)

Base.getindex(::IrrepTable, ::Type{A₄}) = A4Irrep


# Sector iterator
# ---------------
Base.isless(a::A4Irrep, b::A4Irrep) = isless(a.n, b.n)
Base.IteratorSize(::Type{SectorValues{A4Irrep}}) = Base.HasLength()
Base.length(::SectorValues{A4Irrep}) = 4

Base.iterate(v::SectorValues{A4Irrep}, i = 1) = i > length(v) ? nothing : (v[i], i + 1)

@inline function Base.getindex(v::SectorValues{A4Irrep}, i::Int)
    @boundscheck 1 <= i <= length(v) || throw(BoundsError(v, i))
    return A4Irrep(i - 1)
end

findindex(::SectorValues{A4Irrep}, a::A4Irrep) = a.n + 1

# Product iterator
# ----------------

const A4IrrepProdIterator = SectorProductIterator{A4Irrep}
⊗(a::A4Irrep, b::A4Irrep) = SectorProductIterator((a <= b ? (a, b) : (b, a))...)

Base.IteratorSize(::Type{A4IrrepProdIterator}) = Base.HasLength()
Base.length(x::A4IrrepProdIterator) = (x.a == x.b == A4Irrep(3)) ? 4 : 1

function Base.iterate(p::A4IrrepProdIterator, state::Int = 1)
    a, b = p.a, p.b
    if state == 1
        a.n == b.n == 3 && return (A4Irrep(state - 1), state + 1) # 3 ⊗ 3
        b.n == 3 && return (b, 5) # x ⊗ 3 = 3
        return (A4Irrep((a.n + b.n) % 3), 5) # 1d irreps ≡ Z3
    elseif state < 5 # a == b == 3
        return (A4Irrep(state - 1), state + 1)
    else
        return nothing
    end
end

# Topological data
# ----------------
dim(a::A4Irrep) = a.n == 3 ? 3 : 1

function Nsymbol(a::A4Irrep, b::A4Irrep, c::A4Irrep)
    # 3d irreps
    a.n == b.n == 3 && return 1 + (c.n == 3)
    (a.n == 3 || b.n == 3) && return Int(c.n == 3)
    # 1d irreps
    return Int((a.n + b.n) % 3 == c.n)
end


function Fsymbol(a::I, b::I, c::I, d::I, e::I, f::I) where {I <: A4Irrep}
    T = sectorscalartype(I)
    Nabe = Nsymbol(a, b, e)
    Necd = Nsymbol(e, c, d)
    Nbcf = Nsymbol(b, c, f)
    Nafd = Nsymbol(a, f, d)

    Nabe > 0 && Necd > 0 && Nbcf > 0 && Nafd > 0 ||
        return zeros(T, Nabe, Necd, Nbcf, Nafd)

    # fallback through fusiontensor for A4Irrep
    A = fusiontensor(a, b, e)
    B = fusiontensor(e, c, d)[:, :, 1, :]
    C = fusiontensor(b, c, f)
    D = fusiontensor(a, f, d)[:, :, 1, :]
    @tensor F[-1, -2, -3, -4] := conj(D[1, 5, -4]) * conj(C[2, 4, 5, -3]) *
        A[1, 2, 3, -1] * B[3, 4, -2]

    return F
end

# bosonic
function Rsymbol(a::I, b::I, c::I) where {I <: A4Irrep}
    Nabc = Nsymbol(a, b, c)
    R = zeros(sectorscalartype(I), Nabc, Nabc)
    Nabc == 0 && return R
    if a == b == c == A4Irrep(3)
        R[1, 1] = -1
        R[2, 2] = 1
    else
        R[1, 1] = 1
    end
    return R
end

# choice of basis: https://journals.aps.org/rmp/pdf/10.1103/RevModPhys.82.2701
# triplet is a real representation -> can make all representation matrices real
function fusiontensor(a::I, b::I, c::I) where {I <: A4Irrep}
    T = sectorscalartype(I)
    Nabc = Nsymbol(a, b, c)
    C = zeros(T, dim(a), dim(b), dim(c), Nabc)
    isempty(C) && return C

    ω = cis(2π / 3)

    if a.n == b.n == 3 # 3 ⊗ 3
        if c.n != 3 # singlets
            invsqrt3 = 1 / sqrt(3.0)
            for i in 1:3
                j = 4 - mod1(i - c.n - 1, 3)
                C[i, j, 1, 1] = invsqrt3
            end
        else # triplet: eq 38 in above reference
            s2 = 1 / sqrt(2.0)
            s6 = 1 / sqrt(6.0)

            # antisymmetric channel
            C[:, :, 1, 1] .= [0 0 0; 0 0 s2; 0 -s2 0]
            C[:, :, 2, 1] .= [0 s2 0; -s2 0 0; 0 0 0]
            C[:, :, 3, 1] .= [0 0 -s2; 0 0 0; s2 0 0]

            # symmetric channel
            C[:, :, 1, 2] .= [-2 * s6 0 0; 0 0 s6; 0 s6 0]
            C[:, :, 2, 2] .= [0 s6 0; s6 0 0; 0 0 -2 * s6]
            C[:, :, 3, 2] .= [0 0 s6; 0 -2 * s6 0; s6 0 0]
        end
    else
        if a.n != 3 && b.n != 3 # 1d x 1d
            C[1, 1, 1] = one(T)
        elseif a.n == 3 && b.n != 3 # 3 x 1d
            for i in 1:3
                j = mod1(i - b.n, 3)
                C[j, 1, i, 1] = one(T)
            end
        else # 1d x 3: reshape of 3 x 1d
            for i in 1:3
                j = mod1(i - a.n, 3)
                C[1, j, i, 1] = one(T)
            end
        end
    end
    return C
end
