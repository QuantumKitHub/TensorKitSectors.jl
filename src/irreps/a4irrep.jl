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

function fusiontensor(a::I, b::I, c::I) where {I <: A4Irrep}
    T = sectorscalartype(I)
    Nabc = Nsymbol(a, b, c)
    C = zeros(T, dim(a), dim(b), dim(c), Nabc)
    isempty(C) && return C

    ω = cis(2π / 3)

    if a.n == b.n == 3 # 3 ⊗ 3
        if c.n != 3 # singlets
            C = A4Irrep_fusiontensor_3x3_to_1(c.n)
        else
            C = A4Irrep_fusiontensor_3x3_to_3()
        end
    else
        if a.n != 3 && b.n != 3 # 1d x 1d
            C[1, 1, 1] = one(T)
        elseif a.n == 3 && b.n != 3 # 3 x 1d
            C = A4Irrep_fusiontensor_3x1_to_3(b.n)
        else # 1d x 3
            C = reshape(A4Irrep_fusiontensor_3x1_to_3(a.n), 1, 3, 3, 1)
        end
    end
    return C
end

# TODO: check if there's an analytic expression to generate these tensors which satisfy the pentagon equation
function A4Irrep_fusiontensor_3x3_to_3()
    S = zeros(Float64, 3, 3, 3, 2)
    s2 = 1 / sqrt(2.0)
    s6 = 1 / sqrt(6.0)
    r23 = sqrt(2.0 / 3.0)

    im = (2, 1, 1)
    jm = (3, 2, 3)

    for i in 1:3
        S[im[i], jm[i], i, 1] = (i == 3) ? -s2 : s2
        S[jm[i], im[i], i, 1] = (i == 3) ? s2 : -s2
        S[im[i], jm[i], i, 2] += s6
        S[jm[i], im[i], i, 2] += s6
    end
    S[1, 1, 1, 2] = -r23
    S[3, 3, 2, 2] = -r23
    S[2, 2, 3, 2] = -r23
    return S
end

function A4Irrep_fusiontensor_3x3_to_1(n)
    C = zeros(Float64, 3, 3, 1, 1)
    sqrt3 = sqrt(3.0)
    i1s = [1, 2, 3]
    i2s = circshift(reverse(i1s), n + 1)
    is = [i1s i2s]
    for i in 1:3
        C[is[i, 1], is[i, 2], 1, 1] = 1 / sqrt3
    end
    return C
end

function A4Irrep_fusiontensor_3x1_to_3(n)
    C = zeros(Float64, 3, 1, 3, 1)
    is = circshift([1, 2, 3], n)
    for i in 1:3
        C[is[i], 1, i, 1] = 1.0
    end
    return C
end
