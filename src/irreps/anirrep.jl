"""
    struct ANIrrep{N} <: AbstractIrrep{Alternating{N}}
    ANIrrep{N}(n::Integer)
    Irrep[Alternating{N}](n::Integer)

Represents irreps of the alternating group ``A_N``.

## Fields

- `n::Int`: the label of the irrep.

This can take the value 
"""
struct ANIrrep{N} <: AbstractIrrep{Alternating{N}}
    n::Int
    function ANIrrep{N}(n::Int) where {N}
        @assert 0 < N < 5 "ANIrrep is only provided for 0 < N < 5"
        0 <= n <= _npartitions(N) - 1 ||
        throw(DomainError(n, "ANIrrep only has irreps " * lazy"0 <= n <= $(_npartitions(N) - 1)"))
        return new{N}(n)
    end
end

_npartitions(N::Int) = (N == 3) ? 3 : (N == 4) ? 4 : 1

Base.propertynames(x::ANIrrep) = (:n,)

function Base.getproperty(a::ANIrrep{N}, sym::Symbol) where {N}
    if sym === :n
        return getfield(a, :n)
    else
        error("Unknown property $sym")
    end
end

FusionStyle(::Type{ANIrrep{N}}) where {N} = N < 4 ? UniqueFusion() : GenericFusion()
sectorscalartype(::Type{ANIrrep{N}}) where {N} = N < 4 ? Int64 : Float64
Base.isreal(::Type{ANIrrep{N}}) where {N} = true

unit(::Type{ANIrrep{N}}) where {N} = ANIrrep{N}(0)
function dual(a::ANIrrep{N}) where {N}
    N < 3 && return a
    if N == 3
        return ANIrrep{3}((3 - a.n) % 3)
    else # N = 4
        return a.n == 3 ? ANIrrep{4}(3) : ANIrrep{4}((3 - a.n) % 3)
    end
end

Base.hash(a::ANIrrep, h::UInt) = hash(a.n, h)
Base.convert(::Type{ANIrrep{N}}, n::Integer) where {N} = ANIrrep{N}(n)

Base.getindex(::IrrepTable, ::Type{Alternating{N}}) where {N} = ANIrrep{N}

const A4Irrep = ANIrrep{4}

# Sector iterator
# ---------------
Base.isless(a::ANIrrep{N}, b::ANIrrep{N}) where {N} = isless(a.n, b.n)
Base.IteratorSize(::Type{SectorValues{ANIrrep{N}}}) where {N} = Base.HasLength()
Base.length(::SectorValues{ANIrrep{N}}) where {N} = _npartitions(N)

function Base.iterate(v::SectorValues{ANIrrep{N}}, i = 1) where {N}
    return i > length(v) ? nothing : (v[i], i + 1)
end

@inline function Base.getindex(v::SectorValues{ANIrrep{N}}, i::Int) where {N}
    L = length(v)
    @boundscheck 1 <= i <= L || throw(BoundsError(v, i))
    return ANIrrep{N}(i - 1)
end

findindex(::SectorValues{ANIrrep{N}}, a::ANIrrep{N}) where {N} = a.n + 1

# Product iterator
# ----------------

const ANIrrepProdIterator{N} = SectorProductIterator{ANIrrep{N}}
⊗(a::ANIrrep{N}, b::ANIrrep{N}) where {N} = SectorProductIterator((a <= b ? (a, b) : (b, a))...)

Base.IteratorSize(::Type{ANIrrepProdIterator{N}}) where {N} = Base.HasLength()
function Base.length(x::ANIrrepProdIterator{N}) where {N}
    N < 4 && return 1 # abelian
    a, b = x.a, x.b
    if a == b
        return a.n == 3 ? 4 : 1 # 3 x 3 = 1 + 1' + 1'' + 2*3
    else
        return 1
    end
end

# TODO: any way to shorten this?
function Base.iterate(p::ANIrrepProdIterator{N}, state::Int = 1) where {N}
    a, b = p.a, p.b
    if N < 4
        return state > 1 ? nothing : (ANIrrep{N}((a.n + b.n) % _npartitions(N)), state + 1)
    else
        u, u′, u′′, three = ANIrrep{4}(0), ANIrrep{4}(1), ANIrrep{4}(2), ANIrrep{4}(3)
        if state == 1
            # rules with the trivial irrep
            a.n == 0 && return (b, 2)
            b.n == 0 && return (a, 2)
            # 1' rules
            if a.n + b.n == 3 # (1, 2) or (2, 1)
                return (u, 2)
            elseif a.n == b.n != 3 # (1,1) or (2,2)
                return (A4Irrep((2 * a.n) % 3), 2)
            elseif (a.n in (1, 2) && b.n == 3) || (a.n == 3 && b.n in (1, 2))
                # all 1D × 3 = 3
                return (three, 2)
            else # 3 x 3
                return (u, 2) # first in sequence
            end
        elseif state == 2 && a.n == 3 && b.n == 3
            return (u′, 3)
        elseif state == 3 && a.n == 3 && b.n == 3
            return (u′′, 4)
        elseif state == 4 && a.n == 3 && b.n == 3
            return (three, 5)
        else
            return nothing
        end
    end
end

# Topological data
# ----------------
function dim(a::ANIrrep{N}) where {N}
    N < 4 && return 1
    return a.n == 3 ? 3 : 1
end

function Nsymbol(a::ANIrrep{N}, b::ANIrrep{N}, c::ANIrrep{N}) where {N}
    N < 4 && return (c.n == (a.n + b.n) % _npartitions(N)) ? 1 : 0
    # N = 4
    u, u′, u′′, three = ANIrrep{4}(0), ANIrrep{4}(1), ANIrrep{4}(2), ANIrrep{4}(3)
    if a == three
        if b == three # 3 x 3
            return (c == u || c == u′ || c == u′′) ? 1 : 2
        else # 3 x 1D
            return c == three ? 1 : 0
        end
    else
        if b == three # 1D x 3
            return c == three ? 1 : 0
        else # 1D x 1D
            return c.n == (a.n + b.n) % 3 ? 1 : 0
        end
    end
end

function Fsymbol(a::I, b::I, c::I, d::I, e::I, f::I) where {N, I <: ANIrrep{N}}
    T = sectorscalartype(I)
    Nabe = Nsymbol(a, b, e)
    Necd = Nsymbol(e, c, d)
    Nbcf = Nsymbol(b, c, f)
    Nafd = Nsymbol(a, f, d)

    Nabe > 0 && Necd > 0 && Nbcf > 0 && Nafd > 0 ||
        return zeros(sectorscalartype(I), Nabe, Necd, Nbcf, Nafd)

    # fallback through fusiontensor
    A = fusiontensor(a, b, e)
    B = fusiontensor(e, c, d)[:, :, 1, :]
    C = fusiontensor(b, c, f)
    D = fusiontensor(a, f, d)[:, :, 1, :]
    @tensor F[-1, -2, -3, -4] := conj(D[1, 5, -4]) * conj(C[2, 4, 5, -3]) *
        A[1, 2, 3, -1] * B[3, 4, -2]

    return F
end

# bosonic
function Rsymbol(a::I, b::I, c::I) where {N, I <: ANIrrep{N}}
    R = convert(sectorscalartype(I), Nsymbol(a, b, c))
    return R
end

function fusiontensor(a::I, b::I, c::I) where {N, I <: ANIrrep{N}}
    T = sectorscalartype(I)
    Nabc = Nsymbol(a, b, c)
    C = zeros(T, dim(a), dim(b), dim(c), Nabc)
    Nabc == 0 && return C
    N < 4 && return ones(T, 1, 1, 1, 1) # all fusion channels trivial

    ω = cis(2π/3)

    if a.n == b.n == 3 # 3 ⊗ 3
        if c.n != 3 # singlets
            C = _fusiontensor_3x3_to_1(c.n)
        # if c.n == 0
        #     for i in 1:3
        #         C[i,i,1] = 1/sqrt(3)
        #     end
        # elseif c.n == 1
        #     for i in 1:3
        #         C[i,i,1] = ω^(i-1)/sqrt(3)
        #     end
        # elseif c.n == 2
        #     for i in 1:3
        #         C[i,i,1] = ω^(2*(i-1))/sqrt(3)
        #     end
        else
            C = _fusiontensor_3x3_to_3()
            # im = (2, 3, 1) # cyclic pairs
            # jm = (3, 1, 2)
            # delta(i, j) = i == j
            # for m in 1:3
            #     i_m, j_m = im[m], jm[m]
            #     for i in 1:3, j in 1:3
            #         C[m, i, j, 1] = 1/sqrt(2) * (delta(i, i_m) * delta(j, j_m) - delta(i, j_m) * delta(j, i_m))

            #         C[m, i, j, 2] = 1/sqrt(6) * (delta(i, i_m) * delta(j, j_m) + delta(i, j_m) * delta(j, i_m) - 2 * delta(i, m) * delta(j, m))
            #     end
            # end
        end
    
    else
        if a.n != 3 && b.n != 3 # 1d x 1d
            C[1,1,1] = one(T)
        elseif a.n == 3 && b.n != 3 # 3 x 1d
            C = _fusiontensor_3x1_to_3(b.n)
        else # 1d x 3
            C = reshape(_fusiontensor_3x1_to_3(a.n), 1, 3, 3, 1)
        end
    end

    return C
end

# TODO: for some reason the analytic expression doesn't match these results, which is from CategoryData
function _fusiontensor_3x3_to_3()
    S = zeros(Float64, 3, 3, 3, 2)
    s2 = 1 / sqrt(2.0)
    s6 = 1 / sqrt(6.0)
    r23 = sqrt(2.0/3.0)

    im  = (2, 1, 1)
    jm = (3, 2, 3)

    # μ = 1 : antisymmetric off-diagonal entries
    S[im[1], jm[1], 1, 1] =  s2
    S[jm[1], im[1], 1, 1] = -s2
    S[im[2], jm[2], 2, 1] =  s2
    S[jm[2], im[2], 2, 1] = -s2
    S[im[3], jm[3], 3, 1] =  -s2
    S[jm[3], im[3], 3, 1] = s2

    # μ = 2 : diagonal negative + symmetric off-diagonals
    S[1, 1, 1, 2] = -r23
    S[im[1], jm[1], 1, 2] += s6
    S[jm[1], im[1], 1, 2] += s6

    S[3, 3, 2, 2] = -r23
    S[im[2], jm[2], 2, 2] += s6
    S[jm[2], im[2], 2, 2] += s6

    S[2, 2, 3, 2] = -r23
    S[im[3], jm[3], 3, 2] += s6
    S[jm[3], im[3], 3, 2] += s6

    return S
end

function _fusiontensor_3x3_to_1(n::Int)
    C = zeros(Float64, 3, 3, 1, 1)
    sqrt3 = sqrt(3.0)
    ijs = Vector{Tuple{Int,Int}}(undef, 3)
    if n == 0
        ijs[1] = (1, 1)
        ijs[2] = (2, 3)
        ijs[3] = (3, 2)
    elseif n == 1
        ijs[1] = (1, 2)
        ijs[2] = (2, 1)
        ijs[3] = (3, 3)
    elseif n == 2
        ijs[1] = (1, 3)
        ijs[2] = (2, 2)
        ijs[3] = (3, 1)
    end
    for i in 1:3
        C[ijs[i][1], ijs[i][2], 1, 1] = 1 / sqrt3
    end
    return C
end

function _fusiontensor_3x1_to_3(n::Int)
    C = zeros(Float64, 3, 1, 3, 1)
    ijs = [(1, 2, 3), (3, 1, 2), (2, 3, 1)]
    _ijs = ijs[n + 1]
    for i in 1:3
        C[_ijs[i], 1, i, 1] = 1.0
    end
    return C
end