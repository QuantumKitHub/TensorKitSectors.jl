"""
    struct HeisenbergIrrep{N} <: AbstractIrrep{Heisenberg{N}}
    HeisenbergIrrep{N}(n::Int8)
    Irrep[Heisenberg{N}](n::Int8)

Represents irreps of the finite Heisenberg group ``H_N`` for `N` prime.
Irreps are labeled by a single integer `n`: if `N < n < 0`, the irrep is the Schrödinger
representation ``πₖ`` with `k = -n`. Otherwise, the irrep is the character ``χₐ,ᵦ`` with
`a, b = divrem(n, N)`, of which there are N².

## Fields

- `n::Int`: the label of the irrep.
"""
struct HeisenbergIrrep{N} <: AbstractIrrep{Heisenberg{N}}
    n::Int8
    function HeisenbergIrrep{N}(n::Integer) where {N}
        if !isprime(N)
            throw(ArgumentError("N must be a prime number"))
        end
        -N < n < N^2 || throw(ArgumentError("Unknown HeisenbergIrrep{N} with label n = $n."))
        return new{N}(n)
    end
end

FusionStyle(::Type{<:HeisenbergIrrep}) = GenericFusion()
BraidingStyle(::Type{<:HeisenbergIrrep}) = Anyonic()
fusionscalartype(::Type{<:HeisenbergIrrep}) = ComplexF64
braidingscalartype(::Type{<:HeisenbergIrrep}) = ComplexF64
sectorscalartype(::Type{<:HeisenbergIrrep}) = ComplexF64

isschrodinger(a::HeisenbergIrrep) = a.n < 0

unit(::Type{HeisenbergIrrep{N}}) where {N} = HeisenbergIrrep{N}(0)
function dual(x::HeisenbergIrrep{N}) where {N}
    if isschrodinger(x) # Schrödinger irrep: k -> N - k
        return HeisenbergIrrep{N}(-N - x.n)
    else # character: (a, b) -> (-a, -b)
        a, b = divrem(x.n, N)
        return HeisenbergIrrep{N}(mod(-a, N) * N + mod(-b, N))
    end
end

Base.hash(x::HeisenbergIrrep, h::UInt) = hash(x.n, h)
Base.convert(::Type{HeisenbergIrrep{N}}, n::Integer) where {N} = HeisenbergIrrep{N}(n)

Base.getindex(::IrrepTable, ::Type{Heisenberg{N}}) where {N} = HeisenbergIrrep{N}

const Heis3Irrep = HeisenbergIrrep{3}

# Sector iterator
# ---------------
function Base.isless(a::HeisenbergIrrep{N}, b::HeisenbergIrrep{N}) where {N}
    na, nb = a.n, b.n
    if na >= 0 && nb >= 0 # both characters
        return na < nb
    elseif na < 0 && nb < 0 # both Schrödinger: order on k is reversed order on n
        return na > nb
    else # characters come before Schrödinger irreps
        return na >= 0
    end
end
Base.IteratorSize(::Type{SectorValues{<:HeisenbergIrrep}}) = Base.HasLength()
Base.length(::SectorValues{HeisenbergIrrep{N}}) where {N} = N^2 + N - 1

function Base.iterate(v::SectorValues{HeisenbergIrrep{N}}, i = 1) where {N}
    return i > length(v) ? nothing : (v[i], i + 1)
end

@inline function Base.getindex(v::SectorValues{HeisenbergIrrep{N}}, i::Int) where {N}
    L = length(v)
    @boundscheck 1 <= i <= L || throw(BoundsError(v, i))
    if i <= N^2 # first N^2 should be the 1d irreps, then the N-dimensional irreps
        return HeisenbergIrrep{N}(i - 1)
    else
        return HeisenbergIrrep{N}(N^2 - i)
    end
end

findindex(::SectorValues{HeisenbergIrrep{N}}, x::HeisenbergIrrep{N}) where {N} = isschrodinger(x) ? N^2 - x.n : x.n + 1

# Product iterator
# ----------------

const HeisenbergIrrepProdIterator{N} = SectorProductIterator{HeisenbergIrrep{N}}
⊗(a::HeisenbergIrrep{N}, b::HeisenbergIrrep{N}) where {N} = SectorProductIterator((a <= b ? (a, b) : (b, a))...)

Base.IteratorSize(::Type{<:HeisenbergIrrepProdIterator}) = Base.HasLength()
function Base.length(x::HeisenbergIrrepProdIterator{N}) where {N}
    a, b = x.a, x.b
    if isschrodinger(a) && isschrodinger(b) # π ⊗ π
        return iszero(mod(a.n + b.n, N)) ? N^2 : 1 # special case: k1 + k2 = 0 gives N^2 χ's, otherwise gives N π irreps but count it once
    else
        return 1
    end
end

function Base.iterate(p::HeisenbergIrrepProdIterator{N}, state::Int = 1) where {N}
    a, b = p.a, p.b
    if state == 1
        if isschrodinger(a) && isschrodinger(b) # π ⊗ π
            k_new = mod(-(a.n + b.n), N)
            if iszero(k_new) # special case: k1 + k2 = 0 gives the N^2 χ's
                return (unit(typeof(a)), 2) # return unit χ first, then iterate over the other χ's in the next states
            else
                return (HeisenbergIrrep{N}(-k_new), 2) # return this only once, even with multiplicity N
            end
        elseif !isschrodinger(a) && !isschrodinger(b) # χ ⊗ χ
            a1, b1 = divrem(a.n, N)
            a2, b2 = divrem(b.n, N)
            return (HeisenbergIrrep{N}(mod(a1 + a2, N) * N + mod(b1 + b2, N)), 2)
        else # χ ⊗ π or π ⊗ χ
            k = isschrodinger(a) ? -a.n : -b.n
            return (HeisenbergIrrep{N}(-k), 2)
        end
    elseif state <= length(p) # π ⊗ π
        k_new = mod(-(a.n + b.n), N)
        if iszero(k_new) # special case
            return (HeisenbergIrrep{N}(state - 1), state + 1)
        else # all other π irreps
            return nothing
        end
    else
        return nothing
    end
end

# Topological data
# ----------------
dim(a::HeisenbergIrrep{N}) where {N} = isschrodinger(a) ? N : 1

function Nsymbol(a::HeisenbergIrrep{N}, b::HeisenbergIrrep{N}, c::HeisenbergIrrep{N}) where {N}
    a_1d = !isschrodinger(a)
    b_1d = !isschrodinger(b)
    c_1d = !isschrodinger(c)

    # immediate zeroes
    (a_1d && b_1d && !c_1d) && return 0 # χ ⊗ χ -> π
    (a_1d && !b_1d && c_1d) && return 0 # χ ⊗ π -> χ
    (!a_1d && b_1d && c_1d) && return 0 # π ⊗ χ -> χ

    if a_1d && b_1d && c_1d # χ ⊗ χ -> χ
        a1, b1 = divrem(a.n, N)
        a2, b2 = divrem(b.n, N)
        a3, b3 = divrem(c.n, N)
        return Int(mod(a1 + a2, N) == a3 && mod(b1 + b2, N) == b3)
    elseif ((a_1d && !b_1d) || (!a_1d && b_1d)) && !c_1d # χ ⊗ π or π ⊗ χ -> π
        k_ab = a_1d ? -b.n : -a.n
        return Int(k_ab == -c.n)
    else # π ⊗ π
        k_new = mod(-(a.n + b.n), N)
        iszero(k_new) && return Int(c_1d) # special case
        return k_new == -c.n ? N : 0
    end
end

Fsymbol(a::I, b::I, c::I, d::I, e::I, f::I) where {I <: HeisenbergIrrep} =
    Fsymbol_from_fusiontensor(a, b, c, d, e, f)

Rsymbol(a::I, b::I, c::I) where {I <: HeisenbergIrrep} = Rsymbol_from_fusiontensor(a, b, c)

# https://www.rintonpress.com/xxqic8/qic-8-5/0438-0467.pdf eqs 60, 63, 67
function fusiontensor(x::HeisenbergIrrep{N}, y::HeisenbergIrrep{N}, z::HeisenbergIrrep{N}) where {N}
    T = fusionscalartype(HeisenbergIrrep{N})
    Nxyz = Nsymbol(x, y, z)
    dx, dy, dz = dim(x), dim(y), dim(z)
    C = zeros(T, dx, dy, dz, Nxyz)
    isempty(C) && return C

    if !isschrodinger(x) && !isschrodinger(y) # χ ⊗ χ → χ (trivial)
        C[1, 1, 1, 1] = one(T)
        return C
    end

    ω = cispi(2 / N)
    if isschrodinger(x) && isschrodinger(y) # π_k ⊗ π_{k'}
        invsqrtN = T(1 / sqrt(N))
        kx, ky = -x.n, -y.n
        K = mod(kx + ky, N)
        if !iszero(K) # → π_K with multiplicity N
            s = mod(-kx * invmod(ky, N), N)
            for i in 1:dx, j in 1:dy, m in 1:dz
                if iszero(mod(j - m - s * (i - m), N)) # j - m = s(i - m) mod N, s = -k / k' mod N
                    μ = mod(i - j, N) + 1 # C[i,j,m,μ] = 1 if μ - 1 = i - j mod N (pure permutation)
                    C[i, j, m, μ] = one(T)
                end
            end
        else # k + k' = 0: → Σ_{a,b} χ_{a,b}, each with multiplicity 1
            a, b = divrem(z.n, N)
            t = mod(invmod(kx, N) * b, N)
            for i in 1:dx, j in 1:dy
                if iszero(mod(j - i + t, N)) # C[i,j,1,1] = ω^{-a(i-1)} / √N if j = i - t mod N, t = b/k mod N
                    C[i, j, 1, 1] = ω^mod(-a * (i - 1), N) * invsqrtN
                end
            end
        end
    else # exactly one of x, y is χ, the other is π
        if !isschrodinger(x) # χ_{a,b} ⊗ π_k → π_k
            a, b = divrem(x.n, N)
            k = -y.n
            s = mod(invmod(k, N) * b, N)
            for j in 1:dy, m in 1:dz
                if iszero(mod(j - m + s, N)) # C[1,j,m,1] = ω^{a(m-1)} if j = m - s mod N, s = b/k mod N
                    C[1, j, m, 1] = ω^mod(a * (m - 1), N)
                end
            end
        else # π_k ⊗ χ_{a,b} → π_k
            a, b = divrem(y.n, N)
            k = -x.n
            s = mod(invmod(k, N) * b, N)
            for i in 1:dx, m in 1:dz
                if iszero(mod(i - m + s, N)) # C[i,1,m,1] = ω^{a(m-1)} if i = m - s mod N, s = b/k mod N
                    C[i, 1, m, 1] = ω^mod(a * (m - 1), N)
                end
            end
        end
    end
    return C
end
