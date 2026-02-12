"""
    struct HeisenbergIrrep{N} <: AbstractIrrep{Heisenberg{N}}
    HeisenbergIrrep{N}(n::Integer, isodd::Bool=false)
    Irrep[Heisenberg{N}](n::Integer, isodd::Bool=false)

Represents irreps of the finite Heisenberg group ``H_N``.
There are N² one-dimensional irreps labeled by a pair of integers ``(a, b)`` with ``a, b ∈ ℤ_N``,
and N - 1 irreps of dimension N labeled by an integer ``k ∈ ℤ_N \\ {0}``.
These are commonly referred to as the characters ``χₐ,ᵦ`` and the Schrödinger representations ``πₖ``, respectively.

## Properties

- `a::Int`: the first label for the one-dimensional irreps.
- `b::Int`: the second label for the one-dimensional irreps.
- `k::Int`: the label for the N-dimensional irreps.

The last index `k` is only relevant for the N-dimensional irreps, and is set to 0 for the one-dimensional irreps.
Similarly, the first two indices `a` and `b` are only relevant for the one-dimensional irreps, and are set to 0 for the N-dimensional irreps.
"""
struct HeisenbergIrrep{N} <: AbstractIrrep{Heisenberg{N}}
    a::Int
    b::Int
    k::Int
    function HeisenbergIrrep{N}(a::Int, b::Int, k::Int) where {N}
        if iszero(k) # 1d irreps
            a_mod = mod(a, N)
            b_mod = mod(b, N)
            k_mod = 0
        else # Schrödinger irreps
            iszero(a) && iszero(b) || throw(ArgumentError("N-dimensional irreps must have a = b = 0"))
            a_mod = 0
            b_mod = 0
            k_mod = mod(k, N)
        end
        return new{N}(a_mod, b_mod, k_mod)
    end
end
# TODO: consider union of 2 kinds of irreps?

FusionStyle(::Type{HeisenbergIrrep{N}}) where {N} = GenericFusion()
BraidingStyle(::Type{HeisenbergIrrep{N}}) where {N} = NoBraiding() #TODO: make braided again after fusiontensor is fixed
fusionscalartype(::Type{HeisenbergIrrep{N}}) where {N} = ComplexF64
sectorscalartype(::Type{HeisenbergIrrep{N}}) where {N} = ComplexF64

unit(::Type{HeisenbergIrrep{N}}) where {N} = HeisenbergIrrep{N}(0, 0, 0)
function dual(a::HeisenbergIrrep{N}) where {N}
    return if iszero(a.k)
        HeisenbergIrrep{N}(N - a.a, N - a.b, 0)
    else
        HeisenbergIrrep{N}(0, 0, N - a.k)
    end
end

Base.hash(a::HeisenbergIrrep, h::UInt) = hash((a.a, a.b, a.k), h)
Base.convert(::Type{HeisenbergIrrep{N}}, (a, b, k)::NTuple{3, Int}) where {N} = HeisenbergIrrep{N}(a, b, k)

Base.getindex(::IrrepTable, ::Type{Heisenberg{N}}) where {N} = HeisenbergIrrep{N}

const Heis3Irrep = HeisenbergIrrep{3}

function Base.show(io::IO, a::HeisenbergIrrep{N}) where {N}
    if get(io, :typeinfo, nothing) !== typeof(a)
        print(io, type_repr(typeof(a)))
    end
    print(io, "(", a.a, ", ", a.b, ", ", a.k, ")")
    return nothing
end

# Sector iterator
# ---------------
function Base.isless(a::HeisenbergIrrep{N}, b::HeisenbergIrrep{N}) where {N}
    if iszero(a.k)
        return iszero(b.k) ? isless((a.a, a.b), (b.a, b.b)) : true # 1d irreps come before Schrödinger irreps
    else
        return iszero(b.k) ? false : isless(a.k, b.k)
    end
end
Base.IteratorSize(::Type{SectorValues{HeisenbergIrrep{N}}}) where {N} = Base.HasLength()
Base.length(::SectorValues{HeisenbergIrrep{N}}) where {N} = N^2 + N - 1

function Base.iterate(v::SectorValues{HeisenbergIrrep{N}}, i = 1) where {N}
    return i > length(v) ? nothing : (v[i], i + 1)
end

@inline function Base.getindex(v::SectorValues{HeisenbergIrrep{N}}, i::Int) where {N}
    L = length(v)
    @boundscheck 1 <= i <= L || throw(BoundsError(v, i))
    if i <= N^2 # first N^2 should be the 1d irreps, then the N-dimensional irreps
        a = div(i - 1, N) # for fixed a provide all b's, then move to next a
        b = mod(i - 1, N)
        return HeisenbergIrrep{N}(a, b, 0)
    else
        return HeisenbergIrrep{N}(0, 0, i - N^2)
    end
end

function findindex(::SectorValues{HeisenbergIrrep{N}}, a::HeisenbergIrrep{N}) where {N}
    if iszero(a.k) # 1d irreps
        return a.a * N + a.b + 1
    else
        return N^2 + a.k
    end
end

# Product iterator
# ----------------

const HeisenbergIrrepProdIterator{N} = SectorProductIterator{HeisenbergIrrep{N}}
⊗(a::HeisenbergIrrep{N}, b::HeisenbergIrrep{N}) where {N} = SectorProductIterator((a <= b ? (a, b) : (b, a))...)

Base.IteratorSize(::Type{HeisenbergIrrepProdIterator{N}}) where {N} = Base.HasLength()
function Base.length(x::HeisenbergIrrepProdIterator{N}) where {N}
    a, b = x.a, x.b
    if !iszero(a.k) && !iszero(b.k) # π ⊗ π
        return iszero(mod(a.k + b.k, N)) ? N^2 : N # special case: k1 + k2 = 0 gives N^2 χ's, otherwise gives N π irreps
    else
        return 1
    end
end

function Base.iterate(p::HeisenbergIrrepProdIterator{N}, state::Int = 1) where {N}
    a, b = p.a, p.b
    if state == 1
        if !(iszero(a.k) && iszero(b.k)) # π ⊗ π
            k_new = mod(a.k + b.k, N)
            if iszero(k_new) # special case: k1 + k2 = 0 gives the N^2 χ's
                return (unit(typeof(a)), 2) # return unit χ first, then iterate over the other χ's in the next states
            else
                return (HeisenbergIrrep{N}(0, 0, k_new), 2) # return this only once, even with multiplicity N
            end
        elseif iszero(a.k) && iszero(b.k) # χ ⊗ χ
            a_new = mod(a.a + b.a, N)
            b_new = mod(a.b + b.b, N)
            return (HeisenbergIrrep{N}(a_new, b_new, 0), 2)
        else # χ ⊗ π or π ⊗ χ
            return (HeisenbergIrrep{N}(0, 0, a.k + b.k), 2)
        end
    elseif state <= length(p) # π ⊗ π
        k_new = mod(a.k + b.k, N)
        if iszero(k_new) # special case
            return (HeisenbergIrrep{N}(div(state - 1, N), mod(state - 1, N), 0), state + 1)
        else # all other π irreps
            return nothing
        end
    else
        return nothing
    end
end

# Topological data
# ----------------
dim(a::HeisenbergIrrep{N}) where {N} = iszero(a.k) ? 1 : N

#TODO: k1 + k2 = 0 is a special case
function Nsymbol(a::HeisenbergIrrep{N}, b::HeisenbergIrrep{N}, c::HeisenbergIrrep{N}) where {N}
    a_1d = iszero(a.k)
    b_1d = iszero(b.k)
    c_1d = iszero(c.k)

    # immediate zeroes
    (a_1d && b_1d && !c_1d) && return 0 # χ ⊗ χ -> π
    (a_1d && !b_1d && c_1d) && return 0 # χ ⊗ π -> χ
    (!a_1d && b_1d && c_1d) && return 0 # π ⊗ χ -> χ

    if a_1d && b_1d && c_1d # χ ⊗ χ -> χ
        return Int(mod(a.a + b.a, N) == c.a && mod(a.b + b.b, N) == c.b)
    elseif ((a_1d && !b_1d) || (!a_1d && b_1d)) && !c_1d # χ ⊗ π or π ⊗ χ -> π
        return Int(a.k + b.k == c.k)
    else # π ⊗ π
        k_new = mod(a.k + b.k, N)
        iszero(k_new) && return Int(c_1d) # special case
        return mod(k_new, N) == c.k ? N : 0
    end
end

# function Fsymbol(
#     a::HeisenbergIrrep{N}, b::HeisenbergIrrep{N}, c::HeisenbergIrrep{N},
#     d::HeisenbergIrrep{N}, e::HeisenbergIrrep{N}, f::HeisenbergIrrep{N}
#     ) where {N}
#     T = fusionscalartype(HeisenbergIrrep{N})
#     Nabe, Necd, Nbcf, Nafd = Nsymbol
#     iszero(a.k) && iszero(b.k) && iszero(c.k) && return ones(T, Nabe, Necd, Nbcf, Nafd) # 1d ⊗ 1d ⊗ 1d -> 1d is trivial

#     # most nontrivial case: Schrödinger ⊗ Schröd
#     iszero(a.k) && iszero(b.k) && iszero(c.k) && return ones(T, Nabe, Necd, Nbcf, Nafd) # 1d ⊗ 1d ⊗ 1d -> 1d is trivial

#     # most nontrivial case: Schrödinger ⊗ Schrödinger ⊗ Schrödinger
#     if !iszero(a.k) && !iszero(b.k) && !iszero(c.k)
#         if !iszero(d.k) # fuse to Schrödinger
#             # TODO: don't be lazy and figure this out
#         end
#     end
#     return F
# end

Fsymbol(a::I, b::I, c::I, d::I, e::I, f::I) where {I <: HeisenbergIrrep} =
    Fsymbol_from_fusiontensor(a, b, c, d, e, f)

# Rsymbol(a::I, b::I, c::I) where {I <: HeisenbergIrrep} = Rsymbol_from_fusiontensor(a, b, c)

# https://www.rintonpress.com/xxqic8/qic-8-5/0438-0467.pdf eqs 60, 63, 67
function fusiontensor(x::HeisenbergIrrep{N}, y::HeisenbergIrrep{N}, z::HeisenbergIrrep{N}) where {N}
    T = fusionscalartype(HeisenbergIrrep{N})
    Nxyz = Nsymbol(x, y, z)
    dx, dy, dz = dim(x), dim(y), dim(z)
    C = zeros(T, dx, dy, dz, Nxyz)
    isempty(C) && return C

    if iszero(x.k) && iszero(y.k) && iszero(z.k) # χ ⊗ χ -> χ
        C[1, 1, 1, 1] = 1
        return C
    end

    ω = cispi(2 / N)
    if !iszero(x.k) && !iszero(y.k) # π ⊗ π
        k_new = mod(x.k + y.k, N)
        if !iszero(k_new) # fuse to π
            for i in 1:dx, j in 1:dy, m in 1:dz, μ in 1:Nxyz
                C[i, j, m, μ] = T(μ == i) * T(m == (j - m * x.k / k_new)) # missing elements on diagonal in overlap
            end
        else # special case
            for i in 1:dx, j in 1:dy
                C[i, j, 1, 1] = ω^(x.k * z.a * i) * T((j - i) == 1) / sqrt(N - 1) # right for N = 3 at least at level of orthogonality, rest untested
            end
        end
    else
        if iszero(x.k) && !iszero(y.k) # χ ⊗ π
            @info "χ x π: x = $x, y = $y"
            for j in 1:N, m in 1:N
                C[1, j, m, 1] = ω^(- x.a * j) * T(j + x.b/y.k == m) # missing elements on diagonal in overlap
            end
        else # π ⊗ χ
                @info "π x χ: x = $x, y = $y"
            for i in 1:N, m in 1:N
                C[i, 1, m, 1] = ω^(- y.a * i) * T(i + y.b/x.k == m) # missing elements on diagonal in overlap
            end
        end
    end
    return C
end
