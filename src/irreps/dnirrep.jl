"""
    struct DNIrrep{N} <: AbstractIrrep{Dihedral{N}}
    DNIrrep{N}(n::Integer, isodd::Bool=false)
    Irrep[Dihedral{N}](n::Integer, isodd::Bool=false)

Represents irreps of the dihedral group ``D_N = Z_N ⋊ C`` (``Z_N`` and charge conjugation or reflection).

## Properties

- `j::Int`: the value of the ``Z_N`` charge.
- `isodd::Bool`: the representation of charge conjugation.

Combined these take the values ``+0, -0, 1, ..., (N - 1) / 2`` for odd ``N``, and
``+0, -0, 1, ..., N / 2 - 1, +(N/2), -(N/2)`` for even ``N``, where the ``+`` (``-``)
refer to the even (odd) one-dimensional irreps, while the others are two-dimensional.
"""
struct DNIrrep{N} <: AbstractIrrep{Dihedral{N}}
    # store isodd in right bit, use rest for storing the angle
    data::UInt8
    DNIrrep{N}(data::UInt8) where {N} = new{N}(data)
end

function DNIrrep{N}(j::Integer, isodd::Bool = false) where {N}
    @assert 0 < N < 32 "DNIrrep requires 0 < N < 32 to function properly"
    0 <= j <= (N >> 1) ||
        throw(DomainError(j, "DNIrrep only has irreps 0 <= j <= (N ÷ 2)"))
    !isodd || j == 0 || (iseven(N) && (j == (N >> 1))) ||
        throw(DomainError(j, "DNIrrep only has odd irreps when `j == 0` or `iseven(N) && j == N / 2`"))
    return DNIrrep{N}((j % UInt8) << 1 | isodd)
end

Base.propertynames(x::DNIrrep) = (:j, :isodd, :data)

function Base.getproperty(a::DNIrrep{N}, sym::Symbol) where {N}
    if sym === :j
        return getfield(a, :data) >> 1
    elseif sym === :isodd
        return Bool(getfield(a, :data) & true)
    elseif sym === :data
        return getfield(a, :data)
    else
        error("Unknown property $sym")
    end
end

FusionStyle(::Type{DNIrrep{N}}) where {N} = N < 3 ? UniqueFusion() : SimpleFusion()
fusionscalartype(::Type{DNIrrep{N}}) where {N} = Float64
braidingscalartype(::Type{DNIrrep{N}}) where {N} = Float64
sectorscalartype(::Type{DNIrrep{N}}) where {N} = Float64

unit(::Type{DNIrrep{N}}) where {N} = DNIrrep{N}(0, false)
dual(a::DNIrrep) = a

Base.hash(a::DNIrrep, h::UInt) = hash(a.data, h)
Base.convert(::Type{DNIrrep{N}}, (j, n)::Tuple{Integer, Bool}) where {N} = DNIrrep{N}(j, n)

Base.getindex(::IrrepTable, ::Type{Dihedral{N}}) where {N} = DNIrrep{N}

const D3Irrep = DNIrrep{3}
const D4Irrep = DNIrrep{4}

function Base.show(io::IO, a::DNIrrep)
    if get(io, :typeinfo, nothing) !== typeof(a)
        print(io, type_repr(typeof(a)))
    end
    print(io, "(", a.j, ", ", a.isodd, ")")
    return nothing
end

# Sector iterator
# ---------------
Base.isless(a::DNIrrep{N}, b::DNIrrep{N}) where {N} = isless(a.data, b.data)
Base.IteratorSize(::Type{SectorValues{DNIrrep{N}}}) where {N} = Base.HasLength()
Base.length(::SectorValues{DNIrrep{N}}) where {N} = (N >> 1) + (isodd(N) ? 2 : 3)

function Base.iterate(v::SectorValues{DNIrrep{N}}, i = 1) where {N}
    return i > length(v) ? nothing : (v[i], i + 1)
end

@inline function Base.getindex(v::SectorValues{DNIrrep{N}}, i::Int) where {N}
    L = length(v)
    @boundscheck 1 <= i <= L || throw(BoundsError(v, i))
    if i == 1
        return DNIrrep{N}(0x00)
    elseif i == 2
        return DNIrrep{N}(0x01)
    elseif iseven(N) && i == L
        return DNIrrep{N}((UInt8(N >> 1) << 1) | true)
    else
        return DNIrrep{N}(UInt8(i - 2) << 1)
    end
end

function findindex(::SectorValues{DNIrrep{N}}, a::DNIrrep{N}) where {N}
    if a.isodd && a.j > 0
        return Int(a.j) + 3
    elseif a.data < 2
        return Int(a.data) + 1
    else
        return Int(a.j) + 2
    end
end

# Product iterator
# ----------------

const DNIrrepProdIterator{N} = SectorProductIterator{DNIrrep{N}}
⊗(a::DNIrrep{N}, b::DNIrrep{N}) where {N} = SectorProductIterator((a <= b ? (a, b) : (b, a))...)

Base.IteratorSize(::Type{DNIrrepProdIterator{N}}) where {N} = Base.HasLength()
function Base.length(x::DNIrrepProdIterator{N}) where {N}
    a, b = x.a, x.b
    if a.j == 0 || b.j == 0 || (N == 2 * a.j) || (N == 2 * b.j)
        return 1
    end

    Asplits = iszero(b.j - a.j)
    Bsplits = 2 * (a.j + b.j) == N
    return 2 + Asplits + Bsplits
end

function Base.iterate(p::DNIrrepProdIterator{N}, state::Int = 1) where {N}
    a, b = p.a, p.b
    if state == 1
        # A_i x A_j = A_{xor(i,j)}
        if iszero(b.j) # then also iszero(a.j)
            cj = 0
            cisodd = xor(a.isodd, b.isodd)
            return DNIrrep{N}(cj, cisodd), 5
        end

        if iszero(a.j)
            cj = b.j

            # A_i x B_j = B_{xor(i,j)}
            if 2 * cj == N
                cisodd = xor(a.isodd, b.isodd)
                return DNIrrep{N}(cj, cisodd), 5
            end

            # A_i x rho_j = rho_j
            return b, 5
        end

        if a.j == b.j # != 0
            # B_i x B_j = A_{xor(i,j)}
            if 2 * b.j == N
                return DNIrrep{N}(0, xor(a.isodd, b.isodd)), 5
            end

            # rho_i x rho_i = A_1 + A_2 + rho_2i
            return DNIrrep{N}(0, false), 2
        end

        # rho_i x B_j = rho_{j-i}
        if 2 * b.j == N
            return DNIrrep{N}(b.j - a.j, false), 5
        end

        # rho_i x rho_j = rho_{i-j} + rho_{i+j}
        return DNIrrep{N}(b.j - a.j, false), 3
    elseif state == 2
        return DNIrrep{N}(0, true), 3
    elseif state == 3
        cj = a.j + b.j
        if 2 * cj == N
            return DNIrrep{N}(cj, false), 4
        end
        if cj > (N >> 1)
            cj = N - cj
        end
        return DNIrrep{N}(cj, false), 5
    elseif state == 4
        return DNIrrep{N}(N >> 1, true), 5
    else
        return nothing
    end
end

# Topological data
# ----------------
dim(a::DNIrrep{N}) where {N} = ifelse((a.j == 0) | (2 * a.j == N), 1, 2)

function Nsymbol(a::DNIrrep{N}, b::DNIrrep{N}, c::DNIrrep{N}) where {N}
    j_satisfied = (c.j == min(a.j + b.j, N - (a.j + b.j))) | (c.j == max(a.j, b.j) - min(a.j, b.j))
    s_satisfied = (dim(c) != 1) | (dim(a) != 1) | (xor(a.isodd, b.isodd) == c.isodd)
    return j_satisfied & s_satisfied
end

function Fsymbol(a::I, b::I, c::I, d::I, e::I, f::I) where {N, I <: DNIrrep{N}}
    T = fusionscalartype(I)
    (Nsymbol(a, b, e) & Nsymbol(e, c, d) & Nsymbol(b, c, f) & Nsymbol(a, f, d)) || return zero(T)

    # tensoring with units gives 1
    (isunit(a) || isunit(b) || isunit(c)) && return one(T)

    # fallback through fusiontensor
    A = fusiontensor(a, b, e)
    A = reshape(A, size(A, 1), size(A, 2), size(A, 3))
    B = fusiontensor(e, c, d)
    B1 = reshape(view(B, :, :, 1), size(B, 1), size(B, 2))
    C = fusiontensor(b, c, f)
    C = reshape(C, size(C, 1), size(C, 2), size(C, 3))
    D = fusiontensor(a, f, d)
    D1 = reshape(view(D, :, :, 1), size(D, 1), size(D, 2))

    return @tensor conj(D1[1 5]) * conj(C[2 4 5]) * A[1 2 3] * B1[3 4]
end

function Rsymbol(a::I, b::I, c::I) where {N, I <: DNIrrep{N}}
    R = convert(braidingscalartype(I), Nsymbol(a, b, c))
    return ifelse((c.j == 0) & c.isodd & !(a.j == b.j == 0) & !((2 * a.j) == (2 * b.j) == N), -R, R)
end


function fusiontensor(a::I, b::I, c::I) where {N, I <: DNIrrep{N}}
    T = sectorscalartype(I)
    C = zeros(T, dim(a), dim(b), dim(c), 1)
    Nsymbol(a, b, c) || return C

    if c.j == 0
        if a.j == b.j == 0 || (2 * a.j == 2 * b.j == N)
            C[1, 1, 1] = 1
        else # a.j == b.j
            # 0\pm = 1/sqrt(2) (v_i^+ \otimes w_j^- \pm v_i^- \otimes w_j^+)
            C[1, 2, 1] = T(sqrt(2) / 2)
            C[2, 1, 1] = c.isodd ? -C[1, 2, 1] : C[1, 2, 1]
        end
    elseif 2 * c.j == N
        if (a.j == (N >> 1)) | (b.j == (N >> 1))
            C[1, 1, 1] = 1
        else
            C[1, 1, 1] = T(sqrt(2) / 2)
            C[2, 2, 1] = c.isodd ? -C[1, 1, 1] : C[1, 1, 1]
        end
    elseif a.j == 0
        C[1, 1, 1] = 1
        C[1, 2, 2] = a.isodd ? -1 : 1
    elseif b.j == 0
        C[1, 1, 1] = 1
        C[2, 1, 2] = b.isodd ? -1 : 1
    elseif 2 * a.j == N
        C[1, 1, 2] = 1
        C[1, 2, 1] = a.isodd ? -1 : 1
    elseif 2 * b.j == N
        C[1, 1, 2] = 1
        C[2, 1, 1] = b.isodd ? -1 : 1
        # from here on everything is 2D --------------------
    elseif c.j == a.j + b.j
        C[1, 1, 1] = 1
        C[2, 2, 2] = 1
    elseif c.j == N - (a.j + b.j)
        C[1, 1, 2] = 1
        C[2, 2, 1] = 1
    elseif c.j == a.j - b.j
        C[1, 2, 1] = 1
        C[2, 1, 2] = 1
    elseif c.j == b.j - a.j
        C[1, 2, 2] = 1
        C[2, 1, 1] = 1
    end

    return C
end
