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
    n::UInt8
    ANIrrep{N}(data::UInt8) where {N} = new{N}(n)
end

function ANIrrep{N}(j::Integer, isodd::Bool = false) where {N}
    @assert 0 < N < 5 "ANIrrep is only provided for 0 < N < 5"
    0 <= j <= _npartitions(N) ||
        throw(DomainError(j, "ANIrrep only has irreps 0 <= j <= (N ÷ 2)"))
    return ANIrrep{N}(j)
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
sectorscalartype(::Type{ANIrrep{N}}) where {N} = Float64
Base.isreal(::Type{ANIrrep{N}}) where {N} = true

unit(::Type{ANIrrep{N}}) where {N} = ANIrrep{N}(0)
function dual(a::ANIrrep{N}) where {N}
    N < 3 && return a
    if N == 3
        return a.n == 2 ? ANIrrep{3}(1) : a
    else # N = 4
        return a.n == 3 ? ANIrrep{4}(3) : (a.n == 2 ? ANIrrep{4}(1) : a)
    end
end

Base.hash(a::ANIrrep, h::UInt) = hash(a.n, h)
Base.convert(::Type{ANIrrep{N}}, n::Integer) where {N} = ANIrrep{N}(n)

Base.getindex(::IrrepTable, ::Type{Alternating{N}}) where {N} = ANIrrep{N}

const A2Irrep = ANIrrep{2} # trivial
const A3Irrep = ANIrrep{3} # Z3
const A4Irrep = ANIrrep{4}

function Base.show(io::IO, a::ANIrrep)
    if get(io, :typeinfo, nothing) !== typeof(a)
        print(io, type_repr(typeof(a)))
    end
    print(io, "(", a.n, ")")
    return nothing
end

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
        return a == 3 ? 4 : 1
    else
        return 1
    end
end
# continue here
function Base.iterate(p::ANIrrepProdIterator{N}, state::Int = 1) where {N}
    a, b = p.a, p.b
    if state == 1
        # A_i x A_j = A_{xor(i,j)}
        if iszero(b.j) # then also iszero(a.j)
            cj = 0
            cisodd = xor(a.isodd, b.isodd)
            return ANIrrep{N}(cj, cisodd), 5
        end

        if iszero(a.j)
            cj = b.j

            # A_i x B_j = B_{xor(i,j)}
            if 2 * cj == N
                cisodd = xor(a.isodd, b.isodd)
                return ANIrrep{N}(cj, cisodd), 5
            end

            # A_i x rho_j = rho_j
            return b, 5
        end

        if a.j == b.j # != 0
            # B_i x B_j = A_{xor(i,j)}
            if 2 * b.j == N
                return ANIrrep{N}(0, xor(a.isodd, b.isodd)), 5
            end

            # rho_i x rho_i = A_1 + A_2 + rho_2i
            return ANIrrep{N}(0, false), 2
        end

        # rho_i x B_j = rho_{j-i}
        if 2 * b.j == N
            return ANIrrep{N}(b.j - a.j, false), 5
        end

        # rho_i x rho_j = rho_{i-j} + rho_{i+j}
        return ANIrrep{N}(b.j - a.j, false), 3
    elseif state == 2
        return ANIrrep{N}(0, true), 3
    elseif state == 3
        cj = a.j + b.j
        if 2 * cj == N
            return ANIrrep{N}(cj, false), 4
        end
        if cj > (N >> 1)
            cj = N - cj
        end
        return ANIrrep{N}(cj, false), 5
    elseif state == 4
        return ANIrrep{N}(N >> 1, true), 5
    else
        return nothing
    end
end

# Topological data
# ----------------
dim(a::ANIrrep{N}) where {N} = ifelse((a.j == 0) | (2 * a.j == N), 1, 2)

function Nsymbol(a::ANIrrep{N}, b::ANIrrep{N}, c::ANIrrep{N}) where {N}
end

function Fsymbol(a::I, b::I, c::I, d::I, e::I, f::I) where {N, I <: ANIrrep{N}}
    T = sectorscalartype(I)
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

function Rsymbol(a::I, b::I, c::I) where {N, I <: ANIrrep{N}}
    R = convert(sectorscalartype(I), Nsymbol(a, b, c))
    return ifelse((c.j == 0) & c.isodd & !(a.j == b.j == 0) & !((2 * a.j) == (2 * b.j) == N), -R, R)
end

# Fusion tensors for Rep(A4)
# dims: 1 => 1, 1p => 1, 1pp => 1, 3 => 3
# Tensors have shape (dim_out, dim_inA, dim_inB)
#
# Conventions:
# - Triplet basis components ordered (1,2,3).
# - omega = exp(2pi*im/3).
# - Normalization: singlets 1/sqrt(3), symmetric/antisymmetric triplet components 1/sqrt(2).

const ω = cis(2π/3)   # e^{2π i / 3}

# 3 ⊗ 3 -> 1  (scalar)
T_3x3_to_1 = zeros(ComplexF64, 1, 3, 3)
for i in 1:3
    T_3x3_to_1[1,i,i] = 1/sqrt(3)
end

# 3 ⊗ 3 -> 1' 
T_3x3_to_1p = zeros(ComplexF64, 1, 3, 3)
T_3x3_to_1p[1,1,1] = 1/sqrt(3)
T_3x3_to_1p[1,2,2] = ω/sqrt(3)
T_3x3_to_1p[1,3,3] = ω^2/sqrt(3)

# 3 ⊗ 3 -> 1'' 
T_3x3_to_1pp = zeros(ComplexF64, 1, 3, 3)
T_3x3_to_1pp[1,1,1] = 1/sqrt(3)
T_3x3_to_1pp[1,2,2] = ω^2/sqrt(3)
T_3x3_to_1pp[1,3,3] = ω/sqrt(3)

# 3 ⊗ 3 -> 3_Symmetric  (symmetric triplet)
T_3x3_to_3S = zeros(ComplexF64, 3, 3, 3)
# components: (a2 b3 + a3 b2, a3 b1 + a1 b3, a1 b2 + a2 b1)
T_3x3_to_3S[1,2,3] = 1/sqrt(2); T_3x3_to_3S[1,3,2] = 1/sqrt(2)
T_3x3_to_3S[2,3,1] = 1/sqrt(2); T_3x3_to_3S[2,1,3] = 1/sqrt(2)
T_3x3_to_3S[3,1,2] = 1/sqrt(2); T_3x3_to_3S[3,2,1] = 1/sqrt(2)

# 3 ⊗ 3 -> 3_Antisymmetric  (antisymmetric triplet)
T_3x3_to_3A = zeros(ComplexF64, 3, 3, 3)
# components: (a2 b3 - a3 b2, a3 b1 - a1 b3, a1 b2 - a2 b1)
T_3x3_to_3A[1,2,3] =  1/sqrt(2); T_3x3_to_3A[1,3,2] = -1/sqrt(2)
T_3x3_to_3A[2,3,1] =  1/sqrt(2); T_3x3_to_3A[2,1,3] = -1/sqrt(2)
T_3x3_to_3A[3,1,2] =  1/sqrt(2); T_3x3_to_3A[3,2,1] = -1/sqrt(2)

# 1 × anything and anything × 1  (trivial fusion: scalar multiplies vector)
# Represented as identity maps: out_dim = inB_dim (or inA_dim)
# 1 ⊗ 3 -> 3
T_1x3_to_3 = zeros(ComplexF64, 3, 1, 3)
for i in 1:3
    T_1x3_to_3[i,1,i] = 1.0
end
# 3 ⊗ 1 -> 3
T_3x1_to_3 = zeros(ComplexF64, 3, 3, 1)
for i in 1:3
    T_3x1_to_3[i,i,1] = 1.0
end

# 1' ⊗ 3 -> 3  and 3 ⊗ 1' -> 3  (scalar multiplication map)
# (the representation matrix of 1' enters in the action; the fusion tensor itself is the trivial identification)
T_1p_x_3_to_3 = zeros(ComplexF64, 3, 1, 3)
for i in 1:3
    T_1p_x_3_to_3[i,1,i] = 1.0
end
T_3_x_1p_to_3 = zeros(ComplexF64, 3, 3, 1)
for i in 1:3
    T_3_x_1p_to_3[i,i,1] = 1.0
end

# 1' ⊗ 1' -> 1'' ; 1' ⊗ 1'' -> 1 ; 1'' ⊗ 1'' -> 1'
T_1p_x_1p_to_1pp = zeros(ComplexF64, 1, 1, 1); T_1p_x_1p_to_1pp[1,1,1] = 1.0
T_1p_x_1pp_to_1    = zeros(ComplexF64, 1, 1, 1); T_1p_x_1pp_to_1[1,1,1] = 1.0
T_1pp_x_1pp_to_1p  = zeros(ComplexF64, 1, 1, 1); T_1pp_x_1pp_to_1p[1,1,1] = 1.0

function fusiontensor(a::I, b::I, c::I) where {N, I <: ANIrrep{N}}
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
