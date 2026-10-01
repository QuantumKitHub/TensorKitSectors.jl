# Tambara-Yamagami category for the cyclic group ℤ_N
#---------------------------------------------------------------------------------------#
"""
    struct TambaraYamagami{N, K} <: Sector
    TambaraYamagami{N, K}(n::Integer)

Represents the Tambara-Yamagami fusion category built from the cyclic group ``ℤ_N``.
The simple objects are the group elements `0, 1, …, N - 1` of ``ℤ_N``,
together with a single non-invertible object `m`.

The non-trivial fusion rules are given by
```math
g ⊗ h = g + h \\mod N, \\qquad g ⊗ m = m ⊗ g = m, \\qquad m ⊗ m = \\bigoplus_{g ∈ ℤ_N} g.
```

The F-symbols are constructed from a non-degenerate symmetric bicharacter
``χ(g, h) = \\exp(2π i g h / N)``, together with a Frobenius-Schur sign for the
non-invertible object; the latter is fixed as `κ = (-1)^K` through the type parameter
`K::Bool`. For fixed `N`, the two choices `K = false` (`κ = 1`) and `K = true` (`κ = -1`)
generally give distinct fusion categories.

Only the case `N == 2` and `K == false` admits a braiding, as this case coincides with Ising, but this is not currently implemented.

## Fields
- `n::UInt8`: a group element for `0 <= n < N`, or the non-invertible object `m` for `n == N`.

## References
[1] D. Tambara and S. Yamagami, *Tensor categories with fusion rules of self-duality for
    finite abelian groups*, J. Algebra **209**, 692-707 (1998).
[2] M. Barkeshli, P. Bonderson, M. Cheng and Z. Wang, *Symmetry Fractionalization, Defects,
    and Gauging of Topological Phases*, Phys. Rev. B **100**, 115147 (2019),
    [arXiv:1410.4540](https://arxiv.org/abs/1410.4540).
"""
struct TambaraYamagami{N, K} <: Sector
    n::UInt8
    function TambaraYamagami{N, K}(n) where {N, K}
        _check_TY_typeparams(N, K)
        0 <= n <= N || throw(DomainError(n, "TambaraYamagami{$N, $K} labels must satisfy 0 <= n <= $N"))
        return new{N, K}(n)
    end
end
function TambaraYamagami{N, K}(s::Symbol) where {N, K}
    s === :m || throw(ArgumentError("Unknown label $s: use an integer or `:m`"))
    return TambaraYamagami{N, K}(N)
end

# Labels are stored as `UInt8`; restricting to `N <= 128` guarantees that the sum `a.n + b.n`
# of two group labels (at most `2(N - 1) = 254`) never overflows before taking it modulo `N`.
const SMALL_TY_CUTOFF = (typemax(UInt8) + 1) ÷ 2

function _check_TY_typeparams(N, K)
    N isa Int && 1 <= N <= SMALL_TY_CUTOFF || throw(ArgumentError("N must be an Int satisfying 1 <= N <= $SMALL_TY_CUTOFF, got $N::$(typeof(N))"))
    K isa Bool || throw(ArgumentError("K must be a Bool, encoding the Frobenius-Schur indicator (-1)^K, got $K"))
    return nothing
end

"""
    modulus(n::TambaraYamagami{N, K}) -> N
    modulus(::Type{<:TambaraYamagami{N, K}}) -> N

The order of the cyclic group, or the modulus of the charge labels.
"""
modulus(n::TambaraYamagami) = modulus(typeof(n))
modulus(::Type{<:TambaraYamagami{N, K}}) where {N, K} = N

_ism(a::TambaraYamagami) = a.n == modulus(a) # Checks whether a is the non-invertible
_chi(a::I, b::I) where {I <: TambaraYamagami} = cispi(2 * a.n * b.n / modulus(I)) # Non-degenerate symmetric bicharacter on ℤ_N

Base.length(::SectorValues{I}) where {I <: TambaraYamagami} = modulus(I) + 1
Base.IteratorSize(::Type{SectorValues{I}}) where {I <: TambaraYamagami} = HasLength()
Base.@propagate_inbounds function Base.getindex(v::SectorValues{I}, i::Int) where {I <: TambaraYamagami}
    @boundscheck 1 <= i <= length(v) || throw(BoundsError(v, i))
    return I(i - 1)
end
findindex(::SectorValues{I}, c::I) where {I <: TambaraYamagami} = Int(c.n) + 1
Base.IteratorSize(::Type{<:SectorProductIterator{I}}) where {I <: TambaraYamagami} = HasLength()

function Base.length(it::SectorProductIterator{I}) where {I <: TambaraYamagami}
    return (_ism(it.a) && _ism(it.b)) ? modulus(I) : 1
end
function Base.iterate(::SectorValues{I}, i::Int = 0) where {I <: TambaraYamagami}
    return i > modulus(I) ? nothing : (I(i), i + 1)
end
function Base.iterate(it::SectorProductIterator{I}, state::Int = 0) where {I <: TambaraYamagami}
    a, b = it.a, it.b
    am, bm = _ism(a), _ism(b)
    N = modulus(I)
    if am && bm
        state == N && return nothing
        return I(state), state + 1
    else
        state == 0 || return nothing
        c = (am || bm) ? I(N) : I(mod(a.n + b.n, N))
        return c, 1
    end
end

Base.isless(a1::I, a2::I) where {I <: TambaraYamagami} = isless(a1.n, a2.n)
dim(a::TambaraYamagami) = _ism(a) ? sqrt(float(modulus(a))) : 1.0
unit(::Type{I}) where {I <: TambaraYamagami} = I(0)
dual(a::TambaraYamagami) = _ism(a) ? a : typeof(a)(mod(- Int(a.n), modulus(a)))

FusionStyle(::Type{<:TambaraYamagami}) = SimpleFusion()
BraidingStyle(::Type{<:TambaraYamagami}) = NoBraiding()
fusionscalartype(::Type{<:TambaraYamagami}) = ComplexF64

function Nsymbol(a::I, b::I, c::I) where {N, I <: TambaraYamagami{N}}
    am, bm, cm = _ism(a), _ism(b), _ism(c)
    if am && bm
        return !cm
    elseif am || bm
        return cm
    else
        return !cm && c.n == mod(a.n + b.n, N)
    end
end

function Fsymbol(a::I, b::I, c::I, d::I, e::I, f::I) where {N, K, I <: TambaraYamagami{N, K}}
    T = fusionscalartype(I)

    (Nsymbol(a, b, e) && Nsymbol(e, c, d) && Nsymbol(b, c, f) && Nsymbol(a, f, d)) || return zero(T)

    am, bm, cm = _ism(a), _ism(b), _ism(c)

    if am && bm && cm # F^{mmm}_m
        return ((1 - 2 * K) / sqrt(N)) * conj(_chi(e, f))
    elseif !am && bm && !cm # F^{gmh}_{m}
        return _chi(a, c)
    elseif am && !bm && cm # F^{mgm}_{h}
        return _chi(b, d)
    else # F^{abc}_{a+b+c}
        return one(T)
    end
end

function Base.show(io::IO, a::TambaraYamagami)
    print_type = get(io, :typeinfo, nothing) !== typeof(a)
    print_type && print(io, type_repr(typeof(a)), "(")
    print(io, _ism(a) ? ":m" : Int(a.n))
    print_type && print(io, ")")
    return nothing
end
