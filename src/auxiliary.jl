# kronecker product with arbitrary number of dimensions
_kron(A, B, C, D...) = _kron(_kron(A, B), C, D...)
function _kron(A, B)
    sA = size(A)
    sB = size(B)
    s = map(*, sA, sB)
    C = similar(A, promote_type(eltype(A), eltype(B)), s)
    for IA in eachindex(IndexCartesian(), A)
        for IB in eachindex(IndexCartesian(), B)
            I = CartesianIndex(IB.I .+ (IA.I .- 1) .* sB)
            C[I] = A[IA] * B[IB]
        end
    end
    return C
end
# combine the symbols `A` and `B` of the factors `a` and `b` of product sectors; a scalar symbol
# is promoted to an array of the same kind as the other symbol, with the size given by the
# multiplicities of its sectors, unless these are all one: then it just scales the other symbol,
# which is returned as is for a scalar one (symbols are not mutated, and may already be shared
# with a cache); the type of the result thus only depends on the types of the symbols
_kron_promote(A::Number, B::Number, a, b) = A * B
function _kron_promote(A::Number, B::AbstractArray, a, b)
    sz = _symbol_size(a)
    all(isone, sz) || return _kron(fill!(similar(B, typeof(A), sz), A), B)
    return isone(A) && promote_type(typeof(A), eltype(B)) === eltype(B) ? B : A * B
end
function _kron_promote(A::AbstractArray, B::Number, a, b)
    sz = _symbol_size(b)
    all(isone, sz) || return _kron(A, fill!(similar(A, typeof(B), sz), B))
    return isone(B) && promote_type(eltype(A), typeof(B)) === eltype(A) ? A : A * B
end
_kron_promote(A::AbstractArray, B::AbstractArray, a, b) = _kron(A, B)

# Manhattan based distance enumeration: I is supposed to be one-based index
# Points are zero-based coordinates `x` with `0 ≤ x[i] < sz[i]`, ordered first by their
# Manhattan distance `d = sum(x)` to the origin, and then lexicographically.

# size of an unbounded dimension (e.g. the values of an infinite sector type), whose bound can
# then be left out at compile time
struct Unbounded end
_bounded(::Tuple{}) = ()
_bounded(sz::Tuple{Unbounded, Vararg}) = _bounded(Base.tail(sz))
_bounded(sz::Tuple{Int, Vararg}) = (sz[1], _bounded(Base.tail(sz))...)
_maxcoordinate(sz::Int, d::Int) = min(sz - 1, d)
_maxcoordinate(::Unbounded, d::Int) = d

# Number of points `x` with `0 ≤ x[i] < sz[i]` and Manhattan distance `sum(x) < d`.
#
# Without the upper bounds, the number of points `x ≥ 0` with `sum(x) ≤ d - 1` is
# `binomial(d - 1 + N, N)` (stars and bars). The points that violate the upper bound in
# (at least) a subset `S` of the dimensions, i.e. with `x[i] ≥ sz[i]` for `i ∈ S`, are
# counted in the same way after shifting `x[i] → x[i] - sz[i]`, which gives
# `binomial(d - 1 - s + N, N)` with `s = sum(sz[i] for i in S)`. By inclusion-exclusion,
# the number of points that satisfy all bounds is then
#     ∑_S (-1)^|S| binomial(d - 1 - s + N, N)
# where the subsets `S` are enumerated as the bits of `mask`, and terms with `s ≥ d` vanish.
# Unbounded dimensions (size `Unbounded()`) never violate their bound, and are left out of the
# subsets at compile time.
@inline function num_manhattan_points_below(d::Int, sz::NTuple{N, Any}) where {N}
    d ≤ 0 && return 0
    bsz = _bounded(sz)
    num = 0
    for mask in 0:((1 << length(bsz)) - 1)
        s = 0
        sign = 1
        for i in 1:length(bsz)
            if isodd(mask >> (i - 1))
                s += min(bsz[i], d) # clamp to avoid overflow, terms with s ≥ d vanish anyway
                sign = -sign
            end
        end
        s < d && (num += sign * _binomial(d - 1 - s + N, Val(N)))
    end
    return num
end

# binomial(n, N) for n ≥ N, unrolled for small static N
# (this is faster than a lookup table, since the divisions by constants are cheap)
function _binomial(n::Int, ::Val{N}) where {N}
    r = 1
    for k in 1:N
        r = (r * (n - N + k)) ÷ k
    end
    return r
end

# forward mapping from multidimensional to single Manhattan index
@inline localoffset(d, I::Tuple{Int}, sz::Tuple{Any}) = 0
@inline function localoffset(d, I::NTuple{N, Int}, sz::NTuple{N, Any}) where {N}
    # number of points at distance `d` with a smaller first coordinate
    sz′ = Base.tail(sz)
    offset = num_manhattan_points_below(d + 1, sz′) - num_manhattan_points_below(d - I[1] + 2, sz′)
    return offset + localoffset(d - I[1] + 1, Base.tail(I), sz′)
end

to_manhattan_index(::Tuple{}, ::Tuple{}) = 1
to_manhattan_index(I::Tuple{Int}, ::Tuple{Any}) = I[1]
function to_manhattan_index(I::Dims{N}, sz::NTuple{N, Any}) where {N}
    d = sum(I) - N # Manhattan distance to origin for one-based indices
    return 1 + num_manhattan_points_below(d, sz) + localoffset(d, I, sz)
end

# inverse mapping
invertlocaloffset(d, offset, sz::Tuple{Any}) = (d + 1,)
@inline function invertlocaloffset(d, offset, sz::Tuple{Any, Any, Vararg})
    # find the largest first coordinate `x₁` such that the number of points at distance `d`
    # with a smaller first coordinate does not exceed `offset`
    sz′ = Base.tail(sz)
    nd = num_manhattan_points_below(d + 1, sz′)
    lo, hi = 0, _maxcoordinate(sz[1], d)
    while lo < hi
        mid = (lo + hi + 1) >> 1
        if nd - num_manhattan_points_below(d - mid + 1, sz′) ≤ offset
            lo = mid
        else
            hi = mid - 1
        end
    end
    offset -= nd - num_manhattan_points_below(d - lo + 1, sz′)
    return (lo + 1, invertlocaloffset(d - lo, offset, sz′)...)
end

# successor of a multidimensional index in Manhattan order, assuming it exists
# (here an unbounded size can be any size that is never reached)
next_manhattan_index(I::Dims{N}, sz::NTuple{N, Any}) where {N} =
    next_manhattan_index(I, map(s -> s isa Unbounded ? typemax(Int) : s, sz))
function next_manhattan_index(I::Dims{N}, sz::Dims{N}) where {N}
    # find the last position that can be incremented while lowering the distance of the
    # positions after it, and otherwise move to the next layer
    r = 0 # Manhattan distance of the positions after k
    k = 0
    @inbounds for j in N:-1:1
        if r > 0 && I[j] < sz[j]
            k = j
            break
        end
        r += I[j] - 1
    end
    if k > 0
        I = @inbounds Base.setindex(I, I[k] + 1, k)
        r -= 1
    else
        r += 1
    end
    # the smallest index with distance `r` for the positions after k is obtained by
    # filling them up from the back
    @inbounds for j in N:-1:(k + 1)
        x = min(r, sz[j] - 1)
        I = Base.setindex(I, x + 1, j)
        r -= x
    end
    return I
end

manhattan_to_multidimensional_index(index::Int, ::Tuple{}) = ()
manhattan_to_multidimensional_index(index::Int, ::Tuple{Any}) = (index,)
function manhattan_to_multidimensional_index(index::Int, sz::NTuple{N, Any}) where {N}
    # find the layer of the Manhattan distance, i.e. the smallest `d` such that there are at
    # least `index` points with distance `≤ d`, by exponential followed by binary search
    lo, hi = 0, 1
    while num_manhattan_points_below(hi + 1, sz) < index
        lo, hi = hi + 1, 2hi
    end
    while lo < hi
        mid = (lo + hi) >> 1
        if num_manhattan_points_below(mid + 1, sz) < index
            lo = mid + 1
        else
            hi = mid
        end
    end
    # find the position within the layer
    return invertlocaloffset(lo, index - 1 - num_manhattan_points_below(lo, sz), sz)
end
