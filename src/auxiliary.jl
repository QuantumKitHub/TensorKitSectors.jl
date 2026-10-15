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
_kron_promote(A::Number, B::Number, _, _) = A * B
function _kron_promote(A₁, B₁, sz₁, sz₂)
    return _kron(A₁ isa Number ? fill(A₁, sz₁) : A₁, B₁ isa Number ? fill(B₁, sz₂) : B₁)
end

# Manhattan based distance enumeration: I is supposed to be one-based index
# Points are zero-based coordinates `x` with `0 ≤ x[i] < sz[i]`, ordered first by their
# Manhattan distance `d = sum(x)` to the origin, and then lexicographically.

# number of points with Manhattan distance smaller than `d`, computed in closed form through
# inclusion-exclusion over the dimensions where the upper bound `sz[i]` is violated
@inline function num_manhattan_points_below(d::Int, sz::Dims{N}) where {N}
    d ≤ 0 && return 0
    num = 0
    for mask in 0:((1 << N) - 1)
        s = 0
        sign = 1
        for i in 1:N
            if isodd(mask >> (i - 1))
                s += min(sz[i], d) # clamp to avoid overflow, terms with s ≥ d vanish anyway
                sign = -sign
            end
        end
        s < d && (num += sign * _binomial(d - 1 - s + N, Val(N)))
    end
    return num
end

# binomial(n, N) for n ≥ N, unrolled for small static N
@inline function _binomial(n::Int, ::Val{N}) where {N}
    r = 1
    for k in 1:N
        r = (r * (n - N + k)) ÷ k
    end
    return r
end

# number of points with Manhattan distance exactly `d`
@inline function num_manhattan_points(d::Int, sz::Dims)
    return num_manhattan_points_below(d + 1, sz) - num_manhattan_points_below(d, sz)
end

# forward mapping from multidimensional to single Manhattan index
@inline localoffset(d, I::Tuple{Int}, sz::Tuple{Int}) = 0
@inline function localoffset(d, I::NTuple{N, Int}, sz::NTuple{N, Int}) where {N}
    # number of points at distance `d` with a smaller first coordinate
    sz′ = Base.tail(sz)
    offset = num_manhattan_points_below(d + 1, sz′) - num_manhattan_points_below(d - I[1] + 2, sz′)
    return offset + localoffset(d - I[1] + 1, Base.tail(I), sz′)
end

to_manhattan_index(::Tuple{}, ::Tuple{}) = 1
to_manhattan_index(I::Tuple{Int}, ::Tuple{Int}) = I[1]
function to_manhattan_index(I::Dims{N}, sz::Dims{N}) where {N}
    d = sum(I) - N # Manhattan distance to origin for one-based indices
    return 1 + num_manhattan_points_below(d, sz) + localoffset(d, I, sz)
end

# inverse mapping
@inline invertlocaloffset(d, offset, sz::Tuple{Int}) = (d + 1,)
@inline function invertlocaloffset(d, offset, sz::Tuple{Int, Int, Vararg{Int}})
    # find the largest first coordinate `x₁` such that the number of points at distance `d`
    # with a smaller first coordinate does not exceed `offset`
    sz′ = Base.tail(sz)
    nd = num_manhattan_points_below(d + 1, sz′)
    lo, hi = 0, min(sz[1] - 1, d)
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

manhattan_to_multidimensional_index(index::Int, ::Dims{0}) = ()
manhattan_to_multidimensional_index(index::Int, ::Dims{1}) = (index,)
function manhattan_to_multidimensional_index(index::Int, sz::Dims{N}) where {N}
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
