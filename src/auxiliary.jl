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

# forward mapping from multidimensional to single Manhattan index
@inline function num_manhattan_points(d::Int, sz::Dims{1})
    return Int(sz[1] > d)
end

# when sz[k] > d for every k, then x_k < sz[k] - 1 is never violated by any split of d units among N axes (one x_k is maximally d, and sz[k] - 1 ≥ d covers that)
# so the bound drops -> calculate ways to write d as an ordered sum of N nonnegative integers,
# which is the number of compositions of d into N parts
@inline function _compositions(d::Int, ::Val{N}) where {N}
    try
        return binomial(d + N - 1, N - 1)
    catch e
        e isa OverflowError || rethrow()
        return nothing # fallback to recursion
    end
end
@inline function num_manhattan_points(d::Int, sz::Dims{N}) where {N}
    d == 0 && return 1
    if N ≥ 3 && all(>(d), sz)
        c = _compositions(d, Val(N))
        isnothing(c) || return c
    end
    num = 0
    for i in 1:min(sz[1], d + 1) # for N = 2 the loop is already linear
        num += num_manhattan_points(d - i + 1, Base.tail(sz))
    end
    return num
end

@inline localoffset(d, I::Tuple{Int}, sz::Tuple{Int}) = 0
@inline function localoffset(d, I::NTuple{N, Int}, sz::NTuple{N, Int}) where {N}
    offset = 0
    for i in 1:(I[1] - 1)
        offset += num_manhattan_points(d - i + 1, Base.tail(sz))
    end
    offset += localoffset(d - I[1] + 1, Base.tail(I), Base.tail(sz))
    return offset
end

to_manhattan_index(::Tuple{}, ::Tuple{}) = 1
to_manhattan_index(I::Tuple{Int}, ::Tuple{Int}) = I[1]
function to_manhattan_index(I::Dims{N}, sz::Dims{N}) where {N}
    d = sum(I) - N # Manhattan distance to origin for one-based indices
    d == 0 && return 1
    index = 1
    # count all the points with smaller Manhatten distance
    for k in 0:(d - 1)
        index += num_manhattan_points(k, sz)
    end
    return index + localoffset(d, I, sz)
end

# inverse mapping
@inline invertlocaloffset(d, offset, sz::Tuple{Int}) = (d + 1,)
@inline function invertlocaloffset(d, offset, sz::Tuple{Int, Int, Vararg{Int}})
    i₁ = 1
    while i₁ < sz[1]
        jump = num_manhattan_points(d - i₁ + 1, Base.tail(sz))
        if offset < jump
            break
        end
        offset -= jump
        i₁ += 1
    end
    return (i₁, invertlocaloffset(d - i₁ + 1, offset, Base.tail(sz))...)
end

manhattan_to_multidimensional_index(index::Int, ::Dims{0}) = ()
manhattan_to_multidimensional_index(index::Int, ::Dims{1}) = (index,)
function manhattan_to_multidimensional_index(index::Int, sz::Dims{N}) where {N}
    # find the layer of the Manhattan distance
    index == 1 && return ntuple(one, Val(N))
    offset = 1
    d = 1
    while true # reason for the boundscheck
        currentlayer = num_manhattan_points(d, sz)
        if index <= offset + currentlayer
            break
        end
        d += 1
        offset += currentlayer
    end

    # find the position within the layer
    index -= offset
    return invertlocaloffset(d, index - 1, sz)
end

# tabulation of the mapping to avoid repeated calls to manhattan_to_multidimensional_index
# and to_manhattan_index
struct MTable{N}
    multi::Vector{NTuple{N, Int}} # Manhattan index -> multi-index
    lin::Array{Int, N} # multi-index -> Manhattan index
end

function build_table(sz::Dims{N}) where {N}
    n = prod(sz)
    multi = Vector{NTuple{N, Int}}(undef, n)
    k = 0
    for J in CartesianIndices(sz)
        k += 1
        multi[k] = Tuple(J)
    end
    sort!(multi; by = J -> (sum(J), J)) # sort first by distance, matches isless on ProductSector
    lin = Array{Int, N}(undef, sz)
    for i in 1:n
        lin[CartesianIndex(multi[i])] = i
    end
    return MTable{N}(multi, lin)
end

const TABLES = IdDict{DataType, Any}()
const TABLE_LOCK = ReentrantLock() # dictionaries aren't thread-safe for concurrent mutation
# within TensorKit no problem: OhMyThreads allows concurrent calls into e.g. sectors(V)

@inline function mtable(::Type{I}, sz::Dims{N}) where {I, N}
    t = get(TABLES, I, nothing)
    if t === nothing
        t = @lock TABLE_LOCK get!(() -> build_table(sz), TABLES, I)
    end
    return t::MTable{N} # re-establish concrete typing for downstream
end

# refuse to tabulate absurdly large label sets
# GradedSpace NTuple storage is only used for finite sectors
# not sure what a reasonable limit is, might be edited to smaller value
# depending on findings of NTuple vs Dict performance tests
const MAXTABLE = 2^20
