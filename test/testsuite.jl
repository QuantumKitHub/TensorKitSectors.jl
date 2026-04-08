"""
    module SectorTestSuite

Test suite and utilities that ensure a reusable way of verifying the required interface for a `Sector` type.
Framework based on the GPUArrays testsuite.

Downstream packages may include this test suite as follows:

```julia
import TensorKitSectors
testsuite_path = joinpath(
    dirname(dirname(pathof(TensorKitSectors))), # TensorKitSectors root
    "test", "testsuite.jl"
)
include(testsuite_path)

SectorTestSuite.test_sector(MySectorType)
```

Additionally, this test suite exports the following convenience testing utilities:
* [`smallset`](@ref)
* [`randsector`](@ref)
* [`random_fusion`](@ref)
* [`hasfusiontensor`](@ref)
* [`F_unitarity_test`](@ref)
* [`R_unitarity_test`](@ref)
* [`can_fuse`](@ref)
"""
module SectorTestSuite

export smallset, randsector, random_fusion, hasfusiontensor, can_fuse
export F_unitarity_test, R_unitarity_test

using Test
using TestExtras
using TensorKitSectors
using TensorKitSectors: type_repr
using Random
using Base.Iterators: take

const tests = Dict()

"""
    @testsuite name I -> begin
        # test code here
    end

Register a sector testsuite.
The body is executed with a single argument `I`, the concrete `Sector` type under test.
"""
macro testsuite(name, ex)
    safe_name = lowercase(replace(replace(name, " " => "_"), "/" => "_"))
    fn = Symbol("test_$(safe_name)")
    return quote
        $(esc(fn))(I) = $(esc(ex))(I)
        @assert !haskey(tests, $name) "duplicate testsuite name: $name"
        tests[$name] = $fn
    end
end

"""
    test_sectortype(I::Type)

Runs the entire TensorKitSectors test suite on sector type `I`.
"""
function test_sector(I::Type)
    return @testset "$(type_repr(I))" begin
        for (name, fun) in tests
            code = quote
                $fun($I)
            end
            @eval @testset $name $code
        end
    end
end

function smallset(::Type{I}, size::Int = 5, maxdim::Real = 10) where {I <: Sector}
    sectors = collect(Iterators.take(values(I), 10 * size))
    sectors = shuffle!(filter!(s -> dim(s) < maxdim, sectors))
    return resize!(sectors, min(size, length(sectors)))
end

function randsector(::Type{I}) where {I <: Sector}
    s = collect(smallset(I))
    a = Random.rand(s)
    while isunit(a) # don't use trivial label
        a = Random.rand(s)
    end
    return a
end
randsector(::Type{I}) where {I <: Union{Trivial, PlanarTrivial}} = unit(I)

"""
    can_fuse(a::I, b::I) where {I <: Sector}

Returns a boolean indicating whether the fusion of `a` and `b` is allowed,
i.e. whether there exists a sector `c` such that `c ∈ ⊗(a, b)`.
"""
can_fuse(a::I, b::I) where {I <: Sector} = !isempty(⊗(a, b))

"""
    random_fusion(I::Type, N::Int)

Returns a vector of `N` random sectors from `I` that have a non-empty coupled sector.
Compatible with any `Sector` type, including those with `UnitStyle(I) == GenericUnit()`.
"""
function random_fusion(I::Type{<:Sector}, N::Int)
    vec = Vector{I}(undef, N)
    vec[1] = randsector(I)
    N == 1 && return vec
    for i in 2:N
        s = randsector(I)
        sprev = vec[i - 1]
        counter = 0
        while !can_fuse(sprev, s) && counter < 20
            counter += 1
            s = (counter < 20) ? randsector(I) : rightunit(sprev)
        end
        vec[i] = s
    end
    return vec
end

function hasfusiontensor(I::Type{<:Sector})
    try
        u = first(allunits(I))
        fusiontensor(u, u, u)
        return true
    catch e
        if e isa MethodError
            return false
        else
            rethrow(e)
        end
    end
end

"""
    F_unitarity_test(a::I, b::I, c::I; kwargs...) where {I <: Sector}

Tests the unitarity of the F-symbols for the fusion of `a`, `b`, and `c`.
Returns `true` if the F-symbols are unitary, and `false` otherwise.
"""
function F_unitarity_test(a::I, b::I, c::I; kwargs...) where {I <: Sector}
    for d in ⊗(a, b, c)
        es = collect(intersect(⊗(a, b), map(dual, ⊗(c, dual(d)))))
        fs = collect(intersect(⊗(b, c), map(dual, ⊗(dual(d), a))))
        if FusionStyle(I) isa MultiplicityFreeFusion
            @assert length(es) == length(fs)
            F = [Fsymbol(a, b, c, d, e, f) for e in es, f in fs]
        else
            Fblocks = Vector{Any}()
            for e in es, f in fs
                Fs = Fsymbol(a, b, c, d, e, f)
                push!(Fblocks, reshape(Fs, (size(Fs, 1) * size(Fs, 2), size(Fs, 3) * size(Fs, 4))))
            end
            F = hvcat(length(fs), Fblocks...)
        end
        isapprox(F' * F, one(F); kwargs...) || return false
    end
    return true
end

"""
    R_unitarity_test(a::I, b::I; kwargs...) where {I <: Sector}

Tests the unitarity of the R-symbols for the fusion of `a`, `b`, and `c`.
Returns `true` if the R-symbols are unitary, and `false` otherwise.
"""
function R_unitarity_test(a::I, b::I; kwargs...) where {I <: Sector}
    for c in ⊗(a, b)
        R = Rsymbol(a, b, c)
        isapprox(R' * R, one(R); kwargs...) || return false
    end
    return true
end

include("sectors.jl")

end # module SectorTestSuite
