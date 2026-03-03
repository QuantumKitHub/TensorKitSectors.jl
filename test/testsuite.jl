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
* [`hasfusiontensor`](@ref)
"""
module SectorTestSuite

export smallset, randsector, hasfusiontensor, F_unitarity_test

using Test
using TestExtras
using TensorKitSectors
using TensorKitSectors: type_repr
using Random
using StatsBase
using Base.Iterators: take, product

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

function StatsBase.sample(::SectorValues{I}, size::Int) where {I <: Sector}
    Base.IteratorSize(values(I)) === Base.IsInfinite() &&
        throw(ArgumentError("Cannot take random sample of infinite sector values."))
    return sample(collect(values(I)), size; replace = false) # unique sampling
end
function smallset(::Type{I}, size::Int = 5) where {I <: Sector}
    vals = values(I)
    Base.IteratorSize(vals) === Base.IsInfinite() && return take(vals, size)
    return sample(vals, min(size, length(vals)))
end
function smallset(::Type{ProductSector{Tuple{I1, I2}}}) where {I1, I2}
    iter = product(smallset(I1), smallset(I2))
    s = collect(i ⊠ j for (i, j) in iter if dim(i) * dim(j) <= 6)
    return length(s) > 6 ? rand(s, 6) : s
end
function smallset(::Type{ProductSector{Tuple{I1, I2, I3}}}) where {I1, I2, I3}
    iter = product(smallset(I1), smallset(I2), smallset(I3))
    s = collect(i ⊠ j ⊠ k for (i, j, k) in iter if dim(i) * dim(j) * dim(k) <= 6)
    return length(s) > 6 ? rand(s, 6) : s
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

function hasfusiontensor(I::Type{<:Sector})
    try
        fusiontensor(unit(I), unit(I), unit(I))
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

include("sectors.jl")

end # module SectorTestSuite
