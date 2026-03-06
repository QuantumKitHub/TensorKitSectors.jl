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
* [`eval_module`](@ref)
"""
module SectorTestSuite

export smallset, randsector, random_fusion, hasfusiontensor, F_unitarity_test
export eval_module

using Test
using TestExtras
using TensorKitSectors
using TensorKitSectors: type_repr
using Random
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
        while isempty(⊗(sprev, s)) && counter < 20
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
    eval_module(a::I; _module = parentmodule(I)) where {I <: Sector}

Passes the sector `a` through the `show` machinery to generate a string representation, then parses and evaluates the resulting expression in the module `_module`, by default the parent module of the sector.
"""
function eval_module(a::I; _module = parentmodule(I)) where {I <: Sector}
    if occursin("NewSU2Irrep", string(I)) # special case since NewSU2Irrep is defined in test setup, not in TensorKitSectors
        _module = Main
    end
    return Base.eval(_module, Meta.parse(sprint(show, I))) == I &&
        Base.eval(_module, Meta.parse(TensorKitSectors.type_repr(I))) == I &&
        Base.eval(_module, Meta.parse(sprint(show, a))) == a
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
