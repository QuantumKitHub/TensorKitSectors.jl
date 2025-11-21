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

SectorTestSuite.test_sectortype(MySectorType)
```

Additionally, this test suite exports the following convenience testing utilities:
* [`smallset`](@ref)
* [`randsector`](@ref)
* [`hasfusiontensor`](@ref)
"""
module SectorTestSuite

export smallset, randsector, hasfusiontensor

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
        @assert !haskey(tests, $name)
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

"""
    @testinferred [AllowedTypes] ex

Like `Test.@inferred`, but registers failures through a test, rather than an error.
"""
macro testinferred(ex)
    return _inferred(ex, __module__)
end
macro testinferred(ex, allow)
    return _inferred(ex, __module__, allow)
end

# Implementation copied from Test._inferred:
function _inferred(ex, mod, allow = :(Union{}))
    if Meta.isexpr(ex, :ref)
        ex = Expr(:call, :getindex, ex.args...)
    end
    Meta.isexpr(ex, :call)|| error("@testinferred requires a call expression")
    farg = ex.args[1]
    if isa(farg, Symbol) && farg !== :.. && first(string(farg)) == '.'
        farg = Symbol(string(farg)[2:end])
        ex = Expr(
            :call, GlobalRef(Test, :_materialize_broadcasted),
            farg, ex.args[2:end]...
        )
    end
    result = let ex = ex
        quote
            let allow = $(esc(allow))
                allow isa Type || throw(ArgumentError("@inferred requires a type as second argument"))
                $(
                    if any(@nospecialize(a) -> (Meta.isexpr(a, :kw) || Meta.isexpr(a, :parameters)), ex.args)
                        # Has keywords
                        # Create the call expression with escaped user expressions
                        call_expr = :($(esc(ex.args[1]))(args...; kwargs...))
                        quote
                            args, kwargs, result = $(esc(Expr(:call, _args_and_call, ex.args[2:end]..., ex.args[1])))
                            # wrap in dummy hygienic-scope to work around scoping issues with `call_expr` already having `esc` on the necessary parts
                            inftype = $(Expr(:var"hygienic-scope", Base.gen_call_with_extracted_types(mod, Base.infer_return_type, call_expr; is_source_reflection = false), Test))
                        end
                    else
                        # No keywords
                        quote
                            args = ($([esc(ex.args[i]) for i in 2:length(ex.args)]...),)
                            result = $(esc(ex.args[1]))(args...)
                            inftype = Base.infer_return_type($(esc(ex.args[1])), Base.typesof(args...))
                        end
                    end
                )
                rettype = result isa Type ? Type{result} : typeof(result)
                @test rettype <: allow || rettype == Base.typesplit(inftype, allow)
                result
            end
        end
    end
    return Base.remove_linenums!(result)
end

smallset(::Type{I}) where {I <: Sector} = take(values(I), 5)
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

include("sectors.jl")

end # module SectorTestSuite
