"""
    module SectorTestSuite

Lightweight testsuite registration for sector tests, inspired by the GPUArrays
testsuite style. Each logical group of tests is registered under a string key
and can be iterated for every sector type.
"""
module SectorTestSuite

export tests, @testsuite, @testinferred

using Test
using TestExtras
using TensorKitSectors
const TKS = TensorKitSectors

const tests = Dict()

"""
    @testsuite name I -> begin
        # test code here
    end

Register a sector testsuite. The body is executed with a single argument `I`, the concrete `Sector` type under test.
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
Runs the entire TensorKitSectors test suite on sector type `I`
"""
function test(I::Type)
    return @testset "$(TKS.type_repr(I))" begin
        for (name, fun) in tests
            code = quote
                $fun($I)
            end
            @eval @testset $name $code
        end
    end
end

macro testinferred(ex)
    return _inferred(ex, __module__)
end
macro testinferred(ex, allow)
    return _inferred(ex, __module__, allow)
end

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

include("testsetup.jl")
include("sectors.jl")

end # module SectorTestSuite
