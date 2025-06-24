I = IsingBimod
Istr = TensorKitSectors.type_repr(I)
@testset "$Istr sector" begin
    @testset "Basic type properties" begin
        @test eval(Meta.parse(sprint(show, I))) == I
        @test eval(Meta.parse(TensorKitSectors.type_repr(I))) == I
    end

    M = IsingBimod(CatType(3), 0)
    Mop = IsingBimod(CatType(2), 0)
    C0 = IsingBimod(CatType(1), 0)
    C1 = IsingBimod(CatType(1), 1)
    D0 = IsingBimod(CatType(4), 0)
    D1 = IsingBimod(CatType(4), 1)
    C = rand([C0, C1])
    D = rand([D0, D1])
    s = rand((M, Mop, C, D))

    @testset "Basic properties" begin
        @test @constinferred(one(C1)) == @constinferred(leftone(C1)) ==
              @constinferred(rightone(C1))
        @test one(D1) == leftone(D1) == rightone(D1)
        @test one(C1) == leftone(M) == rightone(Mop)
        @test one(D1) == rightone(M) == leftone(Mop)

        @test eval(Meta.parse(sprint(show, s))) == s
        @test @constinferred(hash(s)) == hash(deepcopy(s))
        @constinferred dual(s)
        @test dual(dual(s)) == s
        @constinferred dim(s)
        @constinferred frobeniusschur(s)
        @constinferred convert(IsingAnyon, s)

        @constinferred Bsymbol(C, C, C)
        @constinferred Fsymbol(D, D, D, D, D, D)
    end

    @testset "$Istr: Value iterator" begin
        @test eltype(values(I)) == I
        @test_throws ArgumentError one(I)
        sprev = C0 # first in SectorValues
        for (i, s) in enumerate(values(I))
            @test !isless(s, sprev) # confirm compatibility with sort order
            @test s == @constinferred (values(I)[i])
            @test findindex(values(I), s) == i
            sprev = s
            i >= 10 && break
        end
        @test C0 == first(values(I))
        @test (@constinferred findindex(values(I), C0)) == 1
        for s in smallset(I)
            @test (@constinferred values(I)[findindex(values(I), s)]) == s
        end
    end

    @testset "$Istr Fusion rules and F-symbols" begin
        argerr = ArgumentError("invalid fusion channel")
        # forbidden fusions
        for obs in [(C, D), (D, C), (M, M), (Mop, Mop), (D, M), (M, C), (Mop, D), (C, Mop)]
            @test isempty(⊗(obs...))
            @test_throws argerr Nsymbol(obs..., s)
        end

        # allowed fusions
        for obs in [(C, C), (D, D), (M, Mop), (Mop, M), (C, M), (Mop, C), (M, D), (D, Mop)]
            @test !isempty(⊗(obs...))
        end

        @test Nsymbol(C, C, one(C)) == Nsymbol(D, D, one(D)) == 1
        @test Nsymbol(C, M, M) == Nsymbol(Mop, C, Mop) == 1
        @test Nsymbol(M, D, M) == Nsymbol(D, Mop, Mop) == 1

        @test_throws argerr Nsymbol(M, Mop, D)
        @test_throws argerr Nsymbol(Mop, M, C)
        @test Nsymbol(M, Mop, C) == Nsymbol(Mop, M, D) == 1

        # non-trivial F-symbol checks
        @test Fsymbol(Mop, M, D1, D0, D0, M) == 0 # ℳᵒᵖ x ℳ x 𝒟 → ℳ allowed, but 𝒟-labels mismatch
        @test Fsymbol(Mop, M, D1, D0, D1, M) == 1
        @test Fsymbol(M, Mop, C1, C0, C0, Mop) == 0
        @test Fsymbol(M, Mop, C1, C0, C1, Mop) == 1

        @test Fsymbol(C, M, D, M, M, M) == (C.label * D.label == 0 ? 1 : -1) # 𝒞 x ℳ x 𝒟 → ℳ allowed
        @test_throws argerr Fsymbol(M, Mop, M, Mop, C, D) == 0 # IsingAnyon conversion would give non-zero
        @test_throws argerr Fsymbol(Mop, M, Mop, M, D, C) == 0

        @test Fsymbol(M, Mop, M, M, C, D) ==
              (C.label * D.label == 0 ? inv(sqrt(2)) : -inv(sqrt(2))) # ℳ x ℳᵒᵖ x ℳ → ℳ allowed
        @test Fsymbol(Mop, M, Mop, Mop, D, C) ==
              (C.label * D.label == 0 ? inv(sqrt(2)) : -inv(sqrt(2))) # ℳᵒᵖ x ℳ x ℳᵒᵖ → ℳᵒᵖ allowed

        @test_throws argerr Fsymbol(M, Mop, M, Mop, C, D)
    end

    @testset "$Istr: Unitarity of F-move" begin
        for a in smallset(I), b in smallset(I), c in smallset(I)
            for d in ⊗(a, b, c)
                es = collect(intersect(⊗(a, b), map(dual, ⊗(c, dual(d)))))
                fs = collect(intersect(⊗(b, c), map(dual, ⊗(dual(d), a))))
                @test length(es) == length(fs)
                F = [Fsymbol(a, b, c, d, e, f) for e in es, f in fs]
                @test isapprox(F' * F, one(F); atol=1e-12, rtol=1e-12)
            end
        end
    end
    @testset "$Istr: Pentagon equation" begin
        objects = collect(values(I))
        len(x) = length(collect(values(x))) # not predicting length of IsingBimodIterator
        for a in objects, b in objects
            len(a ⊗ b) > 0 || continue # skip if not compatible
            for c in objects
                len(b ⊗ c) > 0 || continue # skip if not compatible
                for d in objects
                    len(⊗(a, b, c)) > 0 || continue # skip if not compatible
                    @test pentagon_equation(a, b, c, d; atol=1e-12, rtol=1e-12)
                end
            end
        end
    end
end
