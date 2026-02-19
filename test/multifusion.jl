I = IsingBimodule
Istr = TensorKitSectors.type_repr(I)
@testset "$Istr sector" begin
    @testset "Basic type properties" begin
        prodsec = I ⊠ Z2Irrep
        @test UnitStyle(prodsec) isa GenericUnit
        @test FusionStyle(prodsec) isa SimpleFusion
        @test_throws DomainError unit(prodsec)
        @test length(allunits(prodsec)) == 2
    end

    M = IsingBimodule(1, 2, 0)
    Mop = IsingBimodule(2, 1, 0)
    C0 = IsingBimodule(1, 1, 0)
    C1 = IsingBimodule(1, 1, 1)
    D0 = IsingBimodule(2, 2, 0)
    D1 = IsingBimodule(2, 2, 1)
    C = rand([C0, C1])
    D = rand([D0, D1])
    s = rand((M, Mop, C, D))

    @testset "Basic properties" begin
        @test @testinferred(unit(C1)) == @testinferred(leftunit(C1)) ==
            @testinferred(rightunit(C1))
        @test unit(D1) == leftunit(D1) == rightunit(D1)
        @test unit(C1) == leftunit(M) == rightunit(Mop)
        @test unit(D1) == rightunit(M) == leftunit(Mop)

        @test @testinferred(isunit(C0))
        @test isunit(D0)
        @test !isunit(C1) && !isunit(D1) && !isunit(M) && !isunit(Mop)

        @test length(allunits(I)) == 2
        @test allunits(I) == (C0, D0)

        @test !isunit(C0 ⊠ C1)
        @test !isunit(C0 ⊠ D1)
        @test isunit(C0 ⊠ C0)
        @test isunit(D0 ⊠ D0)
        @test isunit(C0 ⊠ D0)
        @test length(allunits(I ⊠ I)) == 4

        @test leftunit(M ⊠ Mop) == C0 ⊠ D0 == rightunit(Mop ⊠ M)
        @testinferred convert(IsingAnyon, s)
    end

    @testset "$Istr: Printing and errors" begin
        @test Base.eval(Meta.parse(sprint(show, C))) == C
        @test Base.eval(Meta.parse(sprint(show, M))) == M
        @test Base.eval(Meta.parse(sprint(show, Mop))) == Mop
        @test Base.eval(Meta.parse(sprint(show, D))) == D
        @test_throws DomainError unit(M)
        @test_throws DomainError unit(Mop)
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

        @test Nsymbol(C, C, unit(C)) == Nsymbol(D, D, unit(D)) == 1
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
end
