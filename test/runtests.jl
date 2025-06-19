using Test
using TestExtras
using Random
# using TensorKit: TensorKitSectors
using TensorKitSectors
using TensorOperations
using Base.Iterators: take, product
using LinearAlgebra: LinearAlgebra

const TKS = TensorKitSectors

include("testsetup.jl")
using .TestSetup
include("newsectors.jl")
using .NewSectors

const sectorlist = (Z2Irrep, Z3Irrep, Z4Irrep, U1Irrep, CU1Irrep, SU2Irrep, NewSU2Irrep,
                    FibonacciAnyon, IsingAnyon, FermionParity,
                    FermionParity ⊠ FermionParity,
                    Z3Irrep ⊠ Z4Irrep, FermionParity ⊠ U1Irrep ⊠ SU2Irrep,
                    FermionParity ⊠ SU2Irrep ⊠ SU2Irrep, NewSU2Irrep ⊠ NewSU2Irrep,
                    NewSU2Irrep ⊠ SU2Irrep, FermionParity ⊠ SU2Irrep ⊠ NewSU2Irrep,
                    FibonacciAnyon ⊠ FibonacciAnyon ⊠ Z2Irrep,
                    TimeReversed{Z2Irrep},
                    TimeReversed{Z3Irrep}, TimeReversed{Z4Irrep},
                    TimeReversed{U1Irrep}, TimeReversed{CU1Irrep}, TimeReversed{SU2Irrep},
                    TimeReversed{FibonacciAnyon}, TimeReversed{IsingAnyon},
                    TimeReversed{FermionParity},
                    TimeReversed{FermionParity ⊠ FermionParity},
                    TimeReversed{Z2Irrep ⊠ Z3Irrep ⊠ Z4Irrep},
                    TimeReversed{Z2Irrep} ⊠ TimeReversed{Z3Irrep} ⊠ TimeReversed{Z4Irrep},
                    TimeReversed{NewSU2Irrep ⊠ NewSU2Irrep},
                    TimeReversed{Z2Irrep ⊠ FibonacciAnyon ⊠ FibonacciAnyon},
                    TimeReversed{NewSU2Irrep ⊠ SU2Irrep},
                    TimeReversed{FermionParity ⊠ U1Irrep ⊠ SU2Irrep},
                    TimeReversed{FermionParity ⊠ SU2Irrep ⊠ SU2Irrep},
                    TimeReversed{FermionParity ⊠ SU2Irrep ⊠ NewSU2Irrep})

@testset "$(TensorKitSectors.type_repr(I))" for I in sectorlist
    @include("sectors.jl")
end

@testset "Deligne product" begin
    sectorlist′ = (Trivial, sectorlist...)
    for I1 in sectorlist′, I2 in sectorlist′
        a = first(smallset(I1))
        b = first(smallset(I2))

        @constinferred a ⊠ b
        @constinferred a ⊠ b ⊠ a
        @constinferred a ⊠ b ⊠ a ⊠ b
        @constinferred I1 ⊠ I2
        @test typeof(a ⊠ b) == I1 ⊠ I2
    end
end

@testset "Issue that came up in #11" begin
    a = Z2Irrep(1) ⊠ NewSU2Irrep(1)
    b = Z2Irrep(0) ⊠ NewSU2Irrep(1)
    c = Z2Irrep(0) ⊠ NewSU2Irrep(1)
    d = Z2Irrep(0) ⊠ NewSU2Irrep(1)
    e = Z2Irrep(0) ⊠ NewSU2Irrep(1)
    f = Z2Irrep(0) ⊠ NewSU2Irrep(1)
    @test size(Fsymbol(a, b, c, d, e, f)) == (0, 1, 1, 0)
    @test size(Bsymbol(a, b, c)) == (0, 0)
    @test size(Rsymbol(a, b, c)) == (0, 0)
    a = NewSU2Irrep(1) ⊠ Z2Irrep(1)
    b = NewSU2Irrep(1) ⊠ Z2Irrep(0)
    c = NewSU2Irrep(1) ⊠ Z2Irrep(0)
    d = NewSU2Irrep(1) ⊠ Z2Irrep(0)
    e = NewSU2Irrep(1) ⊠ Z2Irrep(0)
    f = NewSU2Irrep(1) ⊠ Z2Irrep(0)
    @test size(Fsymbol(a, b, c, d, e, f)) == (0, 1, 1, 0)
    @test size(Bsymbol(a, b, c)) == (0, 0)
    @test size(Rsymbol(a, b, c)) == (0, 0)
end

@testset "IsingBimod sector" begin
    I = IsingBimod
    Istr = TensorKitSectors.type_repr(I)
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

    @testset "Basic properties" begin
        @test @constinferred(one(C1)) == @constinferred(leftone(C1)) ==
              @constinferred(rightone(C1))
        @test one(D1) == leftone(D1) == rightone(D1)
        @test one(C1) == leftone(M) == rightone(Mop)
        @test one(D1) == rightone(M) == leftone(Mop)

        s = rand((M, Mop, C, D))
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

    @testset "Sector $Istr: Value iterator" begin
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

    @testset "Fusion rules and F-symbols" begin
        # forbidden fusions
        for obs in [(C, D), (D, C), (M, M), (Mop, Mop), (D, M), (M, C), (Mop, D), (C, Mop)]
            @test isempty(⊗(obs...))
            @test_throws ArgumentError("invalid fusion channel") Nsymbol(obs..., s)
        end

        # allowed fusions
        for obs in [(C, C), (D, D), (M, Mop), (Mop, M), (C, M), (Mop, C), (M, D), (D, Mop)]
            @test !isempty(⊗(obs...))
        end

        @test Nsymbol(C, C, one(C)) == Nsymbol(D, D, one(D)) == 1
        @test Nsymbol(C, M, M) == Nsymbol(Mop, C, Mop) == 1
        @test Nsymbol(M, D, M) == Nsymbol(Mop, D, Mop) == 1

        @test_throws ArgumentError("invalid fusion channel") Nsymbol(M, Mop, D)
        @test_throws ArgumentError("invalid fusion channel") Nsymbol(Mop, M, C)
        @test Nsymbol(M, Mop, C) == Nsymbol(Mop, M, D) == 1

        # non-trivial F-symbol checks
        @test Fsymbol(Mop, M, D1, D0, D0, M) == 0 # ℳᵒᵖ x ℳ x 𝒟 → ℳ allowed, but 𝒟-labels mismatch
        @test Fsymbol(Mop, M, D1, D0, D1, M) == 1
        @test Fsymbol(M, Mop, C1, C0, C0, Mop) == 0
        @test Fsymbol(M, Mop, C1, C0, C1, Mop) == 1

        @test Fsymbol(C, M, D, M, M, M) == 1 # 𝒞 x ℳ x 𝒟 → ℳ allowed
        @test Fsymbol(M, Mop, M, Mop, C, D) == 0 # IsingAnyon conversion would give non-zero
        @test Fsymbol(Mop, M, Mop, M, D, C) == 0

        @test Fsymbol(M, Mop, M, M, C, D) ==
              (C.label * D.label == 0 ? inv(sqrt(2)) : -inv(sqrt(2))) # ℳ x ℳᵒᵖ x ℳ → ℳ allowed
        @test Fsymbol(Mop, M, Mop, Mop, D, C) ==
              (C.label * D.label == 0 ? inv(sqrt(2)) : -inv(sqrt(2))) # ℳᵒᵖ x ℳ x ℳᵒᵖ → ℳᵒᵖ allowed

        @test_throws ArgumentError("invalid fusion channel") Fsymbol(M, Mop, M, Mop, C, D)
    end
end

@testset "Aqua" begin
    using Aqua: Aqua
    Aqua.test_all(TensorKitSectors)
end

@testset "JET" begin
    using JET: JET
    JET.test_package(TensorKitSectors; target_defined_modules=true)
end
