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

const sectorlist = (
    Trivial, PlanarTrivial,
    Z2Irrep, Z3Irrep, Z4Irrep, Irrep[ℤ{200}], U1Irrep,
    DNIrrep{3}, DNIrrep{4}, DNIrrep{5}, CU1Irrep,
    SU2Irrep, NewSU2Irrep,
    FibonacciAnyon, IsingAnyon, FermionParity,
    FermionParity ⊠ FermionParity,
    Z3Irrep ⊠ Z4Irrep, FermionParity ⊠ U1Irrep ⊠ SU2Irrep,
    FermionParity ⊠ SU2Irrep ⊠ SU2Irrep, NewSU2Irrep ⊠ NewSU2Irrep,
    NewSU2Irrep ⊠ SU2Irrep, FermionParity ⊠ SU2Irrep ⊠ NewSU2Irrep,
    FibonacciAnyon ⊠ FibonacciAnyon ⊠ Z2Irrep,
    Z2Element{0}, Z2Element{1},
    Z3Element{0}, Z3Element{1}, Z3Element{2},
    Z4Element{0}, Z4Element{1}, Z4Element{2},
    Z3Element{1} ⊠ SU2Irrep,
    FibonacciAnyon ⊠ Z4Element{3},
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
    TimeReversed{FermionParity ⊠ SU2Irrep ⊠ NewSU2Irrep},
)

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

        @test @constinferred(length(allunits(I1 ⊠ I2))) == 1
    end
    @test @constinferred(Tuple(SU2Irrep(1) ⊠ U1Irrep(0))) == (SU2Irrep(1), U1Irrep(0))
    @test @constinferred(length(FermionParity(1) ⊠ SU2Irrep(1 // 2) ⊠ U1Irrep(1))) == 3
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

@testset "SectorProduct printing" begin
    a = SU2Irrep(1)
    @test repr(a ⊗ a) == "Irrep[SU₂](1) ⊗ Irrep[SU₂](1)"
    @test repr("text/plain", a ⊗ a) == "Irrep[SU₂](1) ⊗ Irrep[SU₂](1):\n 0\n 1\n 2"

    b = D4Irrep(1)
    @test repr(b ⊗ b) == "Irrep[D₄](1, false) ⊗ Irrep[D₄](1, false)"
    @test repr("text/plain", b ⊗ b) == "Irrep[D₄](1, false) ⊗ Irrep[D₄](1, false):\n (0, false)\n (0, true)\n (2, false)\n (2, true)"
end

@testset "ZNIrrep edge cases" begin
    for N in (127, 128, 129)
        a = Irrep[Cyclic{N}](-1)
        @test typeof(a) == ((N ≤ (typemax(UInt8) + 1) ÷ 2) ? ZNIrrep{N} : LargeZNIrrep{N})
        @test typeof(charge(a)) == Int
        @test charge(dual(a)) == mod(-charge(a), N)
        @test charge(only(a ⊗ a)) == mod(charge(a) + charge(a), N)
    end
end

include("multifusion.jl")

@testset "Aqua" begin
    using Aqua: Aqua
    Aqua.test_all(TensorKitSectors)
end

@testset "JET" begin
    using JET: JET
    JET.test_package(TensorKitSectors; target_defined_modules = true)
end
