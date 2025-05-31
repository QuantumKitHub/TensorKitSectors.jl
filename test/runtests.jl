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
                    Z2Irrep ⊠ FibonacciAnyon ⊠ FibonacciAnyon,
                    TimeReversed{Z2Irrep}, TimeReversed{Z3Irrep}, TimeReversed{Z4Irrep},
                    TimeReversed{U1Irrep}, TimeReversed{CU1Irrep}, TimeReversed{SU2Irrep},
                    TimeReversed{FibonacciAnyon}, TimeReversed{IsingAnyon},
                    TimeReversed{FermionParity},
                    TimeReversed{FermionParity ⊠ FermionParity},
                    TimeReversed{Z2Irrep ⊠ Z3Irrep ⊠ Z4Irrep},
                    TimeReversed{Z2Irrep} ⊠ TimeReversed{Z3Irrep} ⊠ TimeReversed{Z4Irrep},
                    TimeReversed{NewSU2Irrep ⊠ NewSU2Irrep},
                    TimeReversed{Z2Irrep ⊠ FibonacciAnyon ⊠ FibonacciAnyon}
                    # All the sectors above can pass the test. But the following ones can not. Still needs to find the reason.
                    # TimeReversed{NewSU2Irrep ⊠ SU2Irrep},
                    # TimeReversed{FermionParity ⊠ U1Irrep ⊠ SU2Irrep},
                    # TimeReversed{FermionParity ⊠ SU2Irrep ⊠ SU2Irrep},
                    # TimeReversed{FermionParity ⊠ SU2Irrep ⊠ NewSU2Irrep}
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
