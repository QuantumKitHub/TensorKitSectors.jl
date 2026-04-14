using Test
using TestExtras
using TensorKitSectors

include("newsectors.jl")
using .NewSectors

const sectorlist = (
    Trivial, PlanarTrivial,
    Z2Irrep, Z3Irrep, Z4Irrep, Irrep[ℤ{200}], U1Irrep,
    DNIrrep{3}, DNIrrep{4}, DNIrrep{5}, CU1Irrep,
    A4Irrep, SU2Irrep, NewSU2Irrep,
    FibonacciAnyon, IsingAnyon, FermionParity,
    FermionParity ⊠ FermionParity, FibonacciAnyon ⊠ PlanarTrivial,
    Z3Irrep ⊠ Z4Irrep,
    @NamedSector{charge::U1Irrep, spin::SU2Irrep},
    @NamedSector{parity::FermionParity, charge::U1Irrep, spin::SU2Irrep},
    FermionParity ⊠ SU2Irrep ⊠ SU2Irrep, NewSU2Irrep ⊠ NewSU2Irrep,
    NewSU2Irrep ⊠ SU2Irrep, FermionParity ⊠ SU2Irrep ⊠ NewSU2Irrep,
    FibonacciAnyon ⊠ FibonacciAnyon ⊠ Z2Irrep,
    A4Irrep ⊠ Z2Irrep, A4Irrep ⊠ SU2Irrep,
    Z2Element{0}, Z2Element{1},
    Z3Element{0}, Z3Element{1}, Z3Element{2},
    Z4Element{0}, Z4Element{1}, Z4Element{2},
    ZNElement{6, 0}, ZNElement{6, 1}, ZNElement{6, 3},
    Z3Element{1} ⊠ SU2Irrep,
    FibonacciAnyon ⊠ Z4Element{3},
    IsingBimodule, IsingBimodule ⊠ IsingBimodule, IsingBimodule ⊠ Z2Irrep,
    IsingBimodule ⊠ SU2Irrep, IsingBimodule ⊠ FibonacciAnyon,
    TimeReversed{Z2Irrep},
    TimeReversed{Z3Irrep}, TimeReversed{Z4Irrep}, TimeReversed{A4Irrep},
    TimeReversed{U1Irrep}, TimeReversed{CU1Irrep}, TimeReversed{SU2Irrep},
    TimeReversed{FibonacciAnyon}, TimeReversed{IsingAnyon},
    TimeReversed{FermionParity},
    TimeReversed{FermionParity ⊠ FermionParity},
    TimeReversed{Z4Element{2}}, TimeReversed{ZNElement{6, 3}},
    ZNElement{6, 1} ⊠ TimeReversed{FibonacciAnyon},
    TimeReversed{ZNElement{6, 3} ⊠ FibonacciAnyon},
    TimeReversed{Z2Irrep ⊠ Z3Irrep ⊠ Z4Irrep},
    TimeReversed{Z2Irrep} ⊠ TimeReversed{Z3Irrep} ⊠ TimeReversed{Z4Irrep},
    TimeReversed{NewSU2Irrep ⊠ NewSU2Irrep},
    TimeReversed{Z2Irrep ⊠ FibonacciAnyon ⊠ FibonacciAnyon},
    TimeReversed{NewSU2Irrep ⊠ SU2Irrep},
    TimeReversed{FermionParity ⊠ U1Irrep ⊠ SU2Irrep},
    TimeReversed{@NamedSector{parity::FermionParity, charge::U1Irrep, spin::SU2Irrep}},
    TimeReversed{FermionParity ⊠ SU2Irrep ⊠ NewSU2Irrep},
)

include("testsuite.jl")
using .SectorTestSuite

@testset "Sector test suite" verbose = true begin
    for sectortype in sectorlist
        @time SectorTestSuite.test_sector(sectortype)
    end
end

@testset "Intertwiner relation for A4Irrep" begin
    ω = cis(2π / 3)
    T3 = [1 0 0; 0 ω 0; 0 0 ω^2]
    T(a::Int8) = (a == 3) ? T3 : hcat(ω^(a))
    S3 = 1 / 3 * [-1 2 2; 2 -1 2; 2 2 -1]
    S(a::Int8) = (a == 3) ? S3 : hcat(1)
    for a in smallset(A4Irrep), b in smallset(A4Irrep)
        for c in ⊗(a, b)
            C = fusiontensor(a, b, c)
            Ta, Tb, Tc = T(a.n), T(b.n), T(c.n)
            Sa, Sb, Sc = S(a.n), S(b.n), S(c.n)
            for μ in 1:Nsymbol(a, b, c)
                Cmat = reshape(view(C, :, :, :, μ), (dim(a) * dim(b), dim(c)))
                L_T = Cmat' * kron(Ta, Tb) * Cmat
                L_S = Cmat' * kron(Sa, Sb) * Cmat
                @test isapprox(L_T, Tc; atol = 1.0e-12)
                @test isapprox(L_S, Sc; atol = 1.0e-12)
            end
        end
    end
end

@testset "Deligne product" begin
    sectorlist′ = (Trivial, sectorlist...)
    for I1 in sectorlist′, I2 in rand(sectorlist′, 3)
        a = first(smallset(I1))
        b = first(smallset(I2))

        @testinferred a ⊠ b
        @testinferred a ⊠ b ⊠ a
        @testinferred a ⊠ b ⊠ a ⊠ b
        @testinferred I1 ⊠ I2
        @test typeof(a ⊠ b) == I1 ⊠ I2

        if UnitStyle(I1 ⊠ I2) isa SimpleUnit
            @test @testinferred(length(allunits(I1 ⊠ I2))) == 1
            @test @testinferred(unit(I1 ⊠ I2)) == leftunit(a ⊠ b) == rightunit(a ⊠ b)
        end
    end
    @test @testinferred(Tuple(SU2Irrep(1) ⊠ U1Irrep(0))) == (SU2Irrep(1), U1Irrep(0))
    @test @testinferred(length(FermionParity(1) ⊠ SU2Irrep(1 // 2) ⊠ U1Irrep(1))) == 3
end

@testset "NamedSector" begin
    CS = @NamedSector{charge::U1Irrep, spin::SU2Irrep}
    @test CS === NamedSector{@NamedTuple{charge::U1Irrep, spin::SU2Irrep}}
    CS2 = @NamedSector begin
        charge::U1Irrep
        spin::SU2Irrep
    end
    @test CS === CS2

    @test TensorKitSectors.type_repr(CS) == "@NamedSector{charge::Irrep[U₁], spin::Irrep[SU₂]}"
    @test Base.eval(Main, Meta.parse(TensorKitSectors.type_repr(CS))) === CS

    s1 = CS(U1Irrep(1), SU2Irrep(1 // 2))
    s2 = NamedSector(; charge = U1Irrep(1), spin = SU2Irrep(1 // 2))
    s3 = ⊠(; charge = U1Irrep(1), spin = SU2Irrep(1 // 2))
    s4 = CS((1, 1 // 2))
    s5 = CS(1, 1 // 2)
    @test s1 == s2 == s3 == s4 == s5
    @test typeof(s1) == typeof(s2) == typeof(s3) == typeof(s4) == typeof(s5) == CS

    @test s1.charge == U1Irrep(1)
    @test s1.spin == SU2Irrep(1 // 2)
    @test s1[:charge] == U1Irrep(1)
    @test s1[1] == U1Irrep(1)
    @test s1[2] == SU2Irrep(1 // 2)
    @test keys(s1) == (:charge, :spin)
    @test propertynames(s1) == (:sectors, :charge, :spin)

    @test Tuple(s1) == (U1Irrep(1), SU2Irrep(1 // 2))
    @test NamedTuple(s1) == (; charge = U1Irrep(1), spin = SU2Irrep(1 // 2))
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

@testset "Frobenius-Schur" begin
    @test frobenius_schur_phase(SU2Irrep(0)) == 1.0
    @test frobenius_schur_phase(SU2Irrep(1 // 2)) == -1.0
    @test frobenius_schur_phase(Z2Irrep(0)) == 1.0
    @test frobenius_schur_phase(Z2Irrep(1)) == 1.0
    @test frobenius_schur_phase(Z3Irrep(0)) == 1.0
    @test frobenius_schur_phase(Z3Irrep(1)) == 1.0
    @test frobenius_schur_phase(U1Irrep(0)) == 1.0
    @test frobenius_schur_phase(U1Irrep(1)) == 1.0

    @test frobenius_schur_indicator(SU2Irrep(0)) == 1.0
    @test frobenius_schur_indicator(SU2Irrep(1 // 2)) == -1.0
    @test frobenius_schur_indicator(Z2Irrep(0)) == 1.0
    @test frobenius_schur_indicator(Z2Irrep(1)) == 1.0
    @test frobenius_schur_indicator(Z3Irrep(0)) == 1.0
    @test frobenius_schur_indicator(Z3Irrep(1)) == 0.0
    @test frobenius_schur_indicator(U1Irrep(0)) == 1.0
    @test frobenius_schur_indicator(U1Irrep(1)) == 0.0
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

@testset "Topological spin" begin
    @test spin_top(IsingAnyon(:I)) == 0 // 1
    @test spin_top(IsingAnyon(:σ)) == 1 // 16
    @test spin_top(IsingAnyon(:ψ)) == 1 // 2
    @test spin_top(FibonacciAnyon(:τ)) == - 2 // 5
    @test spin_top(FibonacciAnyon(:I)) == 0 // 1
    @test spin_top(Z2Irrep(1)) == 0 // 1
    @test spin_top(TimeReversed{IsingAnyon}(:ψ)) == 1 // 2
    @test spin_top(TimeReversed{IsingAnyon}(:σ)) == - 1 // 16
    @test spin_top(TimeReversed{FibonacciAnyon}(:τ)) == 2 // 5
end

@testset "T vector" begin
    @test Tvector(Z2Irrep) ≈ [1, 1]
    @test Tvector(Z3Irrep) ≈ [1, 1, 1]
    @test Tvector(FermionParity) ≈ [1, -1]
    @test Tvector(Z2Irrep ⊠ Z3Irrep) ≈ ones(6)
    @test Tvector(FibonacciAnyon) ≈ [1, cispi(-4 / 5)]
    @test Tvector(IsingAnyon) ≈ [1, cispi(1 / 8), -1]
    @test Tvector(TimeReversed{IsingAnyon}) ≈ [1, cispi(-1 / 8), -1]
    @test Tvector(TimeReversed{FibonacciAnyon}) ≈ [1, cispi(4 / 5)]

    @test Tvector(FibonacciAnyon ⊠ FibonacciAnyon) ≈ kron(Tvector(FibonacciAnyon), Tvector(FibonacciAnyon))
    @test Tvector(FibonacciAnyon ⊠ TimeReversed{FibonacciAnyon}) ≈ kron(Tvector(FibonacciAnyon), Tvector(TimeReversed{FibonacciAnyon}))
    @test Tvector(FibonacciAnyon ⊠ IsingAnyon) ≈ kron(Tvector(FibonacciAnyon), Tvector(IsingAnyon))
    @test Tvector(IsingAnyon ⊠ TimeReversed{IsingAnyon}) ≈ kron(Tvector(IsingAnyon), Tvector(TimeReversed{IsingAnyon}))

    @test Tvector(A4Irrep) ≈ [1, 1, 1, 1]
end

@testset "Hopf link" begin
    @test hopflink(Z2Irrep(1), Z2Irrep(1)) ≈ 1
    @test hopflink(Z3Irrep(1), Z3Irrep(2)) ≈ 1
    @test hopflink(FermionParity(1), FermionParity(1)) ≈ 1
    @test hopflink(IsingAnyon(:ψ), IsingAnyon(:ψ)) ≈ 1
    @test hopflink(IsingAnyon(:σ), IsingAnyon(:σ)) ≈ 0
    @test hopflink(IsingAnyon(:σ), IsingAnyon(:ψ)) ≈ -sqrt(2)
    @test hopflink(FibonacciAnyon(:τ), FibonacciAnyon(:τ)) ≈ -1
end
@testset "S matrix" begin
    @test Smatrix(Z2Irrep) ≈ ones(2, 2)
    @test Smatrix(Z3Irrep) ≈ ones(3, 3)
    @test Smatrix(FermionParity) ≈ ones(2, 2)
    @test Smatrix(Z2Irrep ⊠ Z3Irrep) ≈ ones(6, 6)
    φ = (1 + sqrt(5)) / 2
    @test Smatrix(FibonacciAnyon) ≈ [1 φ; φ -1]
    @test Smatrix(IsingAnyon) ≈ [1 sqrt(2) 1; sqrt(2) 0 -sqrt(2); 1 -sqrt(2) 1]

    # S matrix is symmetric
    @test transpose(Smatrix(FibonacciAnyon ⊠ FibonacciAnyon)) ≈ Smatrix(FibonacciAnyon ⊠ FibonacciAnyon)
    @test transpose(Smatrix(FibonacciAnyon ⊠ IsingAnyon)) ≈ Smatrix(FibonacciAnyon ⊠ IsingAnyon)
    @test transpose(Smatrix(IsingAnyon ⊠ TimeReversed{IsingAnyon})) ≈ Smatrix(IsingAnyon ⊠ TimeReversed{IsingAnyon})
    @test transpose(Smatrix(TimeReversed{FibonacciAnyon} ⊠ IsingAnyon)) ≈ Smatrix(TimeReversed{FibonacciAnyon} ⊠ IsingAnyon)

    # S matrix is tensor producted when sectors are Deligne tensor producted
    @test Smatrix(FibonacciAnyon ⊠ FibonacciAnyon) ≈ kron(Smatrix(FibonacciAnyon), Smatrix(FibonacciAnyon))
    @test Smatrix(FibonacciAnyon ⊠ IsingAnyon) ≈ kron(Smatrix(FibonacciAnyon), Smatrix(IsingAnyon))
    @test Smatrix(IsingAnyon ⊠ TimeReversed{IsingAnyon}) ≈ kron(Smatrix(IsingAnyon), Smatrix(TimeReversed{IsingAnyon}))
    @test Smatrix(TimeReversed{FibonacciAnyon} ⊠ IsingAnyon) ≈ kron(Smatrix(TimeReversed{FibonacciAnyon}), Smatrix(IsingAnyon))
    @test Smatrix(IsingAnyon ⊠ IsingAnyon ⊠ IsingAnyon) ≈ kron(Smatrix(IsingAnyon), Smatrix(IsingAnyon), Smatrix(IsingAnyon))
end

@testset "Total quantum dimension" begin
    @test sqDim(Z2Irrep) ≈ 2
    @test sqDim(Z3Irrep) ≈ 3
    @test sqDim(FermionParity) ≈ 2
    @test sqDim(Z2Irrep ⊠ Z3Irrep) ≈ 6
    φ = (1 + sqrt(5)) / 2
    @test sqDim(FibonacciAnyon) ≈ 2 + φ
    @test sqDim(IsingAnyon) ≈ 4
    @test sqDim(FibonacciAnyon ⊠ FibonacciAnyon) ≈ (2 + φ)^2
    @test sqDim(IsingAnyon ⊠ TimeReversed{IsingAnyon}) ≈ 16
end

@testset "Topological central charge" begin
    @test c_top(IsingAnyon) == 1 // 2
    @test c_top(TimeReversed{IsingAnyon}) == - 1 // 2
    @test c_top(IsingAnyon ⊠ TimeReversed{IsingAnyon}) == 0 // 1
    @test c_top(FibonacciAnyon) == - 14 // 5
    @test c_top(TimeReversed{FibonacciAnyon}) == 14 // 5
    @test c_top(FibonacciAnyon ⊠ TimeReversed{FibonacciAnyon}) == 0 // 1
    @test c_top(FibonacciAnyon ⊠ FibonacciAnyon ⊠ FibonacciAnyon) == mod(- 3 * 14 // 5 + 4, 8) - 4
    @test c_top(FibonacciAnyon ⊠ FibonacciAnyon) == mod(2 * c_top(FibonacciAnyon) + 4, 8) - 4
    @test c_top(FibonacciAnyon ⊠ IsingAnyon) == mod(c_top(FibonacciAnyon) + c_top(IsingAnyon) + 4, 8) - 4
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
