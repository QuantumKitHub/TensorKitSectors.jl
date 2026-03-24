using TensorOperations
using LinearAlgebra
using TensorKitSectors: TensorKitSectors as TKS

@testsuite "Basic properties" I -> begin
    s = random_fusion(I, 2)
    sc = @testinferred(first(⊗(s...)))
    @test @testinferred(hash(sc)) == hash(deepcopy(sc))
    if UnitStyle(I) isa SimpleUnit
        @test @testinferred(unit(sc)) == @testinferred(unit(I))
    end
    @testinferred dual(sc)
    @testinferred dim(sc)
    @testinferred frobenius_schur_phase(sc)
    @testinferred frobenius_schur_indicator(sc)
    @testinferred(⊗(s..., s...))
    @testinferred Nsymbol(s..., sc)
    @testinferred Asymbol(s..., sc)
    B = @testinferred Bsymbol(s..., sc)
    s2 = random_fusion(I, 3)
    s2c = first(⊗(s2...))
    e, f = first(⊗(s2[1], s2[2])), first(⊗(s2[2], s2[3]))
    @testinferred Fsymbol(s2..., s2c, e, f)
    if BraidingStyle(I) isa HasBraiding
        @testinferred Rsymbol(s..., sc)
    end
end

@testsuite "Element and data types of topological data" I -> begin
    a, b, c = random_fusion(I, 3)
    d, e, f = first(⊗(a, b, c)), first(⊗(a, b)), first(⊗(b, c))
    F = Fsymbol(a, b, c, d, e, f)
    Ftype = @testinferred fusionscalartype(I)
    @test eltype(F) === Ftype
    @test F isa (FusionStyle(I) isa MultiplicityFreeFusion ? Ftype : AbstractArray{Ftype, 4})
    if BraidingStyle(I) isa HasBraiding
        a, b = random_fusion(I, 2)
        R = Rsymbol(a, b, first(⊗(a, b)))
        Rtype = @testinferred braidingscalartype(I)
        @test eltype(R) === Rtype
        if FusionStyle(I) isa MultiplicityFreeFusion
            @test typeof(R * F) <: @testinferred sectorscalartype(I)
            @test R isa Rtype
        else
            @test Base.promote_op(*, eltype(R), eltype(F)) <: @testinferred sectorscalartype(I)
            @test R isa AbstractArray{Rtype, 2}
        end
    else
        if FusionStyle(I) isa MultiplicityFreeFusion
            @test typeof(F) <: @testinferred sectorscalartype(I)
        else
            @test eltype(F) <: @testinferred sectorscalartype(I)
        end
    end
end

@testsuite "Show and parse" I -> begin
    # test in the parent module of where the test is defined
    _module = parentmodule(@__MODULE__)
    _eval(x) = Base.eval(_module, x)
    _sprint(x) = sprint(show, x; context = (:module => _module))
    @test _eval(Meta.parse(_sprint(I))) == I
    @test _eval(Meta.parse(TKS.type_repr(I))) == I
    for a in smallset(I)
        @test _eval(Meta.parse(_sprint(a))) == a
    end
end

@testsuite "Value iterator" I -> begin
    @test eltype(values(I)) == I
    sprev, rest = Iterators.peel(values(I))
    i = 1
    @test (@testinferred findindex(values(I), sprev)) == i
    for s in rest
        i += 1
        @test !isless(s, sprev)
        @test s == @testinferred(values(I)[i])
        @test findindex(values(I), s) == i
        sprev = s
        i >= 10 && break
    end
    for s in smallset(I)
        @test (@testinferred values(I)[findindex(values(I), s)]) == s
    end
end

@testsuite "Invalid fusion channel checks" I -> begin
    for a in smallset(I), b in smallset(I)
        can_fuse(a, b) || @test_throws ArgumentError Nsymbol(a, b, first(smallset(I)))
    end
    for b in smallset(I), c in smallset(I)
        can_fuse(b, c) || @test_throws ArgumentError Nsymbol(b, c, first(smallset(I)))
    end
end

@testsuite "Shapes of topological data" I -> begin
    # shape of data from multiplicities
    for a in smallset(I), b in smallset(I)
        can_fuse(a, b) || continue
        for c in smallset(I)
            can_fuse(b, c) || continue
            for e in ⊗(a, b), f in ⊗(b, c) # at this point can_fuse(e, c) and can_fuse(a, f) are always true
                Nabe, Nbcf = Nsymbol(a, b, e), Nsymbol(b, c, f)
                for d in ⊗(a, b, c)
                    F_size = FusionStyle(I) isa MultiplicityFreeFusion ? () : (Nabe, Nsymbol(e, c, d), Nbcf, Nsymbol(a, f, d))
                    @test size(Fsymbol(a, b, c, d, e, f)) == F_size
                end
            end
        end
    end

    for a in smallset(I), b in smallset(I), c in smallset(I)
        if BraidingStyle(I) isa HasBraiding
            R_size = FusionStyle(I) isa MultiplicityFreeFusion ? () : (Nsymbol(a, b, c), Nsymbol(b, a, c))
            @test size(Rsymbol(a, b, c)) == R_size
        end
    end
end

# https://arxiv.org/pdf/2507.07023v2#=&page=9
@testsuite "Fusion ring properties" I -> begin
    for a in smallset(I), b in smallset(I)
        for c in smallset(I), d in smallset(I) # associativity
            L = sum(e -> Nsymbol(a, b, e) * Nsymbol(e, c, d), intersect(⊗(a, b), ⊗(d, dual(c))); init = zero(Int))
            R = sum(f -> Nsymbol(b, c, f) * Nsymbol(a, f, d), intersect(⊗(b, c), ⊗(dual(a), d)); init = zero(Int))
            @test L == R
        end
        for c in ⊗(a, b) # Frobenius reciprocity
            @test Nsymbol(a, b, c) == Nsymbol(c, dual(b), a) == Nsymbol(dual(c), a, dual(b)) ==
                Nsymbol(dual(b), dual(a), dual(c)) == Nsymbol(b, dual(c), dual(a)) == Nsymbol(dual(a), c, b)
            if BraidingStyle(I) isa HasBraiding
                @test Nsymbol(a, b, c) == Nsymbol(b, a, c)
            end
        end
    end
end

@testsuite "Fusion and dimensions" I -> begin
    for a in smallset(I), b in smallset(I)
        can_fuse(a, b) || continue
        da = dim(a)
        db = dim(b)
        dc = sum(c -> dim(c) * Nsymbol(a, b, c), a ⊗ b)
        @test da * db ≈ dc # needs to be ≈ because of anyons
    end
end

@testsuite "Fusion tensor and Fsymbol" I -> begin
    hasfusiontensor(I) || return nothing
    for a in smallset(I), b in smallset(I), c in smallset(I)
        for e in ⊗(a, b), f in ⊗(b, c)
            for d in intersect(⊗(e, c), ⊗(a, f))
                F1 = Fsymbol(a, b, c, d, e, f)
                F2 = TKS.Fsymbol_from_fusiontensor(a, b, c, d, e, f)
                @test F1 ≈ F2 atol = 1.0e-12 rtol = 1.0e-12
            end
        end
    end
end

@testsuite "Fsymbol and Asymbol" I -> begin
    for a in smallset(I), b in smallset(I)
        for c in ⊗(a, b)
            A1 = Asymbol(a, b, c)
            A2 = TKS.Asymbol_from_Fsymbol(a, b, c)
            @test A1 ≈ A2 atol = 1.0e-12 rtol = 1.0e-12
        end
    end
end

@testsuite "Fsymbol and Bsymbol" I -> begin
    for a in smallset(I), b in smallset(I)
        for c in ⊗(a, b)
            B1 = Bsymbol(a, b, c)
            B2 = TKS.Bsymbol_from_Fsymbol(a, b, c)
            @test B1 ≈ B2 atol = 1.0e-12 rtol = 1.0e-12
        end
    end
end

@testsuite "Fusion tensor and Asymbol" I -> begin
    (BraidingStyle(I) isa Bosonic && hasfusiontensor(I)) || return nothing
    for a in smallset(I), b in smallset(I)
        for c in ⊗(a, b)
            A1 = Asymbol(a, b, c)
            A2 = TKS.Asymbol_from_fusiontensor(a, b, c)
            @test A1 ≈ A2 atol = 1.0e-12 rtol = 1.0e-12
        end
    end
end

@testsuite "Fusion tensor and Bsymbol" I -> begin
    (BraidingStyle(I) isa Bosonic && hasfusiontensor(I)) || return nothing
    for a in smallset(I), b in smallset(I)
        for c in ⊗(a, b)
            B1 = Bsymbol(a, b, c)
            B2 = TKS.Bsymbol_from_fusiontensor(a, b, c)
            @test B1 ≈ B2 atol = 1.0e-12 rtol = 1.0e-12
        end
    end
end

@testsuite "Fsymbol and dim" I -> begin
    for a in smallset(I)
        @test dim(a) ≈ TKS.dim_from_Fsymbol(a) atol = 1.0e-12 rtol = 1.0e-12
    end
end

@testsuite "Fsymbol and frobenius_schur_phase" I -> begin
    for a in smallset(I)
        @test frobenius_schur_phase(a) ≈ TKS.frobenius_schur_phase_from_Fsymbol(a) atol = 1.0e-12 rtol = 1.0e-12
    end
end

@testsuite "Fusion tensor and Rsymbol" I -> begin
    (BraidingStyle(I) isa Bosonic && hasfusiontensor(I)) || return nothing
    for a in smallset(I), b in smallset(I)
        for c in ⊗(a, b)
            R1 = Rsymbol(a, b, c)
            R2 = TKS.Rsymbol_from_fusiontensor(a, b, c)
            @test R1 ≈ R2 atol = 1.0e-12 rtol = 1.0e-12
        end
    end
end

@testsuite "Rsymbol and twist" I -> begin
    BraidingStyle(I) isa HasBraiding || return nothing
    for a in smallset(I)
        @test twist(a) ≈ TKS.twist_from_Rsymbol(a) atol = 1.0e-12 rtol = 1.0e-12
    end
end

@testsuite "Orthogonality of fusiontensors" I -> begin
    hasfusiontensor(I) || return nothing
    for a in smallset(I), b in smallset(I)
        cs = vec(collect(a ⊗ b))
        cgcs = map(c -> fusiontensor(a, b, c), cs)
        for (c, cgc) in zip(cs, cgcs), (c′, cgc′) in zip(cs, cgcs)
            for μ in 1:Nsymbol(a, b, c), ν in 1:Nsymbol(a, b, c′)
                @tensor overlap[mc mc'] := conj(cgc[:, :, :, μ][ma mb mc]) *
                    cgc′[:, :, :, ν][ma mb mc']
                if μ == ν && c == c′
                    @test isapprox(overlap, LinearAlgebra.I; atol = 1.0e-12)
                else
                    @test isapprox(LinearAlgebra.norm(overlap), 0; atol = 1.0e-12)
                end
            end
        end
    end
end

@testsuite "Unitarity of F-move" I -> begin
    for a in smallset(I), b in smallset(I)
        can_fuse(a, b) || continue
        for c in smallset(I)
            can_fuse(b, c) || continue
            @test F_unitarity_test(a, b, c; atol = 1.0e-12, rtol = 1.0e-12)
        end
    end
end

@testsuite "Unitarity of R-move" I -> begin
    BraidingStyle(I) isa HasBraiding || return nothing
    for a in smallset(I), b in smallset(I)
        @test R_unitarity_test(a, b; atol = 1.0e-12, rtol = 1.0e-12)
    end
end

@testsuite "Triangle equation" I -> begin
    for a in smallset(I), b in smallset(I)
        can_fuse(a, b) || continue
        @test triangle_equation(a, b; atol = 1.0e-12, rtol = 1.0e-12)
    end
end

@testsuite "Pentagon equation" I -> begin
    for a in smallset(I), b in smallset(I)
        can_fuse(a, b) || continue
        for c in smallset(I)
            can_fuse(b, c) || continue
            for d in smallset(I)
                can_fuse(c, d) || continue
                @test pentagon_equation(a, b, c, d; atol = 1.0e-12, rtol = 1.0e-12)
            end
        end
    end
end

# https://ncatlab.org/nlab/files/DelaneyModularTensorCategories.pdf#page=9
@testsuite "Hexagon equation" I -> begin
    BraidingStyle(I) isa HasBraiding || return nothing
    for a in smallset(I), b in smallset(I), c in smallset(I)
        @test hexagon_equation(a, b, c; atol = 1.0e-12, rtol = 1.0e-12)
    end
end

# https://quantumkithub.github.io/TensorKit.jl/stable/appendix/categories/#Braidings-and-twists
@testsuite "Ribbon condition" I -> begin
    BraidingStyle(I) isa HasBraiding || return nothing
    for a in smallset(I), b in smallset(I)
        for c in ⊗(a, b)
            R1 = Rsymbol(a, b, c)
            R2 = Rsymbol(b, a, c)
            θa, θb, θc = twist.((a, b, c))
            factor = R1 * θa * θb * R2
            if FusionStyle(I) isa GenericFusion
                @test isapprox(θc * LinearAlgebra.I, factor; atol = 1.0e-12, rtol = 1.0e-12)
            else
                @test isapprox(θc, factor; atol = 1.0e-12, rtol = 1.0e-12)
            end
        end
    end
end

@testsuite "Braiding self-duality condition" I -> begin
    BraidingStyle(I) isa HasBraiding || return nothing
    for a in smallset(I)
        if a == dual(a)
            factor = twist(a) * frobenius_schur_phase(a) * Rsymbol(a, a, unit(a))
            @test factor ≈ one(factor) atol = 1.0e-12 rtol = 1.0e-12
        end
    end
end

@testsuite "Symmetric braiding condition" I -> begin
    BraidingStyle(I) isa SymmetricBraiding || return nothing
    isfermionic = BraidingStyle(I) isa Fermionic
    for a in smallset(I)
        θa = twist(a)
        oneT = one(θa)
        @test isapprox(θa, oneT; atol = 1.0e-12, rtol = 1.0e-12) ||
            (isfermionic && isapprox(θa, -oneT; atol = 1.0e-12, rtol = 1.0e-12))
        for b in smallset(I)
            for c in ⊗(a, b)
                RR = Rsymbol(a, b, c) * Rsymbol(b, a, c)
                @test RR ≈ one(RR) atol = 1.0e-12 rtol = 1.0e-12
            end
        end
    end
end

# https://quantumkithub.github.io/TensorKit.jl/stable/man/fusiontrees/#Manipulations-on-a-fusion-tree
@testsuite "Artin braid equality" I -> begin
    BraidingStyle(I) isa HasBraiding || return nothing
    for a in smallset(I), b in smallset(I), d in smallset(I)
        for f in ⊗(d, a)
            Rdaf, Radf = Rsymbol(d, a, f), Rsymbol(a, d, f)
            for c in ⊗(a, b)
                for e in intersect(⊗(c, d), ⊗(f, b))
                    Rcde, Rdce = Rsymbol(c, d, e), Rsymbol(d, c, e)
                    Fdabefc = Fsymbol(d, a, b, e, f, c)
                    if FusionStyle(I) isa MultiplicityFreeFusion
                        RFR1 = Rcde * conj(Fdabefc) * conj(Rdaf)
                        RFR2 = conj(Rdce) * conj(Fdabefc) * Radf
                    else
                        @tensor RFR1[ν, μ, λ, σ] := Rcde[ν, ρ] * conj(Fdabefc[κ, λ, μ, ρ]) * conj(Rdaf[σ, κ])
                        @tensor RFR2[ν, μ, λ, σ] := conj(Rdce[ν, ρ]) * conj(Fdabefc[κ, λ, μ, ρ]) * Radf[σ, κ]
                    end
                    FRF1, FRF2 = zero(RFR1), zero(RFR2)
                    for g in ⊗(d, b)
                        Fabdecg = Fsymbol(a, b, d, e, c, g)
                        Fadbefg = Fsymbol(a, d, b, e, f, g)
                        Rbdg, Rdbg = Rsymbol(b, d, g), Rsymbol(d, b, g)
                        if FusionStyle(I) isa MultiplicityFreeFusion
                            FRF1 += Fabdecg * Rbdg * conj(Fadbefg)
                            FRF2 += conj(Fabdecg) * conj(Rdbg) * Fadbefg
                        else
                            @tensor FRF1[ν, μ, β, α] += Fabdecg[μ, ν, κ, λ] * Rbdg[κ, θ] * conj(Fadbefg[α, β, θ, λ])
                            @tensor FRF2[ν, μ, β, α] += conj(Fabdecg[μ, ν, κ, λ]) * conj(Rdbg[κ, θ]) * Fadbefg[α, β, θ, λ]
                        end
                    end
                    @test isapprox(RFR1, FRF1; atol = 1.0e-12, rtol = 1.0e-12)
                    @test isapprox(RFR2, FRF2; atol = 1.0e-12, rtol = 1.0e-12)
                end
            end
        end
    end
end
