using TensorOperations
using LinearAlgebra
using TensorKitSectors: TensorKitSectors as TKS

@testsuite "Basic properties" I -> begin
    s = (randsector(I), randsector(I), randsector(I))
    @test Base.eval(Main, Meta.parse(sprint(show, I))) == I
    @test Base.eval(Main, Meta.parse(TensorKitSectors.type_repr(I))) == I
    @test Base.eval(Main, Meta.parse(sprint(show, s[1]))) == s[1]
    @test @testinferred(hash(s[1])) == hash(deepcopy(s[1]))
    @test @testinferred(unit(s[1])) == @testinferred(unit(I))
    @testinferred dual(s[1])
    @testinferred dim(s[1])
    @testinferred frobenius_schur_phase(s[1])
    @testinferred frobenius_schur_indicator(s[1])
    @testinferred Nsymbol(s...)
    @testinferred Asymbol(s...)
    B = @testinferred Bsymbol(s...)
    F = @testinferred Fsymbol(s..., s...)
    @test eltype(F) === @testinferred fusionscalartype(I)
    if BraidingStyle(I) isa HasBraiding
        R = @testinferred Rsymbol(s...)
        @test eltype(R) === @testinferred braidingscalartype(I)
        if FusionStyle(I) === SimpleFusion()
            @test typeof(R * F) <: @testinferred sectorscalartype(I)
        else
            @test Base.promote_op(*, eltype(R), eltype(F)) <: @testinferred sectorscalartype(I)
        end
    else
        if FusionStyle(I) === SimpleFusion()
            @test typeof(F) <: @testinferred sectorscalartype(I)
        else
            @test eltype(F) <: @testinferred sectorscalartype(I)
        end
    end
    @testinferred(s[1] ⊗ s[2])
    @testinferred(⊗(s..., s...))
end

@testsuite "Value iterator" I -> begin
    @test eltype(values(I)) == I
    sprev = unit(I)
    for (i, s) in enumerate(values(I))
        @test !isless(s, sprev)
        @test s == @testinferred(values(I)[i])
        @test findindex(values(I), s) == i
        sprev = s
        i >= 10 && break
    end
    @test unit(I) == first(values(I))
    @test length(allunits(I)) == 1
    @test (@testinferred findindex(values(I), unit(I))) == 1
    for s in smallset(I)
        @test (@testinferred values(I)[findindex(values(I), s)]) == s
    end
end

@testsuite "Fusion and dimensions" I -> begin
    for a in smallset(I), b in smallset(I)
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
        cgcs = map(c -> Array(fusiontensor(a, b, c)), cs) # make dense for sparse arrays to work
        for (c, cgc) in zip(cs, cgcs), (c′, cgc′) in zip(cs, cgcs)
            for μ in 1:Nsymbol(a, b, c), ν in 1:Nsymbol(a, b, c′)
                @tensor overlap[mc mc'] := conj(view(cgc, :, :, :, μ)[ma mb mc]) *
                    view(cgc′, :, :, :, ν)[ma mb mc']
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
    for a in smallset(I), b in smallset(I), c in smallset(I)
        for d in ⊗(a, b, c)
            es = collect(intersect(⊗(a, b), map(dual, ⊗(c, dual(d)))))
            fs = collect(intersect(⊗(b, c), map(dual, ⊗(dual(d), a))))
            if FusionStyle(I) isa MultiplicityFreeFusion
                @test length(es) == length(fs)
                F = [Fsymbol(a, b, c, d, e, f) for e in es, f in fs]
            else
                Fblocks = Vector{Any}()
                for e in es, f in fs
                    Fs = Fsymbol(a, b, c, d, e, f)
                    push!(Fblocks, reshape(Fs, (size(Fs, 1) * size(Fs, 2), size(Fs, 3) * size(Fs, 4))))
                end
                F = hvcat(length(fs), Fblocks...)
            end
            @test isapprox(F' * F, one(F); atol = 1.0e-12, rtol = 1.0e-12)
        end
    end
end

@testsuite "Triangle equation" I -> begin
    for a in smallset(I), b in smallset(I)
        @test triangle_equation(a, b; atol = 1.0e-12, rtol = 1.0e-12)
    end
end

@testsuite "Pentagon equation" I -> begin
    for a in smallset(I), b in smallset(I), c in smallset(I), d in smallset(I)
        @test pentagon_equation(a, b, c, d; atol = 1.0e-12, rtol = 1.0e-12)
    end
end

@testsuite "Hexagon equation" I -> begin
    BraidingStyle(I) isa HasBraiding || return nothing
    for a in smallset(I), b in smallset(I), c in smallset(I)
        @test hexagon_equation(a, b, c; atol = 1.0e-12, rtol = 1.0e-12)
    end
end
