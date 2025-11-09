using TensorOperations
using LinearAlgebra

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
    if BraidingStyle(I) isa HasBraiding
        R = @testinferred Rsymbol(s...)
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

@testsuite "fusion tensor and F-move" I -> begin
    (BraidingStyle(I) isa Bosonic && hasfusiontensor(I)) || return nothing
    for a in smallset(I), b in smallset(I), c in smallset(I)
        for e in ⊗(a, b), f in ⊗(b, c)
            for d in intersect(⊗(e, c), ⊗(a, f))
                X1 = fusiontensor(a, b, e)
                X2 = fusiontensor(e, c, d)
                Y1 = fusiontensor(b, c, f)
                Y2 = fusiontensor(a, f, d)
                @tensor f1[-1, -2, -3, -4] := conj(Y2[a, f, d, -4]) *
                    conj(Y1[b, c, f, -3]) * X1[a, b, e, -1] * X2[e, c, d, -2]
                if FusionStyle(I) isa MultiplicityFreeFusion
                    f2 = fill(Fsymbol(a, b, c, d, e, f) * dim(d), (1, 1, 1, 1))
                else
                    f2 = Fsymbol(a, b, c, d, e, f) * dim(d)
                end
                @test isapprox(f1, f2; atol = 1.0e-12, rtol = 1.0e-12)
            end
        end
    end
end

@testsuite "fusion tensor and F-move and R-move" I -> begin
    (BraidingStyle(I) isa Bosonic && hasfusiontensor(I)) || return nothing
    for a in smallset(I), b in smallset(I)
        for c in ⊗(a, b)
            X1 = permutedims(fusiontensor(a, b, c), (2, 1, 3, 4))
            X2 = fusiontensor(b, a, c)
            l = dim(a) * dim(b) * dim(c)
            R = LinearAlgebra.transpose(Rsymbol(a, b, c))
            sz = (l, convert(Int, Nsymbol(a, b, c)))
            @test reshape(X1, sz) ≈ reshape(X2, sz) * R
        end
    end
end

@testsuite "Orthogonality of fusiontensors" I -> begin
    hasfusiontensor(I) || return nothing
    for a in smallset(I), b in smallset(I)
        cs = vec(collect(a ⊗ b))
        CGCs = map(c -> reshape(fusiontensor(a, b, c), :, dim(c)), cs)
        M = map(Iterators.product(CGCs, CGCs)) do (cgc1, cgc2)
            return LinearAlgebra.norm(cgc1' * cgc2)
        end
        @test isapprox(M' * M, LinearAlgebra.Diagonal(dim.(cs)); atol = 1.0e-12)
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
