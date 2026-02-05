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

@testsuite "Shapes and data types of topological data" I -> begin
    # shape of data from fusion style
    Ftype = fusionscalartype(I)
    TF_actual = Core.Compiler.return_type(Fsymbol, NTuple{6, I})
    TF_imposed = FusionStyle(I) isa MultiplicityFreeFusion ? Ftype : Array{Ftype, 4} # won't work for sparse arrays (CategoryData)
    @test TF_actual == TF_imposed

    if BraidingStyle(I) isa HasBraiding
        Rtype = braidingscalartype(I)
        TR_actual = Core.Compiler.return_type(Rsymbol, NTuple{3, I})
        TR_imposed = if FusionStyle(I) isa MultiplicityFreeFusion
            Rtype
        else
            I <: TimeReversed ? LinearAlgebra.Adjoint{Rtype, Array{Rtype, 2}} : Array{Rtype, 2} # same here
        end
        @test TR_actual == TR_imposed
    end

    # shape of data from multiplicities
    for a in smallset(I), b in smallset(I), c in smallset(I)
        for e in ⊗(a, b), f in ⊗(b, c)
            for d in intersect(⊗(e, c), ⊗(a, f))
                F_size = FusionStyle(I) isa MultiplicityFreeFusion ? () : (Nsymbol(a, b, e), Nsymbol(e, c, d), Nsymbol(b, c, f), Nsymbol(a, f, d))
                @test size(Fsymbol(a, b, c, d, e, f)) == F_size
            end
        end

        if BraidingStyle(I) isa HasBraiding
            Nabc = Nsymbol(a, b, c)
            R_size = FusionStyle(I) isa MultiplicityFreeFusion ? () : (Nabc, Nabc)
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
        end
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

@testsuite "Fusion tensor and F-move" I -> begin
    hasfusiontensor(I) || return nothing
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

@testsuite "Fusion tensor and R-move" I -> begin
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
        cgcs = map(c -> fusiontensor(a, b, c), cs)
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
