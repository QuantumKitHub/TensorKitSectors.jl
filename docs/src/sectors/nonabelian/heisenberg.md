# Representation of the finite Heisenberg Group ``H_N``: `HeisenbergIrrep{N}`

`HeisenbergIrrep{N}` represents irreducible representations (irreps) of the finite Heisenberg group ``H_N`` as a `Sector`.
``H_N`` is the group of ``3× 3`` upper triangular matrices which is identified with the group ``ℤ_N ⋊ ℤ_N^2``; it is only defined here for prime `N`.
This page documents how the type behaves in code (construction, iteration, fusion, and access to topological data).

## Sector type

```@docs; canonical = false
HeisenbergIrrep
```

Every irrep is labeled by a single integer `n::Int8`, but the meaning of `n` depends on its sign:

- If `0 <= n < N^2`, the irrep is the one-dimensional character ``χ_{a,b}`` with `a, b = divrem(n, N)`.
  There are ``N^2`` such characters, one for each pair `(a, b) ∈ 0:N-1`.
- If `-N < n < 0`, the irrep is the ``N``-dimensional Schrödinger representation ``π_k`` with `k = -n`, `k ∈ 1:N-1`.

In total there are ``N^2 + N - 1`` irreps (the trivial ``k = 0`` Schrödinger representation is excluded, since it is reducible and decomposes into the ``N^2`` characters).
The alias `Heis3Irrep = HeisenbergIrrep{3}` is provided for the smallest (and most commonly used) case.

The trivial sector is `HeisenbergIrrep{N}(0)`, i.e. ``χ_{0,0}``.
Duals are

```math
χ_{a,b}^* = χ_{-a \bmod N,\,-b \bmod N}, \qquad π_k^* = π_{N-k}.
```

## Fusion Rules

Characters fuse among themselves as the abelian group ``ℤ_N × ℤ_N``:

```math
χ_{a,b} ⊗ χ_{a',b'} = χ_{(a+a')\bmod N,\,(b+b')\bmod N}.
```

Fusing a character with a Schrödinger irrep leaves the latter unchanged, regardless of the character's labels:

```math
χ_{a,b} ⊗ π_k = π_k ⊗ χ_{a,b} = π_k.
```

This reflects that ``π_k`` is the unique irrep with a given nontrivial central character, so twisting it by any character of the full group produces an isomorphic representation.

Two Schrödinger irreps fuse according to whether their labels cancel modulo `N`:

```math
π_k ⊗ π_{k'} =
\begin{cases}
N ⋅ π_{(k+k') \bmod N}, & k + k' ≢ 0 \pmod N, \\[4pt]
\displaystyle\bigoplus_{a,b=0}^{N-1} χ_{a,b}, & k + k' ≡ 0 \pmod N.
\end{cases}
```

In the first case the ``N``-fold multiplicity accounts for the dimension count ``N × N = N^2`` (``N`` copies of the ``N``-dimensional ``π_{(k+k')\bmod N}``); in the second, the ``N^2``-dimensional product decomposes into all ``N^2`` one-dimensional characters, each with multiplicity one.

Because of this multiplicity-``N`` channel, `FusionStyle(HeisenbergIrrep) = GenericFusion()`.
The [`Nsymbol`](@ref) returns the corresponding integer multiplicities directly from the case distinction above.

Quantum dimensions are the ordinary representation space dimensions:

```math
d_{χ_{a,b}} = 1, \qquad d_{π_k} = N.
```

## Topological Data

`BraidingStyle(HeisenbergIrrep) = Bosonic()`: braiding is symmetric, and `braidingscalartype(HeisenbergIrrep) = ComplexF64` since the Clebsch-Gordan coefficients are generally complex.

The [`Fsymbol`](@ref) and [`Rsymbol`](@ref) are both computed generically from the chosen [`fusiontensor`](@ref) basis (via `Fsymbol_from_fusiontensor` and `Rsymbol_from_fusiontensor`), see below.
On the single-channel sectors (``χ ⊗ χ → χ`` and ``χ ⊗ π → π``) the R-symbol evaluates to the trivial phase `1`.
On the multiplicity-``N`` channel ``π_k ⊗ π_{k'} → π_{(k+k')\bmod N}`` and on the channels ``π_k ⊗ π_{-k} → χ_{a,b}``, the R-symbol carries nontrivial phases built out of `N`-th roots of unity.

## Fusion Tensor and Basis Conventions

`fusiontensor(a, b, c)` returns a rank-4 array of size ``d_a × d_b × d_c × N_c^{ab}``.
The one-dimensional (character) fusion tensors are trivial scalars.
For the Schrödinger irreps, the fusion tensor is built out of the ``N``-th root of unity ``ω = e^{2π i/N}``, following the reference below:

- ``π_k ⊗ π_{k'} → π_{(k+k')\bmod N}`` (``k+k' ≢ 0``): a permutation tensor (with a phase for the case where one of the factors is a character, see below).
- ``π_k ⊗ π_{-k} → χ_{a,b}``: normalized with a factor ``1/\sqrt{N}``, with phase ``ω^{-a i}`` along the surviving diagonal.
- ``χ_{a,b} ⊗ π_k → π_k`` (and the symmetric ``π_k ⊗ χ_{a,b} → π_k``): a permutation tensor dressed with the phase ``ω^{a m}``.

The representation can be checked against the generators `X` and `Z` of ``H_N`` acting as ``X|i\rangle = |i+1\rangle`` and ``Z|i\rangle = ω^{ki}|i\rangle`` on ``π_k``, and as ``X ↦ ω^a``, ``Z ↦ ω^b`` on ``χ_{a,b}``.
For every allowed fusion channel, the fusion tensor acts as an intertwiner between the tensor product representation and the fused representation:

```jldoctest heisenberg_intertwiner; output = false
using TensorKitSectors
using LinearAlgebra: kron, Diagonal
using Test: @test

N = 3
ω = cispi(2 / N)

function X(s::HeisenbergIrrep{N}) where {N}
    s.n < 0 || return hcat(ω^div(s.n, N))
    M = zeros(ComplexF64, N, N)
    for i in 0:(N - 1)
        M[mod(i + 1, N) + 1, i + 1] = 1
    end
    return M
end
Z(s::HeisenbergIrrep{N}) where {N} = s.n < 0 ? Diagonal([ω^(-s.n * i) for i in 0:(N - 1)]) : hcat(ω^mod(s.n, N))

for a in values(HeisenbergIrrep{N}), b in values(HeisenbergIrrep{N})
    for c in a ⊗ b
        C = fusiontensor(a, b, c)
        for μ in 1:Nsymbol(a, b, c)
            Cmat = reshape(view(C, :, :, :, μ), dim(a) * dim(b), dim(c))
            @test Cmat' * kron(X(b), X(a)) * Cmat ≈ X(c)
            @test Cmat' * kron(Z(b), Z(a)) * Cmat ≈ Z(c)
        end
    end
end

# output
```

## Iteration

Iterating over `values(HeisenbergIrrep{N})` first yields the ``N^2`` characters ``χ_{0,0}, χ_{0,1}, …`` (in order of increasing `n`), followed by the Schrödinger irreps ``π_1, …, π_{N-1}``.

## References

For background on the group and the representation basis, see:

- [Heisenberg group](https://en.wikipedia.org/wiki/Heisenberg_group#Heisenberg_group_modulo_an_odd_prime_p)
- [Bacon, D., How a Clebsch-Gordan Transform Helps to Solve the Heisenberg Hidden Subgroup Problem](https://arxiv.org/abs/quant-ph/0612107)
