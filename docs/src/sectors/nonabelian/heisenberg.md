# Representation of the finite Heisenberg Group ``H_N``: `HeisenbergIrrep{N}`

`HeisenbergIrrep{N}` represents irreducible representations (irreps) of the finite Heisenberg group ``H_N`` as a `Sector`.
``H_N`` is the group of ``3\times 3`` upper triangular matrices which is identified with the group ``\mathbb{Z}_N \rtimes \mathbb{Z}_N^2``; it is only defined here for prime `N`.
This page documents how the type behaves in code (construction, iteration, fusion, and access to topological data).

## Sector type

```@docs; canonical = false
HeisenbergIrrep
```

Every irrep is labeled by a single integer `n::Int8`, but the meaning of `n` depends on its sign:

- If `0 <= n < N^2`, the irrep is the one-dimensional character ``\chi_{a,b}`` with `a, b = divrem(n, N)`.
  There are ``N^2`` such characters, one for each pair `(a, b) ∈ 0:N-1`.
- If `-N < n < 0`, the irrep is the ``N``-dimensional Schrödinger representation ``\pi_k`` with `k = -n`, `k ∈ 1:N-1`.

In total there are ``N^2 + N - 1`` irreps (the trivial ``k = 0`` Schrödinger representation is excluded, since it is reducible and decomposes into the ``N^2`` characters).
The alias `Heis3Irrep = HeisenbergIrrep{3}` is provided for the smallest (and most commonly used) case.

The trivial sector is `HeisenbergIrrep{N}(0)`, i.e. ``\chi_{0,0}``.
Duals are

```math
\chi_{a,b}^* = \chi_{-a \bmod N,\,-b \bmod N}, \qquad \pi_k^* = \pi_{N-k}.
```

## Fusion Rules

Characters fuse among themselves as the abelian group ``\mathbb{Z}_N \times \mathbb{Z}_N``:

```math
\chi_{a,b} \otimes \chi_{a',b'} = \chi_{(a+a')\bmod N,\,(b+b')\bmod N}.
```

Fusing a character with a Schrödinger irrep leaves the latter unchanged, regardless of the character's labels:

```math
\chi_{a,b} \otimes \pi_k = \pi_k \otimes \chi_{a,b} = \pi_k.
```

This reflects that ``\pi_k`` is the unique irrep with a given nontrivial central character, so twisting it by any character of the full group produces an isomorphic representation.

Two Schrödinger irreps fuse according to whether their labels cancel modulo `N`:

```math
\pi_k \otimes \pi_{k'} =
\begin{cases}
N \cdot \pi_{(k+k') \bmod N}, & k + k' \not\equiv 0 \pmod N, \\[4pt]
\displaystyle\bigoplus_{a,b=0}^{N-1} \chi_{a,b}, & k + k' \equiv 0 \pmod N.
\end{cases}
```

In the first case the ``N``-fold multiplicity accounts for the dimension count ``N \times N = N^2`` (``N`` copies of the ``N``-dimensional ``\pi_{(k+k')\bmod N}``); in the second, the ``N^2``-dimensional product decomposes into all ``N^2`` one-dimensional characters, each with multiplicity one.

Because of this multiplicity-``N`` channel, `FusionStyle(HeisenbergIrrep) = GenericFusion()`.
The [`Nsymbol`](@ref) returns the corresponding integer multiplicities directly from the case distinction above.

Quantum dimensions are the ordinary representation space dimensions:

```math
d_{\chi_{a,b}} = 1, \qquad d_{\pi_k} = N.
```

## Topological Data

`BraidingStyle(HeisenbergIrrep) = Bosonic()`: braiding is symmetric, and `braidingscalartype(HeisenbergIrrep) = ComplexF64` since the Clebsch-Gordan coefficients are generally complex.

The [`Fsymbol`](@ref) and [`Rsymbol`](@ref) are both computed generically from the chosen [`fusiontensor`](@ref) basis (via `Fsymbol_from_fusiontensor` and `Rsymbol_from_fusiontensor`), see below.
On the single-channel sectors (``\chi \otimes \chi \to \chi`` and ``\chi \otimes \pi \to \pi``) the R-symbol evaluates to the trivial phase `1`.
On the multiplicity-``N`` channel ``\pi_k \otimes \pi_{k'} \to \pi_{(k+k')\bmod N}`` and on the channels ``\pi_k \otimes \pi_{-k} \to \chi_{a,b}``, the R-symbol carries nontrivial phases built out of `N`-th roots of unity.

## Fusion Tensor and Basis Conventions

`fusiontensor(a, b, c)` returns a rank-4 array of size ``d_a \times d_b \times d_c \times N_c^{ab}``.
The one-dimensional (character) fusion tensors are trivial scalars.
For the Schrödinger irreps, the fusion tensor is built out of the ``N``-th root of unity ``\omega = e^{2\pi i/N}``, following the reference below:

- ``\pi_k \otimes \pi_{k'} \to \pi_{(k+k')\bmod N}`` (``k+k' \not\equiv 0``): a permutation tensor (with a phase for the case where one of the factors is a character, see below).
- ``\pi_k \otimes \pi_{-k} \to \chi_{a,b}``: normalized with a factor ``1/\sqrt{N}``, with phase ``\omega^{-a i}`` along the surviving diagonal.
- ``\chi_{a,b} \otimes \pi_k \to \pi_k`` (and the symmetric ``\pi_k \otimes \chi_{a,b} \to \pi_k``): a permutation tensor dressed with the phase ``\omega^{a m}``.

The representation can be checked against the generators `X` and `Z` of ``H_N`` acting as ``X|i\rangle = |i+1\rangle`` and ``Z|i\rangle = \omega^{ki}|i\rangle`` on ``\pi_k``, and as ``X \mapsto \omega^a``, ``Z \mapsto \omega^b`` on ``\chi_{a,b}``.
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

Iterating over `values(HeisenbergIrrep{N})` first yields the ``N^2`` characters ``\chi_{0,0}, \chi_{0,1}, \ldots`` (in order of increasing `n`), followed by the Schrödinger irreps ``\pi_1, \ldots, \pi_{N-1}``.

## References

For background on the group and the representation basis, see:

- [Heisenberg group](https://en.wikipedia.org/wiki/Heisenberg_group#Heisenberg_group_modulo_an_odd_prime_p)
- [Bacon, D., How a Clebsch-Gordan Transform Helps to Solve the Heisenberg Hidden Subgroup Problem](https://arxiv.org/abs/quant-ph/0612107)
