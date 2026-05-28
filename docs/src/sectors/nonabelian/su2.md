# ``SU(2)`` Representations: `SU2Irrep`

`SU2Irrep` represents irreducible representations of the compact group ``SU(2)`` as a `Sector`.
This page documents how the type behaves in code (construction, iteration, fusion, and access to topological data).

## Sector type

The irreducible representations of ``SU(2)`` are labeled by non-negative half-integers (e.g. ``0``, ``\frac{1}{2}``, ``1``, ``\frac{3}{2}``, …).

```@docs; canonical = false
SU2Irrep
```

Fusing two irreps together leads to a direct sum of irreps:

```math
a \otimes b = \bigoplus_{c_j = |j_a - j_b|}^{j_a + j_b} c
```

Since each output appears only once, we have `FusionStyle(SU2Irrep) = SimpleFusion()`.
The [`Nsymbol`](@ref) returns a `Bool` that checks the triangle inequality:

```math
N_c^{ab} = |j_a - j_b| \leq j_c \leq j_a + j_b \land j_a + j_b + j_c \in \mathbb{N}
```

Each irrep has dimension ``d = 2j + 1``.

The [`Fsymbol`](@ref) is computed from [Wigner ``6j``](https://en.wikipedia.org/wiki/6-j_symbol) (Racah-``W``) symbols as

```math
\left(F_{abc}^d\right)_e^f = (-1)^{j_a + j_b + j_c + j_e} \sqrt{d_e d_f}  \begin{Bmatrix}
j_a & j_b & j_d \\
j_c & j_e & j_f
\end{Bmatrix}
```

The `BraidingStyle` is `Bosonic`, and [`Rsymbol`](@ref) is ``\pm 1`` for allowed fusion channels based on the parity of ``j_a + j_b - j_c``.

## Fusion Tensor and Basis Conventions

`fusiontensor(a, b, c)` returns Clebsch–Gordan coefficients as a rank‑4 array of size ``d_a \times d_b \times d_c \times 1``.
We can label the basis by using ``\ket{j, m}``, where the magnetic quantum number takes on values ``m \in \{j, j-1, \ldots, -j\}`` (in that order).

Each irrep acts on the standard basis $\lvert j, m \rangle$ where $m = j, j-1, \dots, -j$.
In that basis, the generators of ``\mathfrak{su}(2)`` are represented by the usual angular momentum operators:

```math
J_z \lvert j, m \rangle = m \lvert j, m \rangle, \quad
J_\pm \lvert j, m \rangle = \sqrt{(j \mp m)(j \pm m + 1)}\, \lvert j, m \pm 1 \rangle.
```

```jldoctest generators; output = false
using TensorKitSectors

function generators(a::SU2Irrep)
    Jp = zeros(dim(a), dim(a))
    Jm = zeros(dim(a), dim(a))
    Jz = zeros(dim(a), dim(a))
    
    for row in axes(Jp, 1), col in axes(Jp, 2)
        m = a.j - col + 1
        if row == col
            Jz[row, col] = m
        elseif row + 1 == col
            Jp[row, col] = sqrt((a.j - m) * (a.j + m + 1))
        elseif row == col + 1
            Jm[row, col] = sqrt((a.j + m) * (a.j - m + 1))
        end
    end

    return Jp, Jm, Jz
end

a = SU2Irrep(1)
Jp, Jm, Jz = generators(a)

# output
([0.0 1.4142135623730951 0.0; 0.0 0.0 1.4142135623730951; 0.0 0.0 0.0], [0.0 0.0 0.0; 1.4142135623730951 0.0 0.0; 0.0 1.4142135623730951 0.0], [1.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 -1.0])
```

The tensor product representation for ``a \otimes b`` is given in terms of the representations ``a`` and ``b`` as:

```math
\mathbf{J}^{(a \otimes b)} = \mathbf{J}^{(a)} \otimes \mathbf{1}^{(b)} + \mathbf{1}^{(a)} \otimes \mathbf{J}^{(b)}
```

The [`fusiontensor`](@ref) supplies the change‑of‑basis coefficients ``C^{J M}_{j_a m_a, j_b m_b}`` that map the uncoupled product basis to the coupled basis:

```math
\lvert J, M \rangle = \sum_{m_a, m_b} C^{J M}_{j_a m_a, j_b m_b}
\lvert j_a, m_a \rangle \otimes \lvert j_b, m_b \rangle.
```

In particular, this block-diagonalizes the generators, and we must have that for every ``\mathbf{J}``, the following holds:

```math
\left(\mathbf{J}^{(a)} \otimes \mathbf{1}^{(b)} + \mathbf{1}^{(a)} \otimes \mathbf{J}^{(b)}\right) \cdot C^c_{ab} =
    C^c_{ab} \cdot \mathbf{J}^{(c)}
```

```jldoctest generators; output = false
using TensorOperations: @tensor
using Test: @test

a = SU2Irrep(1)
b = SU2Irrep(1)

for c in a ⊗ b
    CGC = dropdims(fusiontensor(a, b, c); dims = 4) # drop trivial multiplicity dimension
    for (ga, gb, gc) in zip(generators(a), generators(b), generators(c))
        @tensor lhs[a b; c] := ga[a; a'] * CGC[a' b; c] + gb[b; b'] * CGC[a b'; c]
        @tensor rhs[a b; c] := CGC[a b; c'] * gc[c'; c]
        @test isapprox(lhs, rhs)
    end
end

# output

```

## References

For a quick refresher on the group structure and representation theory (without turning this page into a math course), the following references are useful:

- [Special unitary group](https://en.wikipedia.org/wiki/Special_unitary_group)
- [SU(2)](https://en.wikipedia.org/wiki/SU(2))
- [Representation theory of SU(2)](https://en.wikipedia.org/wiki/Representation_theory_of_SU(2))
- [Clebsch–Gordan coefficients](https://en.wikipedia.org/wiki/Clebsch%E2%80%93Gordan_coefficients)
- [Wigner 6j symbol](https://en.wikipedia.org/wiki/6-j_symbol)
