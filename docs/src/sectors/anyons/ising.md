# Ising Anyons: `IsingAnyon`

`IsingAnyon` represents the Ising modular fusion category.
It has labels `:I`, `:σ`, and `:ψ`; the ASCII aliases `:sigma` and `:psi` are accepted by the constructor.

## Sector type

```@docs; canonical = false
IsingAnyon
```

All three sectors are self-dual, and the unit is `IsingAnyon(:I)`.

## Fusion Rules

The non-trivial fusion rules are

```math
ψ ⊗ ψ = I,\qquad
σ ⊗ ψ = ψ ⊗ σ = σ,\qquad
σ ⊗ σ = I ⊕ ψ.
```

Thus `FusionStyle(IsingAnyon) = SimpleFusion()`.
The quantum dimensions are

```math
d_I = 1,\qquad d_ψ = 1,\qquad d_σ = \sqrt{2}.
```

## Topological Data

The non-trivial associator for four ``σ`` anyons is

```math
F^{σσσ}_{σ} =
\frac{1}{\sqrt{2}}
\begin{pmatrix}
1 & 1\\
1 & -1
\end{pmatrix},
```

in the intermediate basis ``(I,ψ)``.
Additional signs involving the fermion ``ψ`` are encoded by [`Fsymbol`](@ref).

The braiding style is `Anyonic()`.
Important braiding phases are

```math
R^{σσ}_I = e^{-π i/8},\qquad
R^{σσ}_ψ = e^{3π i/8},\qquad
R^{ψψ}_I = -1.
```

There is no fusion tensor as the fusion category does not originate from a group or its representations.

## Iteration and basis conventions
`values(IsingAnyon)` iterates the labels in the order `:I`, `:σ`, `:ψ`.

## References

- [Ising anyon](https://en.wikipedia.org/wiki/Anyon#Ising_anyons)
- [Modular tensor category](https://en.wikipedia.org/wiki/Modular_tensor_category)
