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
\psi \otimes \psi = I,\qquad
\sigma \otimes \psi = \psi \otimes \sigma = \sigma,\qquad
\sigma \otimes \sigma = I \oplus \psi.
```

Thus `FusionStyle(IsingAnyon) = SimpleFusion()`.
The quantum dimensions are

```math
d_I = 1,\qquad d_\psi = 1,\qquad d_\sigma = \sqrt{2}.
```

## Topological Data

The non-trivial associator for four ``\sigma`` anyons is

```math
F^{\sigma\sigma\sigma}_{\sigma} =
\frac{1}{\sqrt{2}}
\begin{pmatrix}
1 & 1\\
1 & -1
\end{pmatrix},
```

in the intermediate basis ``(I,\psi)``.
Additional signs involving the fermion ``\psi`` are encoded by [`Fsymbol`](@ref).

The braiding style is `Anyonic()`.
Important braiding phases are

```math
R^{\sigma\sigma}_I = e^{-\pi i/8},\qquad
R^{\sigma\sigma}_\psi = e^{3\pi i/8},\qquad
R^{\psi\psi}_I = -1.
```

There is no fusion tensor as the fusion category does not originate from a group or its representations.

## Iteration and basis conventions
`values(IsingAnyon)` iterates the labels in the order `:I`, `:σ`, `:ψ`.

## References

- [Ising anyon](https://en.wikipedia.org/wiki/Anyon#Ising_anyons)
- [Modular tensor category](https://en.wikipedia.org/wiki/Modular_tensor_category)
