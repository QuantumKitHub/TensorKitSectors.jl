# Fibonacci Anyons: `FibonacciAnyon`

`FibonacciAnyon` represents the Fibonacci modular fusion category.
It has a trivial anyon `:I` and a non-trivial anyon `:τ` (`:tau` is accepted as an ASCII constructor alias).

## Sector type

```@docs; canonical = false
FibonacciAnyon
```

Both sectors are self-dual, and the unit is `FibonacciAnyon(:I)`.

## Fusion Rules

The only non-trivial fusion rule is

```math
τ ⊗ τ = I ⊕ τ.
```

Thus `FusionStyle(FibonacciAnyon) = SimpleFusion()`.
The quantum dimensions are

```math
d_I = 1,\qquad d_τ = φ = \frac{1 + \sqrt{5}}{2}.
```

## Topological Data

The non-trivial associator appears when all external anyons are ``τ``:

```math
F^{τττ}_{τ} =
\begin{pmatrix}
φ^{-1} & φ^{-1/2}\\
φ^{-1/2} & -φ^{-1}
\end{pmatrix},
```

in the intermediate basis ``(I,τ)``.
All other allowed `Fsymbol` values are `1`.

The braiding style is `Anyonic()`.
For ``τ ⊗ τ``,

```math
R^{ττ}_I = e^{4π i/5},\qquad
R^{ττ}_τ = e^{-3π i/5}.
```

There is no fusion tensor as the fusion category does not originate from a group or its representations.

## Iteration and basis conventions
`values(FibonacciAnyon)` iterates the labels in the order `:I`, `:τ`.

## References

- [Fibonacci anyon](https://en.wikipedia.org/wiki/Anyon#Fibonacci_anyons)
- [Modular tensor category](https://en.wikipedia.org/wiki/Modular_tensor_category)
