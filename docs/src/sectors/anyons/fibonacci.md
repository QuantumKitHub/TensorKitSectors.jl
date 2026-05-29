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
\tau \otimes \tau = I \oplus \tau.
```

Thus `FusionStyle(FibonacciAnyon) = SimpleFusion()`.
The quantum dimensions are

```math
d_I = 1,\qquad d_\tau = \varphi = \frac{1 + \sqrt{5}}{2}.
```

## Topological Data

The non-trivial associator appears when all external anyons are ``\tau``:

```math
F^{\tau\tau\tau}_{\tau} =
\begin{pmatrix}
\varphi^{-1} & \varphi^{-1/2}\\
\varphi^{-1/2} & -\varphi^{-1}
\end{pmatrix},
```

in the intermediate basis ``(I,\tau)``.
All other allowed `Fsymbol` values are `1`.

The braiding style is `Anyonic()`.
For ``\tau \otimes \tau``,

```math
R^{\tau\tau}_I = e^{4\pi i/5},\qquad
R^{\tau\tau}_\tau = e^{-3\pi i/5}.
```

There is no fusion tensor as the fusion category does not originate from a group or its representations.

## Iteration and basis conventions
`values(FibonacciAnyon)` iterates the labels in the order `:I`, `:τ`.

## References

- [Fibonacci anyon](https://en.wikipedia.org/wiki/Anyon#Fibonacci_anyons)
- [Modular tensor category](https://en.wikipedia.org/wiki/Modular_tensor_category)
