# Time-Reversed Sectors: `TimeReversed`

`TimeReversed{I}` represents the time-reversed, or conjugate-braided, version of a sector type `I`.
It keeps the same objects, fusion rules, dimensions, and associators, but reverses the braiding.

## Sector type

```@docs; canonical = false
TimeReversed
timereversed
```

Construct a wrapped sector with `TimeReversed(a)` or `TimeReversed{I}(a)`.
The helper `timereversed(a)` avoids unnecessary wrappers for symmetric braiding categories and unwraps an already time-reversed sector:

```julia
using TensorKitSectors

a = IsingAnyon(:σ)
timereversed(a)               # TimeReversed{IsingAnyon}(:σ)
timereversed(timereversed(a)) # IsingAnyon(:σ)
```

`TimeReversed` is not defined for sectors with `NoBraiding()`.

## Fusion Rules

Fusion is inherited from the original sector type.
The fusion multiplicities, fusion style, units, duals, and quantum dimensions are unchanged:

```math
N_{\overline{c}}^{\overline{a}\,\overline{b}} = N_c^{ab},\qquad
d_{\overline{a}} = d_a,\qquad
\overline{a}^{\,*} = \overline{a^*}.
```

`values(TimeReversed{I})` follows the same order as `values(I)`, with every element wrapped.

## Topological Data

The [`Fsymbol`](@ref), [`Asymbol`](@ref), and [`Bsymbol`](@ref) are inherited directly from the original sector.
The [`Rsymbol`](@ref) is replaced by its adjoint:

```math
R_{\overline{c}}^{\overline{a}\,\overline{b}} = \left(R_c^{ab}\right)^\dagger.
```

For scalar braiding phases this is complex conjugation, so anyonic spins change sign:

```math
θ_{\overline{a}} = \overline{θ_a}.
```

For example, in the conventions of this package,

```math
θ_σ = e^{π i/8}
\quad⟹\quad
θ_{\overline{σ}} = e^{-π i/8}
```

for the Ising `σ` anyon.

The braiding style and scalar types are reported as those of the original sector type.
For bosonic or fermionic symmetric categories, time reversal is physically trivial; `timereversed(a)` therefore returns `a` directly.

## References

- [Time reversal symmetry](https://en.wikipedia.org/wiki/T-symmetry)
