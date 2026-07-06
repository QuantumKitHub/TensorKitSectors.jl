# Fermion Parity: `FermionParity`

`FermionParity` is the two-sector fermionic category ``SVect`` for even and odd parity.
Odd sectors pick up a minus sign when exchanged with another odd sector.

## Sector type

```@docs; canonical = false
FermionParity
```

The labels are `FermionParity(false)` for even parity and `FermionParity(true)` for odd parity.
They are displayed as `FermionParity(0)` and `FermionParity(1)`.
The unit is even parity, and every sector is self-dual.

## Fusion Rules

Fusion is addition modulo two:

```math
p ⊗ q = p ⊕ q.
```

The category has `FusionStyle(FermionParity) = UniqueFusion()`.

## Topological Data

All quantum dimensions are `1`, and [`Fsymbol`](@ref) is trivial on allowed channels.

The braiding style is `Fermionic()`.
The [`Rsymbol`](@ref) is `-1` only when both input sectors are odd:

```math
R^{11}_0 = -1.
```

The twists are

```math
θ_0 = 1,\qquad θ_1 = -1.
```

`fusiontensor` is defined for array construction, but it emits a warning because plain arrays with `FermionParity` labels do not preserve all fermionic signs.

## References

- [Fermion parity](https://en.wikipedia.org/wiki/Fermion_parity)
- [Super vector space](https://en.wikipedia.org/wiki/Super_vector_space)