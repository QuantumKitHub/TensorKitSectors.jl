# Fermion Number: `FermionNumber`

`FermionNumber` is a convenience sector for integer ``U(1)`` charge together with the matching fermion parity.
It is implemented as `U1Irrep ⊠ FermionParity`, with odd charge assigned odd parity.

## Sector type

```@docs; canonical = false
FermionNumber
```

Construct `FermionNumber(n)` from an integer charge `n`.
The underlying product label is

```math
n \mapsto (U(1)\text{ charge } n,\; n \bmod 2).
```

The unit is `FermionNumber(0)`, and duality negates the ``U(1)`` charge while preserving the parity constraint.

## Fusion Rules

Fusion adds charges:

```math
n_1 \otimes n_2 = n_1 + n_2.
```

The fermion parity component is then fixed automatically by the resulting charge.

## Topological data

All sectors have quantum dimension `1`.

The ``U(1)`` part has bosonic representation-category data, while the parity part supplies the fermionic exchange sign.
Thus exchanging two odd-charge sectors contributes a minus sign.

Because this is a product sector, its fusion and topological data is inherited componentwise from [`U1Irrep`](@ref) and [`FermionParity`](@ref).

## References

- [Fermion number](https://en.wikipedia.org/wiki/Fermion_number)