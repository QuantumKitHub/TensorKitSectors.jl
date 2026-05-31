# Fermion Spin: `FermionSpin`

`FermionSpin` is a convenience sector for ``SU(2)`` spin together with the matching fermion parity.
It is implemented as `SU2Irrep ⊠ FermionParity`, with half-integer spin assigned odd parity.

## Sector type

```@docs; canonical = false
FermionSpin
```

Construct `FermionSpin(j)` from an integer or half-integer spin.
The parity is odd exactly when `2j` is odd:

```math
j \mapsto (SU(2)\text{ spin }j,\; 2j \bmod 2).
```

The unit is `FermionSpin(0)`, and all sectors are self-dual because the ``SU(2)`` irreps and parity sectors are self-dual.

## Fusion Rules

The ``SU(2)`` component fuses by angular-momentum addition:

```math
j_1 \otimes j_2 =
\bigoplus_{j = |j_1-j_2|}^{j_1+j_2} j.
```

The parity component is fixed by the output spin.
Quantum dimensions are inherited from `SU2Irrep`:

```math
d_j = 2j + 1.
```

The parity component supplies fermionic signs, so exchanging two half-integer-spin sectors has the odd-parity sign in addition to the `SU2Irrep` representation-category data.

Because this is a product sector, most fusion and topological data is inherited componentwise from [`SU2Irrep`](@ref) and [`FermionParity`](@ref).

## References

- [Spin (physics)](https://en.wikipedia.org/wiki/Spin_(physics))
- [SU(2)](https://en.wikipedia.org/wiki/SU(2))
- [Fermion parity](https://en.wikipedia.org/wiki/Fermion_parity)
