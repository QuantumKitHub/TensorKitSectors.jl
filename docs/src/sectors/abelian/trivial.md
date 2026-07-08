```@meta
DocTestSetup = quote
    using TensorKitSectors
end
```

# Trivial Sector: `Trivial`

`Trivial` is the trivial sector for ordinary vector spaces.
It is the trivial representation of the trivial group, equivalently `Rep[ℤ₁]`.

## Sector type

```@docs; canonical = false
Trivial
```

There is only one label:

```jldoctest
julia> Trivial()
Trivial()
```

## Fusion Rules

Fusion is uniquely trivial:

```math
1 ⊗ 1 = 1.
```

Therefore `FusionStyle(Trivial) = UniqueFusion()` and `Nsymbol(Trivial(), Trivial(), Trivial()) == true`.

## Topological Data

The quantum dimension is `1`.
`fusiontensor` is a `1 × 1 × 1 × 1` array containing `1`.
The associator and braiding are both trivial:

```math
F = 1,\qquad R = 1.
```

The braiding style is `Bosonic()`.

## References

- [Trivial group](https://en.wikipedia.org/wiki/Trivial_group)
