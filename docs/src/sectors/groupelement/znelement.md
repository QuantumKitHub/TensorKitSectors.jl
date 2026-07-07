# Cyclic Group Elements: `ZNElement{N,p}`

`ZNElement{N,p}` represents the elements of the cyclic group ``ℤ_N`` directly as a `Sector`, as opposed to [`ZNIrrep`](@ref) which represents its (one-dimensional) irreps.
Fusion is given by the group multiplication itself, and the associator can carry a nontrivial phase determined by a level-`p` 3-cocycle on ``ℤ_N``.
For special values of `p` this also yields a consistent braiding, turning `ZNElement{N,p}` into a simple model of abelian anyons.
As a fusion category, this is identified with ``\text{Vec}_{ℤ_N}^{ω_p}``, the category of ``ℤ_N``-graded vector spaces with associator the 3-cocycle ``ω_p``.

## Group elements as sectors

```@docs; canonical = false
AbstractGroupElement
```

`ZNElement{N,p}` is currently the only concrete subtype of `AbstractGroupElement{Group}`.

## Sector type

```@docs; canonical = false
ZNElement
```

Elements are labeled by a single integer `n`, stored as `Int8` and only meaningful modulo `N` (which must satisfy `N < 64`).
The concrete type is obtained through `GroupElement[ℤ{N}, p]`, or through the aliases `Z2Element{p}`, `Z3Element{p}`, and `Z4Element{p}` for `N = 2, 3, 4`.
If `p` is omitted, it defaults to `0`, i.e. the trivial cocycle.

The unit is `ZNElement{N,p}(0)`, and duals are given by group inversion:

```math
a^* = a^{-1} = -a \bmod N.
```

## Fusion Rules

Fusion is literally the group operation, so there is a single fusion channel with multiplicity one:

```math
a ⊗ b = (a + b) \bmod N,
```

and `FusionStyle(ZNElement{N,p}) = UniqueFusion()` for every `N` and `p`.
All quantum dimensions are `1`.

## Associativity: the 3-cocycle and F-symbol

The [`Fsymbol`](@ref) is not assumed trivial: it is given by the 3-cocycle

```math
ω_p(a,b,c) = \exp\left(\frac{2π i\, p\, a\,\big(b + c - [b+c]_N\big)}{N^2}\right), \qquad [x]_N := x \bmod N,
```

which is a representative of a class in ``H^3(ℤ_N, U(1)) ≅ ℤ_N``, indexed by the level `p ∈ 0:N-1`.
Concretely,

```math
F^{abc}_{e,d,f} = ω_p(a,b,c) \quad \text{when } e = a+b,\ d = e+c,\ f = b + c \pmod N,
```

and zero otherwise (since fusion is unique, all the other symbols such as [`Asymbol`](@ref), [`Bsymbol`](@ref) and `frobenius_schur_phase` are likewise built directly out of `ω_p`).
For `p = 0` the cocycle is trivial and `ZNElement{N,0}` reduces to the ordinary ``ℤ_N``-graded vector spaces.
The [`pentagon_equation`](@ref) holds for every `N` and every `p ∈ 0:N-1`, i.e. every `ZNElement{N,p}` is a consistent fusion category.

## Topological Data (Braiding)

Braiding is only defined in two cases:

- `p = 0`: `BraidingStyle(ZNElement{N,0}) = Bosonic()`, and the [`Rsymbol`](@ref) is trivial (`Rsymbol(a,b,c) = Nsymbol(a,b,c)`).
- `2p ≡ N`: `BraidingStyle(ZNElement{N,p}) = Anyonic()`, and

  `````math
  R^{ab}_{a+b} = \exp\left(\frac{2\pi i\, p\, a\, b}{N^2}\right),
  `````

  giving abelian anyonic braiding statistics for the ``ℤ_N`` charges. For example, `N = 4, p = 2` (used in the example below) yields topological spins ``θ_a = \exp(iπ a^2 / 4)``.

For any other value of `p`, `BraidingStyle(ZNElement{N,p}) = NoBraiding()`.

## Iteration

`values(ZNElement{N,p})` iterates the labels `0, 1, …, N-1` in increasing order, exactly as for [`ZNIrrep`](@ref).

## Examples

```jldoctest znelement_consistency
using TensorKitSectors

a, b = Z3Element(1), Z3Element(2)
only(a ⊗ b) == unit(typeof(a))

# output
true
```

```jldoctest znelement_consistency
I = Z4Element{2} # 2p == N, so this admits a braiding
BraidingStyle(I)

# output
Anyonic()
```

## References

- [Cyclic group](https://en.wikipedia.org/wiki/Cyclic_group)
- [Group cohomology](https://en.wikipedia.org/wiki/Group_cohomology)