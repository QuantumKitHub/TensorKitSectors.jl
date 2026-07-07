# Named Sectors: `NamedSector`

`NamedSector` represents the same Deligne tensor product as [`ProductSector`](@ref), but stores its component sectors in a *named* tuple so that they can be accessed by name rather than by position.
This is convenient when a product mixes several symmetries (e.g. a charge and a spin) and positional indexing would be error-prone.

## Sector type

```@docs; canonical = false
NamedSector
@NamedSector
```

A named sector is constructed with keyword arguments to `⊠`/`NamedSector`, or from a type built with the [`@NamedSector`](@ref) macro:

```julia
using TensorKitSectors

s = ⊠(; charge = U1Irrep(1), spin = SU2Irrep(1//2)) # (charge=1, spin=1/2)
s.charge # 1
s.spin   # 1/2

I = @NamedSector{charge::U1Irrep, spin::SU2Irrep}
I(U1Irrep(1), SU2Irrep(1//2)) # same as above
```


## Fusion Rules

`NamedSector` behaves identically to the underlying `ProductSector`: fusion, duality, and units act componentwise, and every piece of topological data is forwarded to the corresponding `ProductSector` of the same components.
Concretely, `NamedSector` converts to `ProductSector` (via `ProductSector(s)`) and delegates [`⊗`](@ref), [`Nsymbol`](@ref), [`dual`](@ref), [`unit`](@ref), and [`allunits`](@ref).
See the [`ProductSector`](@ref) page for the componentwise fusion rules.

## Topological Data

[`Fsymbol`](@ref), [`Rsymbol`](@ref), [`Asymbol`](@ref), [`Bsymbol`](@ref), [`fusiontensor`](@ref), and [`dim`](@ref) are all computed by forwarding to the underlying `ProductSector`.
Likewise, [`FusionStyle`](@ref), [`BraidingStyle`](@ref), `UnitStyle`, and the scalar-type traits are inherited from the product of the components.

## Iteration

`values(I)` for a `NamedSector` type `I` iterates the same Cartesian product as the corresponding `ProductSector`, wrapping each result back into the named type.
Indexing and `findindex` use the same ordering.

## References

- [Deligne tensor product](https://ncatlab.org/nlab/show/Deligne+tensor+product)
