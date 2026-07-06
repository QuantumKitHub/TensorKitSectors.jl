# Product Sectors: `ProductSector`

`ProductSector` represents the Deligne tensor product of sector categories.
It is the standard way to combine independent symmetry or anyon labels, such as charge and parity, or two different anyon theories. For the bosonic group or representation categories, this Deligne product corresponds to the ordinary product of groups or direct product of representations.

## Sector type

```@docs; canonical = false
ProductSector
⊠
deligneproduct
```

A product sector stores its component sectors in a tuple.
The recommended constructor is the Deligne product operator `⊠` or `boxtimes`.
The more verbose option is `deligneproduct`:

```julia
using TensorKitSectors

a = Z2Irrep(1) ⊠ FibonacciAnyon(:τ) # (Irrep[ℤ₂](1) ⊠ FibonacciAnyon(:τ))
deligneproduct(Z2Irrep(1), FibonacciAnyon(:τ)) # same as above
```

The same operator works on sector types:

```julia
using TensorKitSectors

I = Z2Irrep ⊠ FibonacciAnyon
a = I(1, :τ)
```

Products are flattened, and `Trivial` factors are removed:

```math
(a \boxtimes b) \boxtimes c = a \boxtimes b \boxtimes c,\qquad
a \boxtimes \mathbf{1} = a.
```

For representation categories, type-level products are displayed as irreps of the product group.
For example, `Z2Irrep ⊠ Z3Irrep` behaves as `Irrep[ℤ₂ × ℤ₃]`.
Product sectors are used to define several convenience aliases, such as [`FermionNumber`](@ref) and [`FermionSpin`](@ref).

## Fusion Rules

Fusion is componentwise.
For product sectors

```math
a = (a_1,\ldots,a_n),\qquad b = (b_1,\ldots,b_n),
```

the fusion outputs are

```math
a \otimes b = 
\bigoplus_{c_i \in a_i \otimes b_i} N^{ab}_c
(c_1,\ldots,c_n).
```

The fusion multiplicities multiply:

```math
N^{ab}_c = \prod_i N^{a_i b_i}_{c_i}.
```

The fusion style is the combined style of the components.
In particular, a product of unique-fusion sectors remains unique, while adding a simple or generic component promotes the product accordingly.

Duals, units, and quantum dimensions are also componentwise:

```math
(a_1,\ldots,a_n)^* = (a_1^*,\ldots,a_n^*),\qquad
d_{(a_1,\ldots,a_n)} = \prod_i d_{a_i}.
```

If a component category has multiple units, `allunits` returns the product of the component unit sets.

## Topological Data

The topological data is the tensor product of the component data.
For multiplicity-free channels this reduces to ordinary multiplication:

```math
F_{\boxtimes_i a_i,\boxtimes_i b_i,\boxtimes_i c_i}^{\boxtimes_i d_i}
= \prod_i F_{a_i b_i c_i}^{d_i},
\qquad
R_{\boxtimes_i a_i,\boxtimes_i b_i}^{\boxtimes_i c_i}
= \prod_i R_{a_i b_i}^{c_i}.
```

When multiplicity spaces are present, the implementation forms the corresponding Kronecker products so that the product multiplicity basis is ordered componentwise.
The same convention is used for `Fsymbol`, `Rsymbol`, `Asymbol`, `Bsymbol`, and `fusiontensor`.

Braiding style and scalar types are combined from the component categories.

For products containing [`FermionParity`](@ref), `fermionparity` is the XOR of the component parities.

## Iteration

`values(I1 ⊠ I2)` iterates over the Cartesian product of `values(I1)` and `values(I2)`.
Indexing and `findindex` use the same componentwise ordering.
For infinite component sectors, the product iterator can also be infinite.

## References

- [Deligne tensor product](https://ncatlab.org/nlab/show/Deligne+tensor+product)