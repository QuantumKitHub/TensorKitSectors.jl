```@meta
CollapsedDocStrings = true
```

# Optional Methods

The following methods have default implementations but can be overridden for performance or to provide additional functionality.

## Quantum Dimensions

The quantum dimension of a sector is a fundamental invariant that determines the behavior of fusions and braiding.
The default implementation extracts the dimension from [`Fsymbol`](@ref) via the quantum dimension formula.
For many common sectors, however, the dimension is known directly from representation theory.
In these cases, it can be beneficial to overload `dim` to bypass computing F-symbols, either for performance reasons or to enforce tighter output types.
For example, dimensions of irreducible representations of groups are always integers.

```math
d_a = \left| \frac{1}{(F_{a \bar{a} a}^a)^{1_a}_{_a1} } \right|
```

```@docs; canonical = false
dim
sqrtdim
invsqrtdim
```

## Frobenius-Schur Indicators

The Frobenius-Schur indicator and phase characterize the self-duality properties of sectors.
The indicator distinguishes real, complex, and quaternionic representations.
The phase is the category-theoretic version that appears in line bending operations.

```math
\kappa_a = \text{sign}\left( (F_{a \bar{a} a}^a)^{1_a}_{_a1} \right)
```

```@docs; canonical = false
frobenius_schur_indicator
frobenius_schur_phase
```

## Scalar Type

The scalar type declares what number type is used in the F-symbols and R-symbols for a given sector.
This is automatically inferred from the return types of those symbols, but can be explicitly declared for efficiency.

```@docs; canonical = false
sectorscalartype
```

## Topological Data Symbols

The [`Asymbol`](@ref), [`Bsymbol`](@ref) and [`twist`](@ref) are derived from F- and R-symbols but are often occurring combinations.
The A-symbol and B-symbol relate different ways of bending strands, while the twist is the topological spin (quantum dimension phase) of a sector.

The A-symbol ``A^{ab}_c`` relates splitting and fusion vertices:
```math
A^{ab}_c = \sqrt{\frac{d_c}{d_b}} \overline{\kappa_a} F^{\bar{a} a b}_b
```

The B-symbol ``B^{ab}_c`` relates splitting and fusion vertices:
```math
B^{ab}_c = \sqrt{\frac{d_c}{d_a}} F^{a b \bar{b}}_a R^{ab}_a F^{a \bar{b} b}_a
```

The twist ``\theta_a`` of a sector is the topological spin phase, computed as the trace of the R-matrix for braiding a sector with itself:
```math
\theta_a = \frac{1}{d_a} \sum_{b \in a \otimes a} d_b \text{tr}(R^{aa}_b)
```

```@docs; canonical = false
Asymbol
Bsymbol
twist
```

## Fusion Basis

The fusion tensor provides explicit matrix elements for the tensor product of representations.
It is a rank-4 array whose components are the Clebsch-Gordan coefficients for fusing sector `a` and `b` into `c`.
The fusion tensor is not uniquely determined by topological data alone, instead the topological data can be extracted from it when it is available.
However, there is not always a concrete representation of these tensors in terms of simple `Array` objects, and the exact representation is not important for TensorKit.jl.
For this reason, it is optional: TensorKit can work with the topological data alone.
Note however that they are still required whenever we want to convert symmetric tensors to and from dense arrays.

```@docs; canonical = false
fusiontensor
```
