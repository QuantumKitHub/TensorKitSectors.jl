```@meta
CollapsedDocStrings = true
```

# Optional Methods

Most of the following methods have default implementations derived from the required methods, but can be overridden for performance or to provide additional functionality.
The exception is [`fusiontensor`](@ref), which has no default: it is genuinely optional and, when needed, must be supplied explicitly (see the Fusion Basis section below).

## Quantum Dimensions

The quantum dimension of a sector is a fundamental invariant that determines the behavior of fusions and braiding.
The default implementation extracts the dimension from [`Fsymbol`](@ref) via the quantum dimension formula.

```math
d_a = \left| \frac{1}{(F_{a \bar{a} a}^a)^{1_a}_{_a1} } \right|
```

For many common sectors, however, the dimension is known directly from representation theory.
In these cases, it can be beneficial to overload [`dim`](@ref) to bypass computing F-symbols, either for performance reasons or to enforce tighter output types.
For example, dimensions of irreducible representations of groups are always integers.

```@docs; canonical = false
dim
sqrtdim
invsqrtdim
```

## Frobenius-Schur Indicators

The Frobenius-Schur indicator and phase characterize the self-duality properties of sectors.
The Frobenius-Schur phase ``κ_a`` is obtained from the F-symbol as

```math
κ_a = \text{sign}\left( (F_{a \bar{a} a}^a)^{1_a}_{_a1} \right)
```

and is the category-theoretic quantity that appears in line-bending operations.
The Frobenius-Schur indicator ``ν_a`` coincides with ``κ_a`` for self-dual sectors (``a = \bar{a}``) and is ``0`` otherwise; it distinguishes real, complex, and quaternionic representations.

```@docs; canonical = false
frobenius_schur_indicator
frobenius_schur_phase
```

## Scalar Type

Various utility functions exist for determining what number type is used in various parts of the topological data.

```@docs; canonical = false
fusionscalartype
braidingscalartype
dimscalartype
sectorscalartype
```

!!! note
    While there is a fallback definition that tries to determine the result from computing the functions on the unit sector,
    it is often a good idea to define this method explicitly to avoid depending on compiler heuristics to constant-fold these calls.

## Topological Data Symbols

The [`Asymbol`](@ref), [`Bsymbol`](@ref) and [`twist`](@ref) are derived from F- and R-symbols but are often occurring combinations.
The A-symbol and B-symbol relate different ways of bending strands, while the twist is the topological spin (quantum dimension phase) of a sector.

The A-symbol ``A^{ab}_c`` relates splitting and fusion vertices:
```math
A^{ab}_c = \sqrt{\frac{d_a d_b}{d_c}} \, \overline{κ_a (F_{\bar{a} a b}^b)^1_c}
```

The B-symbol ``B^{ab}_c`` relates splitting and fusion vertices:
```math
B^{ab}_c = \sqrt{\frac{d_a d_b}{d_c}} (F_{a b \bar{b}}^a)^c_1
```

The twist ``θ_a`` of a sector is the topological spin phase, computed as the trace of the R-matrix for braiding a sector with itself:
```math
θ_a = \frac{1}{d_a} \sum_{b ∈ a ⊗ a} d_b \text{tr}(R^{aa}_b)
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
