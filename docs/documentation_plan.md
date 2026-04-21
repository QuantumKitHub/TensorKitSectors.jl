# Documentation Improvement Plan for TensorKitSectors.jl

This document outlines the plan to address issues #55, #56, and #57 regarding documentation improvements.

## Overview

The goal is to comprehensively document:
1. The `Sector` interface for implementing new sector types
2. The conventions (labels, fusion rules, basis choices) for all concrete sector implementations

## Issue #56: Sector Interface Documentation

**Goal:** Create clear documentation explaining what the `Sector` interface is and how to implement it.

### Create new page: `docs/src/interface.md`

#### Section 1: Introduction
- What is a Sector? (simple objects in unitary fusion categories)
- Connection to representation theory and fusion categories
- How sectors relate to symmetric tensors and graded vector spaces
- Brief mention that sectors implement the topological data of fusion categories

#### Section 2: Required Methods

Group methods by functionality with clear documentation, starting with iteration:

**Iteration over Sector Values:**
- `Base.iterate(::SectorValues{I}, state...)`
- `Base.IteratorSize(::Type{SectorValues{I}})` - HasLength, IsInfinite, or SizeUnknown
- If `HasLength()`: also implement `length`, `getindex`, and `findindex`
- Explain `SectorValues{I}` as a singleton type for iterating over sector values
- This is fundamental: sectors must be iterable to be useful

**Identity and Duality:**
- `unit(::Type{I})` or `allunits(::Type{I})` - the unit/identity object(s)
- `dual(a::I)` - the dual/conjugate object
- `leftunit(a::I)` and `rightunit(a::I)` (for multifusion categories)
- Note that `one(a)` and `conj(a)` work as aliases

**Fusion Structure:**
- `⊗(a::I, b::I)` - returns iterable of unique fusion outputs (each c appears once)
- `Nsymbol(a::I, b::I, c::I)` - fusion multiplicities (number of times c appears)
- `FusionStyle(::Type{I})` - UniqueFusion, SimpleFusion, or GenericFusion
- Explain relationship between FusionStyle and form of Nsymbol (Bool vs Integer)

**Topological Data:**
- `Fsymbol(a, b, c, d, e, f)` - F-symbol/recoupling coefficients (associator)
  - Scalar for UniqueFusion/SimpleFusion
  - Rank-4 array for GenericFusion
- `Rsymbol(a, b, c)` - R-symbol/braiding
  - Scalar for UniqueFusion/SimpleFusion
  - Matrix for GenericFusion
- `BraidingStyle(::Type{I})` - NoBraiding, Bosonic, Fermionic, Anyonic
- Explain what each braiding style means

**Utility Methods:**
- `Base.isless(a::I, b::I)` - canonical ordering (important for representing spaces)
- `Base.hash(a::I, h::UInt)` - hashing (sectors used as dictionary keys)

For each method: provide signature, purpose, requirements, and examples.

#### Section 3: Optional Methods

List with explanations of when/why to implement:
- `dim(a::I)` - quantum dimension (extracted from Fsymbol by default)
- `sqrtdim(a::I)`, `invsqrtdim(a::I)` - cached square roots for efficiency
- `frobenius_schur_indicator(a::I)` - returns 1, 0, or -1 (from Fsymbol by default)
- `frobenius_schur_phase(a::I)` - returns ±1 (from Fsymbol by default)
- `sectorscalartype(::Type{I})` - scalar type of F/R symbols (auto-inferred by default)
- `Bsymbol(a, b, c)` - B-symbol (computed from F-symbol by default)
- `twist(a)` - twist/ribbon element (computed from R-symbol by default)
- `fusiontensor(a, b, c)` - explicit basis for fusion spaces (important for applications)
- `Base.getindex(::SectorValues{I}, i::Int)` - access i-th sector value by index
- `findindex(::SectorValues{I}, c::I)` - find index of a sector value
- `Base.convert(::Type{I}, ...)` - convenience for constructing sectors

Explain that these have default implementations but can be overridden for performance.

#### Section 4: Traits and Styles

**FusionStyle:**
- `UniqueFusion()` - single fusion output (e.g., ℤₙ, U₁)
  - Nsymbol returns Bool
  - F-symbols and R-symbols are scalars
- `SimpleFusion()` - multiple outputs, multiplicity-free (e.g., SU₂)
  - Nsymbol returns Bool (0 or 1)
  - F-symbols and R-symbols are scalars
- `GenericFusion()` - multiple outputs with multiplicities > 1 (e.g., SU₃)
  - Nsymbol returns Integer
  - F-symbols are rank-4 arrays, R-symbols are matrices
  - Requires multiplicity labels μ

Abstract supertypes: `MultipleFusion` (includes Simple and Generic), `MultiplicityFreeFusion` (includes Unique and Simple)

**BraidingStyle:**
- `NoBraiding()` - no braiding structure (planar diagrams only)
- `Bosonic()` - symmetric braiding, R^{ab}_c R^{ba}_c = 1, all twists = +1
- `Fermionic()` - symmetric braiding with non-trivial twists (±1)
- `Anyonic()` - general non-symmetric braiding (braid group, not permutation group)

`SymmetricBraiding` is abstract supertype of Bosonic and Fermionic

**UnitStyle:**
- `SimpleUnit()` - unique unit element (standard case)
- `GenericUnit()` - multiple units (multifusion categories only)

Note: Most users will only need SimpleUnit

#### Section 5: Implementation Guidelines

**Helper Type: SectorProductIterator**
- Instead of returning tuple/array from ⊗, can return `SectorProductIterator(a, b)`
- Implement iterator interface for this type
- Type-stable and memory efficient
- Pretty printing shows as `a ⊗ b`

**Key Considerations:**
- Prefer `isbitstype` types for sectors when possible (performance)
- Ensure consistency: pentagon equation for F-symbols, hexagon equation for R-symbols
- Type stability: use SectorProductIterator for varying number of fusion outputs
- For `IteratorSize(SectorValues{I})`:
  - `HasLength()` for small finite sets → efficient GradedSpace representation
  - `IsInfinite()` or `SizeUnknown()` for large/infinite sets
- Scalar type of topological data affects tensor operations (real vs complex)

**Shape of Topological Data (determined by FusionStyle):**
- **UniqueFusion / SimpleFusion:** Nsymbol, F-symbols, and R-symbols are **scalars**
  - Nsymbol returns Bool
  - Fsymbol returns a single scalar value
  - Rsymbol returns a single scalar value
- **GenericFusion:** Nsymbol, F-symbols, and R-symbols include **multiplicity indices**
  - Nsymbol returns Integer (multiplicity > 1 possible)
  - Fsymbol returns rank-4 array with dimensions [da, db, dc] × [multiplicity_out]
    - Indices: `[F^{abc}_{d}]^{f}_{e}` where e, f label output/input multiplicities
  - Rsymbol returns rank-2 array (matrix) with dimensions [multiplicity_in] × [multiplicity_out]
    - Indices: `[R^{ab}_{c}]^{μ'}_{μ}` where μ, μ' label the multiplicities

**Registering with IrrepTable:**
- For group irreps, register via `Base.getindex(::IrrepTable, ::Type{G}) = IrrepType`
- Enables convenient `Irrep[G](...)` syntax
- Also support Unicode group symbols (ℤ₂, SU₂, etc.)

**Common Patterns:**
- Abelian groups: typically UniqueFusion, trivial F-symbols (all 1)
- Non-abelian groups: SimpleFusion or GenericFusion, non-trivial but structured F-symbols
- Anyons: SimpleFusion with irrational quantum dimensions

## Issue #57: Sector Implementation Conventions

**Goal:** Document specific conventions for each concrete sector type, especially labels and `fusiontensor` basis choices.

### Create overview page: `docs/src/sectors.md`

#### Opening Section: Group Representations vs Non-Group Sectors

Provide context before diving into specific implementations:

**Group Representations (Irreps):**
- Most sectors arise as irreducible representations (irreps) of groups
- Irreps are subtypes of `AbstractIrrep{G} <: Sector` where `G <: Group`
- Examples: ℤₙ, U(1), SU(2), dihedral groups, etc.
- Characterized by group-theoretic structure
- F-symbols and R-symbols derived from representation theory
- All have `BraidingStyle = Bosonic()`
- `AbelianIrrep`: subtypes for abelian groups (UniqueFusion, trivial F-symbols)

**Non-Group Sectors (Anyons):**
- Arise from more general (unitary, ribbon) fusion categories
- Examples: Fibonacci anyons, Ising anyons
- Not representations of any group
- Have intrinsic topological data (non-trivial F and R symbols)
- Characterized by `BraidingStyle = Anyonic()`
- Relevant for topological phases of matter and topological quantum computing
- Braiding described by braid group (not permutation group)

**Key Differences:**
- **Fusion:** Group irreps follow representation theory; anyons have more general fusion rules
- **Quantum dimensions:** Group irreps have positive integer dimensions; anyons can have irrational dimensions (e.g., golden ratio φ for Fibonacci)
- **Associators:** Group irreps typically have simple F-symbols; anyons encode intrinsic topological structure
- **Braiding:** Group reps use permutation group; anyons use braid group
- **Symmetric braiding:** For groups, R^{ab}_c R^{ba}_c = 1; not generally true for anyons

**Additional Categories:**
- **Planar (NoBraiding):** Fusion categories without braiding structure
- **Multi-fusion:** Categories where unit object is not simple (advanced topic)

This documentation focuses on the practical implementation details for each concrete sector type.

#### Navigation to Specific Pages

Link to the detailed documentation pages for each category of sectors.

### Create separate pages for each sector type:

#### Abelian Group Representations

- `docs/src/sectors/abelian/trivial.md` - `Trivial`
- `docs/src/sectors/abelian/zn.md` - `ZNIrrep{N}` and `LargeZNIrrep{N}` (including Z2, Z3, Z4)
- `docs/src/sectors/abelian/u1.md` - `U1Irrep`

#### Non-Abelian Group Representations

- `docs/src/sectors/nonabelian/su2.md` - `SU2Irrep`
- `docs/src/sectors/nonabelian/cu1.md` - `CU1Irrep`
- `docs/src/sectors/nonabelian/dn.md` - `DNIrrep{N}` (including D3, D4)
- `docs/src/sectors/nonabelian/a4.md` - `A4Irrep`

#### Anyonic Sectors

- `docs/src/sectors/anyons/planartrivial.md` - `PlanarTrivial`
- `docs/src/sectors/anyons/fibonacci.md` - `FibonacciAnyon`
- `docs/src/sectors/anyons/ising.md` - `IsingAnyon`

#### Fermionic Sectors

- `docs/src/sectors/fermions/parity.md` - `FermionParity`
- `docs/src/sectors/fermions/number.md` - `FermionNumber`
- `docs/src/sectors/fermions/spin.md` - `FermionSpin`

#### Composite Sectors

- `docs/src/sectors/composite/product.md` - `ProductSector`
- `docs/src/sectors/composite/timereversed.md` - `TimeReversed`

### Documentation Template

For each sector type, use a consistent documentation template:

```markdown
## SectorName

Brief one-line description.

### Type Definition
Show the struct definition with field explanations.

### Construction
How to create instances, including:
- Direct constructor
- `Irrep[G](...)` syntax (for group irreps)
- `Base.convert` convenience methods
- Unicode aliases

### Labels
How sectors are identified and what the parameters mean.

### Physical Interpretation
What this symmetry represents in physical systems (when applicable).

### Fusion Rules
Mathematical formula and implementation for ⊗ operation.

### Quantum Dimensions
Formula for `dim(a)` and any special values.

### Topological Data
- **F-symbols:** Key values, formulas, or references
- **R-symbols:** Key values, formulas, or references  
- **Braiding style:** Bosonic/Fermionic/Anyonic
- **Twists:** Values of `twist(a)`

### Basis Conventions
**Critical for fusiontensor:**
- State ordering within representations
- For group reps: generator matrices in the chosen basis
- Phase conventions (e.g., Condon-Shortley)
- References to external packages (e.g., WignerSymbols.jl)

### Iteration
- Order of `values(SectorType)` (for finite/small types)
- How `getindex` and `findindex` map indices to sectors

### Code Examples
```julia
# Concrete, runnable examples demonstrating:
# - Construction
# - Fusion
# - Quantum dimensions
# - Duality
# - Special properties
```

### Implementation Notes
- Storage details (e.g., UInt8 vs UInt for ZNIrrep)
- Performance considerations
- Related types
- External dependencies
```

Each page should be self-contained but cross-reference related sectors.

### Sectors to Document

#### Abelian Group Representations

**`Trivial` (`docs/src/sectors/abelian/trivial.md`):**
- `Trivial` - trivial representation of trivial group

**Cyclic Groups (`docs/src/sectors/abelian/zn.md`):**
- `ZNIrrep{N}` / `LargeZNIrrep{N}` - cyclic groups ℤₙ
  - Special cases: `Z2Irrep`, `Z3Irrep`, `Z4Irrep`
  - Label: integer n ∈ {0, 1, ..., N-1}
  - Fusion: (n₁ + n₂) mod N
  - Basis: all irreps are 1D, trivial basis

**Continuous Abelian (`docs/src/sectors/abelian/u1.md`):**
- `U1Irrep` - U(1) representations
  - Label: charge m ∈ ℤ or ½ℤ (HalfInt)
  - Fusion: m₁ + m₂
  - Iteration order: 0, 1/2, -1/2, 1, -1, 3/2, -3/2, ...
  - Basis: 1D representations

#### Non-Abelian Group Representations

**SU(2) (`docs/src/sectors/nonabelian/su2.md`):**
- `SU2Irrep` - SU(2) irreps
  - Label: spin j ∈ {0, 1/2, 1, 3/2, ...}
  - Fusion: Clebsch-Gordan series |j₁ - j₂| ⊕ ... ⊕ j₁ + j₂
  - Quantum dimension: 2j + 1
  - **Basis:** States ordered as m = j, j-1, ..., -j+1, -j
  - F-symbols via Racah W-coefficients
  - R-symbols: (-1)^(j₁+j₂-j₃)
  - **Generator basis (j=1/2):** Pauli matrices (specify normalization)
  - Clebsch-Gordan convention: follows WignerSymbols.jl (Condon-Shortley phases)

**U(1) ⊂ SU(2) (`docs/src/sectors/nonabelian/cu1.md`):**
- `CU1Irrep` - U(1) ⊂ SU(2) weight spaces
  - Label: (j, m) pairs
  - Fusion: based on SU(2) fusion with m₁ + m₂ constraint

**Dihedral Groups (`docs/src/sectors/nonabelian/dn.md`):**
- `DNIrrep{N}` - Dihedral group Dₙ
  - Special cases: `D3Irrep`, `D4Irrep`
  - Label conventions (differ for N even/odd)
  - Character table references

**Alternating Group (`docs/src/sectors/nonabelian/a4.md`):**
- `A4Irrep` - Alternating group A₄
  - Four irreps: three 1D, one 3D
  - Specific labels and fusion rules

#### Anyonic Sectors

**PlanarTrivial (`docs/src/sectors/anyons/planartrivial.md`):**
- Trivial anyon sector without braiding

**Fibonacci Anyons (`docs/src/sectors/anyons/fibonacci.md`):**
- `FibonacciAnyon`
  - Labels: `:I` (identity) and `:τ` (tau)
  - Fusion: τ ⊗ τ = I ⊕ τ
  - Quantum dimension: dim(τ) = φ (golden ratio)
  - **Key F-symbols:** F^{τ τ τ}_{τ I I}, F^{τ τ τ}_{τ τ τ}, etc.
  - **Key R-symbols:** R^I_{τ,τ}, R^τ_{τ,τ} with specific phases
  - Physical context: universal topological quantum computing

**Ising Anyons (`docs/src/sectors/anyons/ising.md`):**
- `IsingAnyon`
  - Labels: `:I`, `:σ`, `:ψ`
  - Fusion rules: ψ ⊗ ψ = I, σ ⊗ σ = I ⊕ ψ, σ ⊗ ψ = σ
  - Quantum dimensions: dim(σ) = √2, others = 1
  - F-symbols and R-symbols (key values)
  - Physical context: Majorana fermions, ν=5/2 quantum Hall

#### Fermionic Sectors

**Fermion Parity (`docs/src/sectors/fermions/parity.md`):**
- `FermionParity` (fℤ₂)
  - Labels: 0 (even/bosonic), 1 (odd/fermionic)
  - Fusion: like ℤ₂ but with fermionic braiding
  - R-symbol: -1 for 1 ⊗ 1

**Fermion Number (`docs/src/sectors/fermions/number.md`):**
- `FermionNumber` (fU₁)
  - Labels: integer fermion number
  - Fusion: additive
  - R-symbol: (-1)^(m₁·m₂)

**Fermion Spin (`docs/src/sectors/fermions/spin.md`):**
- `FermionSpin` (fSU₂)
  - Labels: spin j with fermion parity from 2j
  - Fusion: like SU₂ with fermionic statistics

#### Composite Sectors

**Product Sectors (`docs/src/sectors/composite/product.md`):**
- `ProductSector` (A × B)
  - Combining multiple symmetries
  - Notation: ⊠ operator
  - Fusion: component-wise

**Time-Reversed Sectors (`docs/src/sectors/composite/timereversed.md`):**
- `TimeReversed{I}`
  - Anti-unitary symmetries
  - Conjugate braiding

### Summary Table

Include a comparison table (on the overview `sectors.md` page) with columns:
- Sector Type
- Symmetry Group/Category
- Fusion Style
- Braiding Style
- Infinite Sectors?
- Use Cases

## Additional Documentation Improvements

### Enhance `docs/src/index.md`

Add to the existing content:

1. **Quick Start section** (2-3 examples)
   - Creating sectors
   - Fusion product
   - Checking quantum dimensions

2. **Documentation Guide section**
   - Links to Interface documentation
   - Links to Sector conventions
   - Link to Library reference

3. Keep introduction concise but add context about:
   - Physical applications (quantum mechanics, condensed matter, quantum computing)
   - Mathematical background (fusion categories)

### Update `docs/make.jl`

Add new pages to navigation in logical order:
```julia
pages = [
    "Home" => "index.md",
    "Sector Interface" => "interface.md",
    "Sector Types" => [
        "Overview" => "sectors.md",
        "Abelian Groups" => [
            "Trivial" => "sectors/abelian/trivial.md",
            "ℤₙ (Cyclic)" => "sectors/abelian/zn.md",
            "U₁" => "sectors/abelian/u1.md",
        ],
        "Non-Abelian Groups" => [
            "SU₂" => "sectors/nonabelian/su2.md",
            "CU₁" => "sectors/nonabelian/cu1.md",
            "Dₙ (Dihedral)" => "sectors/nonabelian/dn.md",
            "A₄ (Alternating)" => "sectors/nonabelian/a4.md",
        ],
        "Anyonic Sectors" => [
            "Fibonacci" => "sectors/anyons/fibonacci.md",
            "Ising" => "sectors/anyons/ising.md",
        ],
        "Fermionic Sectors" => [
            "Fermion Parity" => "sectors/fermions/parity.md",
            "Fermion Number" => "sectors/fermions/number.md",
            "Fermion Spin" => "sectors/fermions/spin.md",
        ],
        "Composite Sectors" => [
            "Product" => "sectors/composite/product.md",
            "Time-Reversed" => "sectors/composite/timereversed.md",
        ],
    ],
    "Library" => "lib.md",
]
```

## Key Principles

### Specificity for Basis Conventions
Be very explicit about:
- State ordering in representations (e.g., m = j, j-1, ..., -j for SU₂)
- Phase conventions (e.g., Condon-Shortley for Clebsch-Gordan)
- Generator matrices in standard basis (for testability)
- References to external conventions (e.g., WignerSymbols.jl)

### Testability
- Provide generator matrices where applicable so users can verify
- Example: For `SU2Irrep` j=1/2, give explicit Pauli matrices
- Users can check their code against these

### Completeness
- Every exported sector type must be documented
- Include both common and specialized sectors

### Balance
- Interface doc: conceptual and tutorial-like
- Sectors doc: reference-style with specific conventions
- Library doc: API reference (already exists)

## Implementation Order

1. Create `interface.md` with complete Sector interface documentation
2. Create `sectors.md` overview page with group vs non-group explanation and summary table
3. Create individual sector documentation pages in this order:
   - Abelian groups: `trivial.md`, `zn.md`, `u1.md`
   - Non-abelian groups: `su2.md`, `cu1.md`, `dn.md`, `a4.md`
   - Anyons: `fibonacci.md`, `ising.md`
   - Fermions: `parity.md`, `number.md`, `spin.md`
   - Composite: `product.md`, `timereversed.md`
4. Enhance `index.md` with Quick Start and navigation
5. Update `make.jl` to include all new pages
6. Review and test all code examples

## Open Questions

Before implementation, clarify:

1. **Page naming:** Should conventions be in `sectors.md` or `conventions.md`?
2. **Generator matrices detail level:** How many examples? (e.g., just j=1/2 for SU₂, or also j=1?)
3. **Priority sectors:** Any specific sector types needing more detail?
4. **Mathematical depth:** Include background on pentagon/hexagon equations, or keep purely practical?
5. **Code examples:** Should all examples be runnable, or is pseudocode acceptable?
6. **External references:** Level of detail for linking to external conventions (WignerSymbols.jl, literature)?
