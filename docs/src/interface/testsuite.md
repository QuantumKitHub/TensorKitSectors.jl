```@meta
CollapsedDocStrings = true
```

# Testing a Sector Implementation

Implementing a new `Sector` subtype means providing a family of mutually consistent methods: fusion rules, F- and R-symbols, dimensions, dualities, an ordering, ...
These have to satisfy several non-trivial categorical identities (the pentagon and hexagon equations, unitarity of the F- and R-move, ...), but must also be compatible with TensorKit.jl as a data structure.
TensorKitSectors.jl ships a reusable test suite, `SectorTestSuite`, that checks exactly these properties for any `Sector` subtype.
It is defined within the package's own tests, but can be accessed by downstream packages or users testing their own sector implementations locally.

Here, we explain how to use the test suite, what it checks, and how to add new tests.
Currently, the suite is designed to run all tests for a single sector type at a time, as every single test within the test suite necessarily must pass for compatibility with TensorKit.jl.

## Running the test suite

Since the test suite is not part of the installed package, it is loaded as follows:

```julia
import TensorKitSectors
testsuite_path = joinpath(
    dirname(dirname(pathof(TensorKitSectors))), # TensorKitSectors root
    "test", "testsuite.jl"
)
include(testsuite_path)

SectorTestSuite.test_sector(MySectorType)
```

`test_sector` runs one `@testsuite` per property (see [What gets tested](@ref) below) and reports which ones fail for `MySectorType`.
The suite relies on [`TestExtras.jl`](https://github.com/Jutho/TestExtras.jl) for `@testinferred`, which checks both the return value and the type stability of an expression, so downstream packages should add `TestExtras` as a test dependency.


## What gets tested

Broadly, the sector tests fall into three categories:

 Category  | Checks 
---|---
Interface & type stability | 1) Basic properties: the required methods exist, are type-stable, and return the documented types <br> 2) Show and parse: `Sector`s are printed in a parseable manner <br> 3) Value iterator: `Sector`s are ordered within `values(I)` in a consistent and expected manner (through `findindex`)
Category-theoretic consistency | The algebraic identities a unitary (braided) fusion category must satisfy: pentagon and hexagon equations, unitarity of the F- and R-move, triangle equation, ribbon condition, self-duality of the braiding, symmetric braiding condition, Artin braid equality
Cross-consistency of derived quantities | Whenever a quantity can be computed both directly and through a generic fallback derived from other data (e.g. [`dim`](@ref) versus the internal `dim_from_Fsymbol`, or [`Bsymbol`](@ref) versus `Bsymbol_from_fusiontensor`), the two must agree

Tests that only apply to a subset of sectors are skipped automatically based on traits: braiding-related testsets pass untested unless `BraidingStyle(I) isa HasBraiding`, the fusion-tensor comparisons are skipped when `fusiontensor` has no method for `I`, and so on.
A test therefore never fails simply because a sector chooses not to implement an optional method.

## Test utilities

Most tests need a handful of representative sectors rather than the full (possibly infinite) set of values of `I`.
`SectorTestSuite` exports a few helpers to construct those:

- `smallset(I, size=5, maxdim=10)`: a small, shuffled sample of sectors of type `I` with dimension below `maxdim`, biased to include a non-abelian sector when `FusionStyle(I) isa MultipleFusion` so tests don't silently degenerate to the abelian case.
- `randsector(I)`: draws a single, non-unit sector at random from `smallset(I)`.
- `random_fusion(I, N)`: draws `N` random sectors such that every consecutive pair can be fused, i.e. builds a valid fusion series `a ⊗ b ⊗ ...`.

Besides this, other utility functions are:
- `can_fuse(a, b)`: whether `a ⊗ b` is non-empty.
  Allows to test only non-trivial topological data, e.g. the F-symbols of fusion vertices that actually exist.
- `hasfusiontensor(I)`: whether [`fusiontensor`](@ref) has an implementation for `I`.

Additionally, the following tests are provided as standalone functions that can be used outside of the test suite:
- `F_unitarity_test(a, b, c; kwargs...)` / `R_unitarity_test(a, b; kwargs...)`: check unitarity of the F- and R-move, forwarding `kwargs` to `isapprox`.

These are all exported by `SectorTestSuite` and are also useful when writing additional, sector-specific tests beyond what the shared suite covers.

## Adding a new test to the test suite

New checks are registered with the `@testsuite` macro, which takes a name and a function of the sector type `I`:

```julia
@testsuite "My new property" I -> begin
    # test body
end
```
