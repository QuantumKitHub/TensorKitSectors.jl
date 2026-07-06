# Sector Types

This page provides an overview of the concrete sector types implemented in TensorKitSectors.jl.

## Summary Table

| Sector Type              | Group/Category                | Fusion Style | Braiding  | Infinite? | Common Use Cases                                          |
|--------------------------|--------------------------------|--------------|-----------|-----------|-------------------------------------------------------------|
| [`Trivial`](@ref)        | Trivial group                  | Unique       | Bosonic   | No        | No symmetry                                                |
| [`ZNIrrep`](@ref)        | Rep[ℤₙ]               | Unique       | Bosonic   | No        | Clock models, discrete symmetries                          |
| [`ZNElement`](@ref)      | Vec[ℤₙ]  | Unique       | Varies    | No        | Group-graded categories, discrete torsion  |
| [`U1Irrep`](@ref)        | Rep[U₁]                            | Unique       | Bosonic   | Yes       | Particle number, charge conservation                       |
| [`SU2Irrep`](@ref)       | Rep[SU₂]                           | Simple       | Bosonic   | Yes       | Spin systems, angular momentum                             |
| [`CU1Irrep`](@ref)       | Rep[U₁ ⋊ ℤ₂]                       | Simple       | Bosonic   | Yes       | Particle-hole symmetry, O(2)                               |
| [`DNIrrep`](@ref)        | Rep[Dₙ]                   | Simple       | Bosonic   | No        | Molecular/crystal symmetries                               |
| [`A4Irrep`](@ref)        | Rep[A₄]                | Generic      | Bosonic   | No        | Tetrahedral symmetry                                       |
| [`HeisenbergIrrep`](@ref)| Rep[H_N]  | Generic      | Bosonic   | No        | Weyl-Heisenberg symmetry, projective representations |
| [`FibonacciAnyon`](@ref) | Fibonacci category              | Simple       | Anyonic   | No        | Topological quantum computing                              |
| [`IsingAnyon`](@ref)     | Ising category                  | Simple       | Anyonic   | No        | Majorana fermions, ν=5/2 QHE                               |
| [`FermionParity`](@ref)  | fℤ₂                             | Unique       | Fermionic | No        | Fermion parity conservation                                |
| [`FermionNumber`](@ref)  | fU₁                             | Unique       | Fermionic | Yes       | Fermion number conservation                                |
| [`FermionSpin`](@ref)    | fSU₂                            | Simple       | Fermionic | Yes       | Fermions with spin symmetry                                |
| [`ProductSector`](@ref)  | Product categories              | Varies       | Varies    | Varies    | Multiple simultaneous symmetries                           |
| [`TimeReversed`](@ref)   | Inverted braiding               | Varies       | Varies    | Varies    | Time-reversal or other orientation-reversing symmetries    |

## Other packages

TensorKitSectors.jl provides the architecture for implementing new sector types, but does not implement all possible sectors. Other packages which implement additional sector types and their topological data include:
- [`SUNRepresentations.jl`](https://github.com/QuantumKitHub/SUNRepresentations.jl): Implements irreducible representations of SU(N) for arbitrary N
- [`CategoryData.jl`](https://github.com/QuantumKitHub/CategoryData.jl): Provides a variety of fusion categories up to rank 7, based on the [`AnyonWiki`](https://anyonwiki.github.io/)
- [`QWignerSymbols.jl`](https://github.com/QuantumKitHub/QWignerSymbols.jl): Provides the irreps of q-deformed SU(2)
