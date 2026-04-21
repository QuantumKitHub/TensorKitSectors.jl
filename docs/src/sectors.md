# Sector Types

This page provides an overview of the concrete sector types implemented in TensorKitSectors.jl.

## Summary Table

| Sector Type              | Group/Category           | Fusion Style | Braiding  | Infinite? | Common Use Cases                      |
|--------------------------|--------------------------|--------------|-----------|-----------|---------------------------------------|
| [`Trivial`](@ref)        | Trivial group            | Unique       | Bosonic   | No        | No symmetry                           |
| [`ZNIrrep`](@ref)        | ℤₙ (cyclic)              | Unique       | Bosonic   | No        | Clock models, discrete symmetries     |
| [`U1Irrep`](@ref)        | U(1)                     | Unique       | Bosonic   | Yes       | Particle number, charge conservation  |
| [`SU2Irrep`](@ref)       | SU(2)                    | Simple       | Bosonic   | Yes       | Spin systems, angular momentum        |
| [`CU1Irrep`](@ref)       | U(1) ⋊ ℤ₂                | Simple       | Bosonic   | Yes       | Particle-hole symmetry, O(2)          |
| [`DNIrrep`](@ref)        | Dₙ (dihedral)            | Simple       | Bosonic   | No        | Molecular/crystal symmetries          |
| [`A4Irrep`](@ref)        | A₄ (alternating)         | Generic      | Bosonic   | No        | Tetrahedral symmetry                  |
| [`FibonacciAnyon`](@ref) | Fibonacci category       | Simple       | Anyonic   | No        | Topological quantum computing         |
| [`IsingAnyon`](@ref)     | Ising category           | Simple       | Anyonic   | No        | Majorana fermions, ν=5/2 QHE          |
| [`FermionParity`](@ref)  | fℤ₂                      | Unique       | Fermionic | No        | Fermion parity conservation           |
| [`FermionNumber`](@ref)  | fU₁                      | Unique       | Fermionic | Yes       | Fermion number conservation           |
| [`FermionSpin`](@ref)    | fSU₂                     | Simple       | Fermionic | Yes       | Fermions with spin symmetry           |
| [`ProductSector`](@ref)  | Product categories       | Varies       | Varies    | Varies    | Multiple simultaneous symmetries      |
| [`TimeReversed`](@ref)   | Inverted braiding        | Varies       | Varies    | Varies    | ??                                    |
