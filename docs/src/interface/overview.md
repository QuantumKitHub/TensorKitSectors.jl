# Sector Interface

## Introduction

A `Sector` represents the isomorphism classes of simple objects in unitary and pivotal fusion categories.
In the context of tensor networks and representation theory, sectors label the different *symmetry charges* or *quantum numbers* that organize graded vector spaces.

Mathematically, sectors encode the topological data of fusion categories:
- **Simple objects**: Irreducible representations or anyonic excitations
- **Fusion rules**: How sectors combine (``a ⊗ b → \bigoplus_c N^{ab}_c c``)
- **Associativity**: F-symbols (6j symbols) encoding recoupling transformations
- **Braiding**: R-symbols encoding exchange statistics (for braided categories)

Sectors can be used as labels for graded vector spaces, which represent vector spaces decomposed into different sectors:
```math
V = \bigoplus_{a \in \text{Sectors}} \mathbb{C}^{n_a} \otimes V_a
```
where ``n_a`` is the multiplicity of sector ``a`` and ``V_a`` is the associated vector space.

This section explains the required and optional methods needed to create a new sector type.
Once this interface is fulfilled, TensorKit.jl should be able to create symmetric tensors that correspond to the implemented sector.

## Organization

The interface documentation is divided into several pages:

- **[Required Methods](required.md)**: The minimum interface needed for a valid sector type
- **[Optional Methods](optional.md)**: Additional methods with default implementations that can be specialized
- **[Traits and Styles](traits.md)**: Compile-time properties that control behavior and optimizations
- **[Implementation Guidelines](guidelines.md)**: Practical advice and helper types for implementing sectors
