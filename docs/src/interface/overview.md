# Sector Interface

## Introduction

A `Sector` is a label for the different *symmetry charges* or *quantum numbers* that organize graded vector spaces.
In particular, these are the labels use to decompose vector spaces as:

```math
V = \bigoplus_{a \in \text{Sectors}} \mathbb{C}^{n_a} \otimes V_a
```

where ``n_a`` is the multiplicity of sector ``a`` and ``V_a`` is the associated vector space.

Sectors encode the structural rules TensorKit needs:
- **Objects**: the labels themselves (irreps, anyons, …)
- **Fusion rules**: which labels appear in ``a \otimes b \rightarrow \bigoplus_c N^{ab}_c c``
- **Associativity**: how different parenthesizations are related
- **Braiding**: how labels behave under exchange (if supported)

More rigorously, a `Sector` represents the isomorphism classes of simple objects in unitary and pivotal fusion categories.
This package defines an interface for accessing the topological data that is associated to these categories.
This section explains the required and optional methods needed to create a new sector type.
Once this interface is fulfilled, TensorKit.jl will create symmetric tensors that correspond to the implemented sector.

## Organization

The interface documentation is divided into several pages:

- **[Required Methods](required.md)**: The minimum interface needed for a valid sector type
- **[Optional Methods](optional.md)**: Additional methods with default implementations that can be specialized
- **[Traits and Styles](traits.md)**: Compile-time properties that control behavior and optimizations
- **[Implementation Guidelines](guidelines.md)**: Practical advice and helper types for implementing sectors
