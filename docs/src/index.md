# TensorKitSectors.jl

A Julia package for working with objects in fusion categories.

This package provides functionality for defining objects in fusion categories, along with their topological data.
This includes the fusion rules, the associators, and the braiding.
In particular, this is the data that is needed to define (symmetric) tensors, which are defined over vector spaces graded by these objects.
For the full functionality, we refer to [TensorKit.jl](https://github.com/QuantumKitHub/TensorKit.jl) and [its documentation](https://quantumkithub.github.io/TensorKit.jl/stable/).

## Installation

Install via the package manager, by entering the Pkg REPL mode with `]` and running:

```julia-repl
pkg> add TensorKitSectors
```

## Where to go next

- **[Sector Interface](@ref)** — the methods a sector type must (and may) implement, and the traits that control its behavior. Start here to implement a new sector.
- **[Sector Types](@ref)** — the concrete sectors shipped with this package (groups, anyons, fermions, and composites), with a summary table and a page per type.
