# TensorKitSectors

A Julia package for working with objects in fusion categories.

| **Build Status** | **PkgEval** | **Coverage** | **Quality assurance** |
|:----------------:|:------------:|:------------:|:---------------------:|
| [![CI][ci-img]][ci-url] | [![PkgEval][pkgeval-img]][pkgeval-url] | [![Codecov][codecov-img]][codecov-url] | [![Aqua QA][aqua-img]][aqua-url] |

[ci-img]: https://github.com/QuantumKitHub/TensorKitSectors.jl/actions/workflows/CI.yml/badge.svg
[ci-url]: https://github.com/QuantumKitHub/TensorKitSectors.jl/actions/workflows/CI.yml

[pkgeval-img]: https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/T/TensorKitSectors.svg
[pkgeval-url]: https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/T/TensorKitSectors.html

[codecov-img]: https://codecov.io/gh/QuantumKitHub/TensorKitSectors.jl/branch/main/graph/badge.svg
[codecov-url]: https://codecov.io/gh/QuantumKitHub/TensorKitSectors.jl

[aqua-img]: https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg
[aqua-url]: https://github.com/JuliaTesting/Aqua.jl

This package provides functionality for defining objects in fusion categories, along with their topological data.
This includes the fusion rules, the associators, and the braiding.
In particular, this is the data that is needed to define (symmetric) tensors, which are defined over vector spaces graded by these objects.
For the full functionality, we refer to [TensorKit.jl](https://github.com/QuantumKitHub/TensorKit.jl) and [its documentation](https://quantumkithub.github.io/TensorKit.jl/stable/).

Install via the package manager.
