module TensorKitSectors

# exports
# -------
export Sector, Group, AbstractIrrep, AbstractGroupElement
export Irrep, Element

export Nsymbol, Fsymbol, Rsymbol, Asymbol, Bsymbol
export sectorscalartype
export dim, sqrtdim, invsqrtdim, frobeniusschur, twist, fusiontensor, dual
export otimes, deligneproduct, times
export FusionStyle, UniqueFusion, MultipleFusion, SimpleFusion, GenericFusion,
    MultiplicityFreeFusion
export BraidingStyle, NoBraiding, HasBraiding, SymmetricBraiding, Bosonic, Fermionic, Anyonic
export SectorSet, SectorValues, findindex
export rightone, leftone

export triangle_equation, pentagon_equation, hexagon_equation

export Trivial, Z2Irrep, Z3Irrep, Z4Irrep, ZNIrrep, U1Irrep, SU2Irrep, CU1Irrep
export ZNElement, Z2Element, Z3Element, Z4Element
export ProductSector, TimeReversed
export FermionParity, FermionNumber, FermionSpin
export PlanarTrivial, FibonacciAnyon, IsingAnyon
export IsingBimodule

# unicode exports
# ---------------
export ⊠, ⊗, ×
export ℤ, ℤ₂, ℤ₃, ℤ₄, U₁, SU, SU₂, CU₁
export fℤ₂, fU₁, fSU₂

# imports
# -------
using Base: SizeUnknown, HasLength, IsInfinite
using Base: HasEltype, EltypeUnknown
using Base.Iterators: product, filter
using Base: tuple_type_head, tuple_type_tail

using LinearAlgebra: tr
using TensorOperations
using HalfIntegers
using WignerSymbols

# includes
# --------
include("auxiliary.jl")
include("sectors.jl")
include("trivial.jl")
include("groups.jl")
include("irreps.jl")    # irreps of symmetry groups, with bosonic braiding
include("groupelements.jl") # group elements with cocycles, no braiding
include("product.jl")   # direct product of different sectors
include("fermions.jl")  # irreps with defined fermionparity and fermionic braiding
include("anyons.jl")    # non-group sectors
include("multifusion.jl") # multifusion example, namely Rep Z2 ⊕ Rep Z2 ≅ Ising

# precompile
# ----------
include("precompile.jl")

function __precompile__()
    for I in (
            Trivial, Z2Irrep, Z3Irrep, Z4Irrep, ZNIrrep, U1Irrep, SU2Irrep, CU1Irrep,
            FermionParity, FermionNumber, FermionSpin, PlanarTrivial, FibonacciAnyon,
            IsingAnyon, TimeReversed{IsingAnyon}, TimeReversed{FibonacciAnyon},
        )
        precompile_sector(I)
    end
    return
end

end # module TensorKitSectors
