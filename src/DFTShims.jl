__precompile__(true)
module DFTShims
using Unitful

include("UnitfulHartree.jl")
using .UnitfulHartree
const UH = UnitfulHartree

include("Dispatch.jl")
using .Dispatch

include("Traits.jl")
using .Traits

include("ArrayInitialization.jl")
using .ArrayInitialization

export UnitfulHartree, Unitful, @u_str
export DFTArray, DFTAxisArray, dft_similar, Dispatch
export has_axis, is_spin_polarized, SpinCategory, SpinDegenerate, components
export ColinearSpin, ColinearSpinFirst, ColinearSpinPreferLast, FunctionalCategory, LDA, GGA

end # module
