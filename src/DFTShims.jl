__precompile__(true)
module DFTShims
using Unitful

include("UnitfulHartree.jl")
using .UnitfulHartree
const UH = UnitfulHartree
export UnitfulHartree, Unitful, @u_str

include("Dispatch.jl")
using .Dispatch
using .Dispatch: DFTArray, DFTAxisArray
export DFTArray, DFTAxisArray, Dispatch

include("Traits.jl")
using .Traits: ColinearSpin, ColinearSpinFirst, ColinearSpinLast, ColinearSpinPreferLast,
               SpinDegenerate, SpinCategory, FunctionalCategory, GGA, LDA, components,
               is_spin_polarized, has_axis
export ColinearSpin, ColinearSpinFirst, ColinearSpinLast, ColinearSpinPreferLast,
       SpinDegenerate, SpinCategory, FunctionalCategory, GGA, LDA, components,
       is_spin_polarized, has_axis


include("ArrayInitialization.jl")
using .ArrayInitialization

end # module
