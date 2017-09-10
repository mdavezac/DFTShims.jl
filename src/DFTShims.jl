__precompile__(true)
module DFTShims
using Unitful
using AxisArrays

if Pkg.installed("AxisArrays") < v"0.1.5"
    include("AxisArrays.v0.1.x.jl")
end


include("UnitfulHartree.jl")
using .UnitfulHartree
const UH = UnitfulHartree
export UnitfulHartree, Unitful, @u_str

include("Dispatch.jl")
using .Dispatch
using .Dispatch: DFTArray, DFTAxisArray
export DFTArray, DFTAxisArray, Dispatch

include("Traits.jl")
using .Traits: ColinearSpin, ColinearSpinFirst, ColinearSpinLast, SpinDegenerate,
               SpinCategory, FunctionalCategory, GGA, LDA, components, is_spin_polarized,
               has_axis, concretize_type, SpinAware
export ColinearSpin, ColinearSpinFirst, ColinearSpinLast, SpinDegenerate, SpinCategory,
       FunctionalCategory, GGA, LDA, components, is_spin_polarized, has_axis,
       concretize_type, SpinAware

include("ConstantArrays.jl")
using .ConstantArrays

include("ArrayInitialization.jl")
using .ArrayInitialization: wrap
export wrap

end # module
