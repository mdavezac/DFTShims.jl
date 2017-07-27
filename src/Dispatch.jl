module Dispatch

export DFTArray, DFTAxisArray
using AxisArrays

macro lintpragma(s) end
@lintpragma("Ignore use of undeclared variable T")
@lintpragma("Ignore use of undeclared variable N")


""" An axis array with dense memory allocation """
const DFTAxisArray = AxisArray{T, N, <: DenseArray{T, N}} where N where T
""" Set of array types used in LibXC and AtomicDFT """
const DFTArray =
    Union{DenseArray{T, N}, AxisArray{T, N, <: DenseArray{T, N}}} where N where T 

""" Dimensions and units of interest for DFT """
const Dimensioned = (
    :Density => :ρ, :Potential => :Eₕ
)
include("DispatchDimensions.jl")
include("DispatchHartree.jl")

end
