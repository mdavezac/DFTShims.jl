module Dispatch

export DFTArray, DFTAxisArray
using AxisArrays
using Unitful

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
    :Density                         => :ρ,
    :ContractedDensityGradient       => :σ,
    :FirstDensityDerivative          => :∂ϵ_∂ρ,
    :FirstGradientDerivative         => :∂ϵ_∂σ,
    :SecondDensityDerivative         => :∂²ϵ_∂ρ²,
    :SecondGradientDerivative        => :∂²ϵ_∂σ²,
    :SecondDensityGradientDerivative => :∂²ϵ_∂ρ∂σ,
    :ThirdDensityDerivative          => :∂³ϵ_∂ρ³,
    :ThirdGradientDerivative         => :∂³ϵ_∂σ³,
    :ThirdDensity2GradientDerivative => :∂³ϵ_∂ρ²∂σ,
    :ThirdDensityGradient2Derivative => :∂³ϵ_∂ρ∂σ²,
    :Potential                       => :Eₕ
)
include("DispatchDimensions.jl")
include("DispatchHartree.jl")
end
