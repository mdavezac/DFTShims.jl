using DFTShims
using Base.Test
using Unitful
using AxisArrays

const DD = DFTShims.Dispatch.Dimensions
const DH = DFTShims.Dispatch.Hartree
const UH = DFTShims.UnitfulHartree

@testset "Dispatch" begin
    include("dispatch.jl")
end
@testset "Traits" begin
    include("Traits.jl")
end
@testset "Initialization" begin
    include("Initialization.jl")
end
@testset "Constant arrays" begin
    include("ConstantArrays.jl")
end
nothing
