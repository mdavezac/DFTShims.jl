using DFTShims
using Base.Test
using Unitful

const DD = DFTShims.Dispatch.Dimensions
const DH = DFTShims.Dispatch.Hartree
const UH = DFTShims.UnitfulHartree

@testset "DFTshims" begin
    @testset "Dispatch" begin
        @test 1u"ρ" isa DD.Scalars.ρ
        @test 1u"ρ" isa DH.Scalars.ρ

        @test DH.Scalars.ρ <: DD.Scalars.ρ
        @test DH.Scalars.ρ{Float32} <: DD.Scalars.ρ
        @test DH.Scalars.ρ{Float32} <: DD.Scalars.ρ{Float32}
        @test DH.Scalars.ρ{Float32} === DD.Scalars.ρ{Float32, typeof(UH.ρ)}

        @test [1u"ρ", 2.0u"ρ"] isa DFTShims.DFTArray
        @test [1u"ρ", 2.0u"ρ"] isa DD.DensityArray
        @test [1u"ρ", 2.0u"ρ"] isa DD.Arrays.ρ
        @test [1u"ρ", 2.0u"ρ"] isa DD.Arrays.ρ{Float64}
        @test [1u"ρ", 2.0u"ρ"] isa DH.DensityDenseArray{Float64, 1}
        @test [1u"ρ", 2.0u"ρ"] isa DH.DenseArrays.ρ{Float64, 1}

        @test DH.DensityDenseArray{Int8, 2} <: DD.DensityDenseArray{Int8, 2, typeof(u"ρ")}
    end
end
