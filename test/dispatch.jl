using DFTShims
using Base.Test
using Unitful
using AxisArrays

const DD = DFTShims.Dispatch.Dimensions
const DH = DFTShims.Dispatch.Hartree
const UH = DFTShims.UnitfulHartree

@testset "Scalars" begin
    @test 1u"ρ" isa DD.Scalars.ρ
    @test 1u"nm^-3" isa DD.Scalars.ρ
    @test 1u"ρ" isa DH.Scalars.ρ
    @test !(1u"nm^-3" isa DH.Scalars.ρ)

    @test DH.Scalars.ρ <: DD.Scalars.ρ
    @test DH.Scalars.ρ{Float32} <: DD.Scalars.ρ
    @test DH.Scalars.ρ{Float32} <: DD.Scalars.ρ{Float32}
    @test DH.Scalars.ρ{Float32} === DD.Scalars.ρ{Float32, typeof(UH.ρ)}
    @test DH.Scalars.ρ{Float32} <: DD.Scalars.All{Float32}
    @test !(DH.Scalars.ρ{Int32} <: DD.Scalars.All{Float32})
    @test !(DH.Scalars.ρ <: DD.Scalars.All{Float32})

    f = (::DH.Scalars.ρ{<: Integer}) -> :hartree_with_integer
    (::typeof(f))(::DH.Scalars.ρ) = :hartree_with_any_number
    (::typeof(f))(::DD.Scalars.ρ{<: Integer}) = :density_with_integer
    (::typeof(f))(::DD.Scalars.ρ) = :density_with_any_number
    @test f(1u"ρ") == :hartree_with_integer
    @test f(1.0u"ρ") == :hartree_with_any_number
    @test f(1u"nm^-3") == :density_with_integer
    @test f(1.0u"nm^-3") == :density_with_any_number
end


@testset "Arrays" begin
    @test [1u"ρ", 2.0u"ρ"] isa DFTShims.DFTArray
    @test [1u"ρ", 2.0u"ρ"] isa DD.DensityArray
    @test [1u"ρ", 2.0u"ρ"] isa DD.Arrays.ρ
    @test [1u"ρ", 2.0u"ρ"] isa DD.Arrays.ρ{Float64}
    @test [1u"ρ", 2.0u"ρ"] isa DH.DensityDenseArray{Float64, 1}
    @test [1u"ρ", 2.0u"ρ"] isa DH.DenseArrays.ρ{Float64, 1}
    @test [1u"nm^-3", 2.0u"nm^-3"] isa DD.DenseArrays.ρ{Float64, 1}
    @test !([1u"nm^-3", 2.0u"nm^-3"] isa DH.DenseArrays.ρ{Float64, 1})

    @test DH.DensityDenseArray{Int8, 2} === DD.DensityDenseArray{Int8, 2, typeof(u"ρ")}
    @test DH.Arrays.ρ{Float32} <: DD.Arrays.All{Float32}
    @test !(DH.Arrays.ρ{Int32} <: DD.Arrays.All{Float32})
    @test !(DH.Arrays.ρ <: DD.Arrays.All{Float32})
    @test DH.Arrays.ρ{Float32, 2} <: DD.Arrays.All{Float32, 2}
    @test [1u"nm^-3", 2.0u"nm^-3"] isa DD.Arrays.ρ{Float64, 1}
    @test !([1u"nm^-3", 2.0u"nm^-3"] isa DD.Arrays.ρ{Float64, 2})

    f = (::DH.Arrays.ρ{<: Integer}) -> :hartree_with_integer
    (::typeof(f))(::DH.Arrays.ρ) = :hartree_with_any_number
    (::typeof(f))(::DD.Arrays.ρ{<: Integer}) = :density_with_integer
    (::typeof(f))(::DD.Arrays.ρ) = :density_with_any_number
    (::typeof(f))(::DD.AxisArrays.ρ) = :axis_density_with_any_number
    @test f([1u"ρ"]) == :hartree_with_integer
    @test f([1.0u"ρ"]) == :hartree_with_any_number
    @test f([1u"nm^-3"]) == :density_with_integer
    @test f([1.0u"nm^-3"]) == :density_with_any_number
    @test f(AxisArray([1.0u"nm^-3"])) == :axis_density_with_any_number
end
