using DFTShims
using AxisArrays
using Unitful

spin = AxisArray(zeros((10, 2)), Axis{:radius}(0u"m":1u"m":9u"m"), Axis{:spin}((:+, :-)))
nospin = AxisArray(zeros((10, 2)), Axis{:radius}(0u"m":1u"m":9u"m"), Axis{:n}((1, 2)))

@testset "has_axis" begin
    @test has_axis(spin, :spin)
    @test has_axis(spin, Axis{:spin})
    @test has_axis(typeof(spin), :spin)
    @test has_axis(typeof(spin), Axis{:spin})
    @test !has_axis(nospin, :spin)
    @test !has_axis(nospin, Axis{:spin})
    @test !has_axis(typeof(nospin), :spin)
    @test !has_axis(typeof(nospin), Axis{:spin})
end

@testset "spin polarization" begin
    @test is_spin_polarized(spin)
    @test !is_spin_polarized(nospin)

    @test is_spin_polarized(typeof(spin))
    @test !is_spin_polarized(typeof(nospin))

    @test @inferred(PolarizationCategory(spin)) === Polarized
    @test @inferred(PolarizationCategory(nospin)) === Unpolarized
end

@testset "components" begin
    @test @inferred(components(UnitfulHartree.ρ, Unpolarized)) == (:ρ,)
    @test all(components(UnitfulHartree.ρ, Polarized) .== (:α, :β))
    @test length(components(Dispatch.Dimensions.Scalars.∂ϵ_∂σ, Unpolarized)) == 1
    @test length(@inferred(components(UnitfulHartree.∂ϵ_∂σ, Polarized))) == 3
    @test length(@inferred(components(UnitfulHartree.∂²ϵ_∂ρ², Polarized))) == 3
    @test length(@inferred(components(UnitfulHartree.∂²ϵ_∂ρ∂σ, Unpolarized))) == 1
    @test length(@inferred(components(UnitfulHartree.∂²ϵ_∂ρ∂σ, Polarized))) == 6
    @test length(@inferred(components(UnitfulHartree.∂²ϵ_∂σ², Polarized))) == 6
    @test length(@inferred(components(UnitfulHartree.∂³ϵ_∂ρ³, Polarized))) == 4
    @test length(@inferred(components(UnitfulHartree.∂³ϵ_∂ρ²∂σ, Polarized))) == 9
    @test length(@inferred(components(UnitfulHartree.∂³ϵ_∂ρ∂σ², Polarized))) == 12
    @test length(@inferred(components(UnitfulHartree.∂³ϵ_∂σ³, Polarized))) == 10
end

@testset "functional category" begin
    @test FunctionalCategory(UnitfulHartree.ρ) === LDA
    @test FunctionalCategory(UnitfulHartree.∂ϵ_∂σ) === GGA
    @test FunctionalCategory(UnitfulHartree.∂²ϵ_∂ρ²) === LDA
    @test FunctionalCategory(UnitfulHartree.∂²ϵ_∂ρ∂σ) === GGA
    @test FunctionalCategory(Dispatch.Dimensions.Scalars.∂²ϵ_∂σ²) === GGA
    @test FunctionalCategory(UnitfulHartree.∂³ϵ_∂ρ³) === LDA
    @test FunctionalCategory(UnitfulHartree.∂³ϵ_∂ρ²∂σ) === GGA
    @test FunctionalCategory(UnitfulHartree.∂³ϵ_∂ρ∂σ²) === GGA
    @test FunctionalCategory(UnitfulHartree.∂³ϵ_∂σ³) === GGA
end
