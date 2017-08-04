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

    @test Polarization(spin) === Polarized
    @test Polarization(nospin) === Unpolarized
end

@testset "components" begin
    @test components(UnitfulHartree.ρ, Unpolarized) == (:ρ,)
    @test all(components(UnitfulHartree.ρ, Polarized) .== (:α, :β))
    @test length(components(Dispatch.Dimensions.Scalars.∂ϵ_∂σ, Unpolarized)) == 1
    @test length(components(UnitfulHartree.∂ϵ_∂σ, Polarized)) == 3
    @test length(components(UnitfulHartree.∂²ϵ_∂ρ², Polarized)) == 3
    @test length(components(UnitfulHartree.∂²ϵ_∂ρ∂σ, Unpolarized)) == 1
    @test length(components(UnitfulHartree.∂²ϵ_∂ρ∂σ, Polarized)) == 6
    @test length(components(UnitfulHartree.∂²ϵ_∂σ², Polarized)) == 6
    @test length(components(UnitfulHartree.∂³ϵ_∂ρ³, Polarized)) == 4
    @test length(components(UnitfulHartree.∂³ϵ_∂ρ²∂σ, Polarized)) == 9
    @test length(components(UnitfulHartree.∂³ϵ_∂ρ∂σ², Polarized)) == 12
    @test length(components(UnitfulHartree.∂³ϵ_∂σ³, Polarized)) == 10
end

@testset "functional category" begin
    @test functional_category(UnitfulHartree.ρ) === LDA
    @test functional_category(UnitfulHartree.∂ϵ_∂σ) === GGA
    @test functional_category(UnitfulHartree.∂²ϵ_∂ρ²) === LDA
    @test functional_category(UnitfulHartree.∂²ϵ_∂ρ∂σ) === GGA
    @test functional_category(Dispatch.Dimensions.Scalars.∂²ϵ_∂σ²) === GGA
    @test functional_category(UnitfulHartree.∂³ϵ_∂ρ³) === LDA
    @test functional_category(UnitfulHartree.∂³ϵ_∂ρ²∂σ) === GGA
    @test functional_category(UnitfulHartree.∂³ϵ_∂ρ∂σ²) === GGA
    @test functional_category(UnitfulHartree.∂³ϵ_∂σ³) === GGA
end
