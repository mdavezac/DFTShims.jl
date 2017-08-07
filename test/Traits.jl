using DFTShims
using AxisArrays
using Unitful

spin = AxisArray(zeros((10, 2)), Axis{:radius}(0u"m":1u"m":9u"m"), Axis{:spin}((:+, :-)))
nospin = AxisArray(zeros((10, 2)), Axis{:radius}(0u"m":1u"m":9u"m"), Axis{:n}((1, 2)))
const polarized = ColinearSpin()
const unpolarized = SpinDegenerate()

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
    @test @inferred is_spin_polarized(spin)
    @test ! @inferred is_spin_polarized(nospin)

    @test @inferred is_spin_polarized(typeof(spin))
    @test ! @inferred is_spin_polarized(typeof(nospin))

    @test @inferred(SpinCategory(spin)) === polarized
    @test @inferred(SpinCategory(nospin)) === unpolarized
end

@testset "components" begin
    @test @inferred(components(UnitfulHartree.ρ, unpolarized)) == (:ρ,)
    @test all(components(UnitfulHartree.ρ, polarized) .== (:α, :β))
    @test length(components(Dispatch.Dimensions.Scalars.∂ϵ_∂σ, unpolarized)) == 1
    @test length(@inferred(components(UnitfulHartree.∂ϵ_∂σ, polarized))) == 3
    @test length(@inferred(components(UnitfulHartree.∂²ϵ_∂ρ², polarized))) == 3
    @test length(@inferred(components(UnitfulHartree.∂²ϵ_∂ρ∂σ, unpolarized))) == 1
    @test length(@inferred(components(UnitfulHartree.∂²ϵ_∂ρ∂σ, polarized))) == 6
    @test length(@inferred(components(UnitfulHartree.∂²ϵ_∂σ², polarized))) == 6
    @test length(@inferred(components(UnitfulHartree.∂³ϵ_∂ρ³, polarized))) == 4
    @test length(@inferred(components(UnitfulHartree.∂³ϵ_∂ρ²∂σ, polarized))) == 9
    @test length(@inferred(components(UnitfulHartree.∂³ϵ_∂ρ∂σ², polarized))) == 12
    @test length(@inferred(components(UnitfulHartree.∂³ϵ_∂σ³, polarized))) == 10
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
