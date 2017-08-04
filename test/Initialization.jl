using DFTShims
using AxisArrays
using Unitful

const DH = DFTShims.Dispatch.Hartree
const Dρ = DH.Scalars.ρ
const SIZES = 10, 2
const AXES = Axis{:radius}(1:SIZES[1]), Axis{:bb}((:α, :β))

@testset "Unpolarized" begin
    ρ = zeros(Dρ{Int64}, false, SIZES...)
    @test typeof(ρ) <: AxisArray{Dρ{Int64}}
    @test size(ρ) == SIZES

    ρ = zeros(Dρ{Int64}, Unpolarized, SIZES[1], AXES[1], SIZES[2])
    @test typeof(ρ) <: AxisArray{Dρ{Int64}}
    @test size(ρ) == SIZES
    @test axes(ρ, 1) == AXES[1]

    ρ = zeros(Dρ{Int64}, Unpolarized, SIZES..., AXES[1], AXES[2])
    @test typeof(ρ) <: AxisArray{Dρ{Int64}}
    @test size(ρ) == SIZES
    @test axes(ρ, 1) == AXES[1]
    @test axes(ρ, 2) == AXES[2]

    @test_throws ArgumentError zeros(Dρ{Int64}, Unpolarized, SIZES..., AXES..., AXES[1])

    @inferred zeros(Dρ{Int64}, Unpolarized, SIZES)
    @inferred zeros(Dρ{Int64}, Unpolarized, SIZES, AXES)
    @inferred zeros(Dρ{Int64}, Unpolarized, SIZES, AXES[1:1])
end

@testset "Polarized" begin
    @testset "axis names" begin
        comps = components(Dρ{Int64}, Polarized)
        extra_axis = Axis{:a}((1, 2))
        spin_axis = Axis{:spin}((:α, :β)) 
        defaults(a) = AxisArrays.default_axes(zeros((SIZES..., 2)), a)
        pa(a) = @inferred DFTShims.ArrayInitialization.polarized_axis(comps, defaults(a), a)

        @test pa((AXES..., extra_axis)) == defaults((AXES..., extra_axis))
        @test pa(AXES) == (AXES..., spin_axis)
        @test pa(AXES[1:1]) == (defaults(AXES[1:1])[1:end - 1]..., spin_axis)
        @test pa(AXES[1:0]) == (defaults(AXES[1:0])[1:end - 1]..., spin_axis)
    end

    ρ = zeros(Dρ{Int64}, true, SIZES...)
    @test typeof(ρ) <: AxisArray{Dρ{Int64}}
    @test size(ρ) == (SIZES..., 2)
    @test axes(ρ, 3) == Axis{:spin}((:α, :β))

    ρ = zeros(Dρ{Int64}, Polarized, SIZES..., AXES...)
    @test typeof(ρ) <: AxisArray{Dρ{Int64}}
    @test size(ρ) == (SIZES..., 2)
    @test axes(ρ, 1) == AXES[1]
    @test axes(ρ, 2) == AXES[2]
    @test axes(ρ, 3) == Axis{:spin}((:α, :β))

    ρ = zeros(Dρ{Int64}, Polarized, SIZES..., AXES..., Axis{:aa}((1, 2)))
    @test typeof(ρ) <: AxisArray{Dρ{Int64}}
    @test size(ρ) == (SIZES..., 2)
    @test axes(ρ, 1) == AXES[1]
    @test axes(ρ, 2) == AXES[2]
    @test axes(ρ, 3) == Axis{:aa}((1, 2))

    const D∂³ϵ_∂ρ∂σ² = DH.Scalars.∂³ϵ_∂ρ∂σ²
    ∂³ϵ_∂ρ∂σ² = zeros(D∂³ϵ_∂ρ∂σ²{Int16}, Polarized, SIZES..., AXES...)
    @test size(∂³ϵ_∂ρ∂σ²) == (SIZES..., 12)
    @test axes(∂³ϵ_∂ρ∂σ², 1) == AXES[1]
    @test axes(∂³ϵ_∂ρ∂σ², 2) == AXES[2]

    @inferred zeros(Dρ{Int64}, Polarized, SIZES)
    @inferred zeros(Dρ{Int64}, Polarized, SIZES, AXES)
    @inferred zeros(Dρ{Int64}, Polarized, SIZES, AXES[1:1])
end
