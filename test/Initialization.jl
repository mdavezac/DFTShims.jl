using DFTShims
using AxisArrays
using Unitful

const DH = DFTShims.Dispatch.Hartree
const Dρ = DH.Scalars.ρ
const SIZES = 10, 2
const AXES = Axis{:radius}(1:SIZES[1]), Axis{:bb}((:α, :β))
const polarized = ColinearSpin()
const unpolarized = SpinDegenerate()

@testset "SpinDegenerate" begin
    ρ = zeros(Dρ{Int64}, false, SIZES...)
    @test typeof(ρ) <: AxisArray{Dρ{Int64}}
    @test size(ρ) == SIZES

    ρ = zeros(Dρ{Int64}, unpolarized, SIZES[1], AXES[1], SIZES[2])
    @test typeof(ρ) <: AxisArray{Dρ{Int64}}
    @test size(ρ) == SIZES
    @test axes(ρ, 1) == AXES[1]

    ρ = zeros(Dρ{Int64}, unpolarized, SIZES..., AXES[1], AXES[2])
    @test typeof(ρ) <: AxisArray{Dρ{Int64}}
    @test size(ρ) == SIZES
    @test axes(ρ, 1) == AXES[1]
    @test axes(ρ, 2) == AXES[2]

    @test_throws ArgumentError zeros(Dρ{Int64}, unpolarized, SIZES..., AXES..., AXES[1])

    @inferred zeros(Dρ{Int64}, unpolarized, SIZES)
    @inferred zeros(Dρ{Int64}, unpolarized, SIZES, AXES)
    @inferred zeros(Dρ{Int64}, unpolarized, SIZES, AXES[1:1])
end

@testset "Axis Manipulations" begin
    a = AxisArray(zeros(2, 3))
    const add_spin_axis = DFTShims.ArrayInitialization.add_spin_axis
    const replace_spin_axis = DFTShims.ArrayInitialization.replace_spin_axis
    @test @inferred(add_spin_axis(polarized, (2, 3), 5)) == (2, 3, 5)
    @test @inferred(add_spin_axis(polarized, AXES, AXES[1])) == (AXES..., AXES[1])
    @test @inferred(add_spin_axis(ColinearSpinFirst(), AXES, AXES[1])) ==
        (AXES[1], AXES...)

    saxes = AXES..., Axis{:spin}((:u, :d))
    actual = @inferred replace_spin_axis(Dρ, polarized, (4, 3, 2), saxes)
    @test actual == ((4, 3, 2), (AXES..., Axis{:spin}((:α, :β))))
    actual = @inferred replace_spin_axis(Dρ, ColinearSpinFirst(), (4, 3, 2), saxes)
    @test actual == ((2, 4, 3), (Axis{:spin}((:α, :β)), AXES...))
    actual = @inferred replace_spin_axis(Dρ, ColinearSpinLast(), (4, 3, 2), saxes)
    @test actual == ((4, 3, 2), (AXES..., Axis{:spin}((:α, :β))))

    saxes = Axis{:spin}((:u, :d)), AXES... 
    actual = @inferred replace_spin_axis(Dρ, polarized, (2, 3, 4), saxes)
    @test actual == ((2, 3, 4), (Axis{:spin}((:α, :β)), AXES...))
    actual = @inferred replace_spin_axis(Dρ, ColinearSpinFirst(), (2, 3, 4), saxes)
    @test actual == ((2, 3, 4), (Axis{:spin}((:α, :β)), AXES...))
    actual = @inferred replace_spin_axis(Dρ, ColinearSpinLast(), (2, 3, 4), saxes)
    @test actual == ((3, 4, 2), (AXES..., Axis{:spin}((:α, :β))))
end

@testset "ColinearSpin" begin
    @testset "axis names" begin
        comps = components(Dρ{Int64}, polarized)
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

    ρ = zeros(Dρ{Int64}, polarized, SIZES..., AXES...)
    @test typeof(ρ) <: AxisArray{Dρ{Int64}}
    @test size(ρ) == (SIZES..., 2)
    @test axes(ρ, 1) == AXES[1]
    @test axes(ρ, 2) == AXES[2]
    @test axes(ρ, 3) == Axis{:spin}((:α, :β))

    ρ = zeros(Dρ{Int64}, polarized, SIZES..., AXES..., Axis{:aa}((1, 2)))
    @test typeof(ρ) <: AxisArray{Dρ{Int64}}
    @test size(ρ) == (SIZES..., 2)
    @test axes(ρ, 1) == AXES[1]
    @test axes(ρ, 2) == AXES[2]
    @test axes(ρ, 3) == Axis{:aa}((1, 2))

    const D∂³ϵ_∂ρ∂σ² = DH.Scalars.∂³ϵ_∂ρ∂σ²
    ∂³ϵ_∂ρ∂σ² = zeros(D∂³ϵ_∂ρ∂σ²{Int16}, polarized, SIZES..., AXES...)
    @test size(∂³ϵ_∂ρ∂σ²) == (SIZES..., 12)
    @test axes(∂³ϵ_∂ρ∂σ², 1) == AXES[1]
    @test axes(∂³ϵ_∂ρ∂σ², 2) == AXES[2]

    @inferred zeros(Dρ{Int64}, polarized, SIZES)
    @inferred zeros(Dρ{Int64}, polarized, SIZES, AXES)
    @inferred zeros(Dρ{Int64}, polarized, SIZES, AXES[1:1])
end
