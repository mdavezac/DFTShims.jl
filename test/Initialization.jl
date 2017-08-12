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

    @test_throws(ArgumentError,
                 zeros(Dρ{Int64}, polarized, SIZES..., AXES..., Axis{:aa}((1, 2))))

    const D∂³ϵ_∂ρ∂σ² = DH.Scalars.∂³ϵ_∂ρ∂σ²
    ∂³ϵ_∂ρ∂σ² = zeros(D∂³ϵ_∂ρ∂σ²{Int16}, polarized, SIZES..., AXES...)
    @test size(∂³ϵ_∂ρ∂σ²) == (SIZES..., 12)
    @test axes(∂³ϵ_∂ρ∂σ², 1) == AXES[1]
    @test axes(∂³ϵ_∂ρ∂σ², 2) == AXES[2]

    @inferred zeros(Dρ{Int64}, polarized, SIZES)
    @inferred zeros(Dρ{Int64}, polarized, SIZES, AXES)
    @inferred zeros(Dρ{Int64}, polarized, SIZES, AXES[1:1])
end

@testset "From template array" begin
    ρ = zeros(Dρ{Int32}, true, SIZES..., AXES...)
    @test all(zeros(DH.Scalars.∂²ϵ_∂σ²{Int64}, ρ) .== 0u"∂²ϵ_∂σ²")
    @test eltype(zeros(DH.Scalars.∂²ϵ_∂σ²{Int64}, ρ)) == typeof(0u"∂²ϵ_∂σ²")
    @test is_spin_polarized(zeros(DH.Scalars.∂²ϵ_∂σ²{Int64}, ρ))
    @test length(zeros(DH.Scalars.∂²ϵ_∂σ²{Int64}, ρ)[Axis{:spin}]) == 6
    @test find(x -> typeof(x) <: Axis{:spin}, 
               axes(zeros(DH.Scalars.∂²ϵ_∂σ²{Int64}, ρ))) == [length(AXES) + 1]
    @test axes(zeros(DH.Scalars.∂²ϵ_∂σ²{Int64}, ρ))[1:end - 1] == AXES
    @test axes(zeros(DH.Scalars.∂²ϵ_∂σ²{Int64}, ColinearSpinFirst(), ρ))[2:end] == AXES

    @test !is_spin_polarized(zeros(DH.Scalars.∂²ϵ_∂σ²{Int64}, SpinDegenerate(), ρ))
    @test ndims(zeros(DH.Scalars.∂²ϵ_∂σ²{Int64}, SpinDegenerate(), ρ)) == length(SIZES)
    @test axes(zeros(DH.Scalars.∂²ϵ_∂σ²{Int64}, SpinDegenerate(), ρ)) == AXES

    ρ = zeros(Dρ{Int64}, SpinDegenerate(), ρ)
    @test axes(zeros(DH.Scalars.∂²ϵ_∂σ²{Int64}, ColinearSpinFirst(), ρ))[2:end] == AXES
    @test is_spin_polarized(zeros(DH.Scalars.∂²ϵ_∂σ²{Int64}, ColinearSpinFirst(), ρ))
end

@testset "Permute spin axis" begin
    ρₙ = zeros(Dρ{Int32}, true, SIZES..., AXES...)
    ρₙ[:] = (1:length(ρₙ)) * oneunit(eltype(ρₙ))
    @test permutedims(ρₙ, ColinearSpinLast()) === ρₙ
    @test permutedims(ρₙ, ColinearSpinPreferLast()) === ρₙ
    ρ₀  = permutedims(ρₙ, ColinearSpinFirst())
    @test axes(ρ₀) == (axes(ρₙ, Axis{:spin}), AXES...)

    ρ = permutedims(ρ₀, ColinearSpinLast())
    @test axes(ρ) == (AXES..., axes(ρₙ, Axis{:spin}))
    @test permutedims(ρ₀, ColinearSpinPreferLast()) === ρ₀
end
