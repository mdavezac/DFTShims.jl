using DFTShims: ColinearSpin, SpinDegenerate, is_spin_polarized
using AxisArrays
using Unitful
import AxisArrays: axes

const DH = DFTShims.Dispatch.Hartree
const Dρ = DH.Scalars.ρ
const Dϵ = DH.Scalars.ϵ
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

    ρ = @inferred reinterpret(Dρ{Int64}, SpinDegenerate(), [1 1 1; 2 2 2])
    @test typeof(ρ) <: AxisArray
    @test eltype(ρ) == Dρ{Int64}
    @test !is_spin_polarized(ρ)

    @test_throws ArgumentError zeros(Dρ{Int64}, unpolarized, SIZES..., AXES..., AXES[1])

    @inferred zeros(Dρ{Int64}, unpolarized, SIZES)
    @inferred zeros(Dρ{Int64}, unpolarized, SIZES, AXES)
    @inferred zeros(Dρ{Int64}, unpolarized, SIZES, AXES[1:1])

    ϵ = @inferred zeros(Dϵ{Int64}, ColinearSpin(), SIZES)
    @test is_spin_polarized(ϵ)
    @test :spin ∈ axisnames(ϵ)
    @test size(ϵ) == (SIZES..., 2)

    @test @inferred(SpinCategory(ρ)) == SpinDegenerate()
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
    @test @inferred(SpinCategory(ρ)) == ColinearSpin()

    ρ = @inferred zeros(Dρ{Int64}, polarized, AXES)
    @test typeof(ρ) <: AxisArray{Dρ{Int64}}
    @test size(ρ) == (SIZES..., 2)
    @test axes(ρ, 1) == AXES[1]
    @test axes(ρ, 2) == AXES[2]
    @test axes(ρ, 3) == Axis{:spin}((:α, :β))

    @test_throws(ArgumentError,
                 zeros(Dρ{Int64}, polarized, SIZES..., AXES..., Axis{:aa}((1, 2))))

    const D∂³ϵ_∂ρ∂σ² = DH.Scalars.∂³ϵ_∂ρ∂σ²
    ∂³ϵ_∂ρ∂σ² = @inferred zeros(D∂³ϵ_∂ρ∂σ²{Int16}, polarized, AXES)
    @test size(∂³ϵ_∂ρ∂σ²) == (SIZES..., 12)
    @test axes(∂³ϵ_∂ρ∂σ², 1) == AXES[1]
    @test axes(∂³ϵ_∂ρ∂σ², 2) == AXES[2]

    ρ = @inferred reinterpret(Dρ{Int64}, ColinearSpinFirst(), [1 1 1; 2 2 2])
    @test typeof(ρ) <: AxisArray
    @test eltype(ρ) == Dρ{Int64}
    @test is_spin_polarized(ρ)
    @test @inferred(SpinCategory(ρ)) == ColinearSpinFirst()

    ∂ϵ_∂ρ = @inferred rand(typeof(1u"eV*nm^3"), ColinearSpin{length(AXES)}(), AXES)
    @test axes(∂ϵ_∂ρ, length(AXES)) isa Axis{:spin}
    @test ndims(∂ϵ_∂ρ) == length(AXES) + 1

    @test_throws(ArgumentError,
                 reinterpret(Dρ{Int64}, ColinearSpinLast(), [1 1 1; 2 2 2]))

    @inferred zeros(Dρ{Int64}, polarized, SIZES)
    @inferred zeros(Dρ{Int64}, polarized, AXES)
    @inferred rand(typeof(1u"eV*nm^3"), ColinearSpin{length(AXES)}(), AXES)
    @inferred zeros(Dρ{Int64}, polarized, SIZES, AXES)
    @inferred zeros(Dρ{Int64}, polarized, SIZES, AXES[1:1])
    @inferred zeros(Dρ{Int64}, polarized, SIZES, AXES[1:1])
end

@testset "From template array" begin
    ρ = zeros(Dρ{Int32}, true, SIZES..., AXES...)
    @test all(zeros(ρ, DH.Scalars.∂²ϵ_∂σ², SpinAware()) .== 0u"∂²ϵ_∂σ²")
    @test !(typeof(zeros(ρ, DH.Scalars.∂²ϵ_∂σ², SpinAware()).data) <: AxisArray)
    @test eltype(zeros(ρ, DD.Scalars.∂²ϵ_∂σ², SpinAware())) == typeof(Int32(0)u"∂²ϵ_∂σ²")
    @test eltype(zeros(ρ, DH.Scalars.∂²ϵ_∂σ²{Int64}, SpinAware())) == typeof(0u"∂²ϵ_∂σ²")
    @test is_spin_polarized(zeros(ρ, DH.Scalars.∂²ϵ_∂σ²{Int64}, SpinAware()))
    @test length(zeros(ρ, DH.Scalars.∂²ϵ_∂σ²{Int64}, SpinAware())[Axis{:spin}]) == 6
    @test find(x -> typeof(x) <: Axis{:spin},
               axes(zeros(ρ, DH.Scalars.∂²ϵ_∂σ²{Int64}, SpinAware()))) == [length(AXES) + 1]
    @test axes(zeros(ρ, DH.Scalars.∂²ϵ_∂σ²{Int64}, SpinAware()))[1:end - 1] == AXES
    @test axes(zeros(ρ, DH.Scalars.∂²ϵ_∂σ², ColinearSpinFirst()))[2:end] == AXES

    @test !is_spin_polarized(zeros(ρ, DH.Scalars.∂²ϵ_∂σ²{Int64}, SpinDegenerate()))
    @test !(typeof(zeros(DH.Scalars.∂²ϵ_∂σ²{Int64}, SpinDegenerate()).data) <: AxisArray)
    @test ndims(zeros(ρ, DH.Scalars.∂²ϵ_∂σ²{Int64}, SpinDegenerate())) == length(SIZES)
    @test axes(zeros(ρ, DH.Scalars.∂²ϵ_∂σ²{Int64}, SpinDegenerate())) == AXES

    ϵ = @inferred zeros(ρ, Dϵ, SpinDegenerate())
    @test !is_spin_polarized(ϵ)
    @test :spin ∉ axisnames(ϵ)
    @test size(ϵ) == SIZES
    @test !(typeof(ϵ.data) <: AxisArray)

    ρ = zeros(ρ, Dρ{Int64}, SpinDegenerate())
    @test axes(zeros(ρ, DD.Scalars.∂²ϵ_∂σ², ColinearSpinFirst()))[2:end] == AXES
    @test is_spin_polarized(zeros(ρ, DH.Scalars.∂²ϵ_∂σ², ColinearSpinFirst()))

    ϵ = @inferred zeros(ρ, Dϵ, ColinearSpin())
    @test is_spin_polarized(ϵ)
    @test :spin ∈ axisnames(ϵ)
    @test size(ϵ) == (SIZES..., 2)

    @inferred zeros(ρ, DH.Scalars.∂²ϵ_∂σ², SpinAware())
    @inferred zeros(ρ, DH.Scalars.∂²ϵ_∂σ²{Int32}, SpinAware())
    @inferred zeros(ρ, DD.Scalars.∂²ϵ_∂σ², ColinearSpin())
end

@testset "Convert between arrays" begin
    ρₙ = zeros(Dρ{Int32}, true, SIZES..., AXES...)
    ρₙ[:] = (1:length(ρₙ)) * oneunit(eltype(ρₙ))
    @test @inferred(convert(ColinearSpinLast(), ρₙ)) === ρₙ
    @test !(typeof(convert(ColinearSpinLast(), ρₙ).data) <: AxisArray)
    ρ₀  = @inferred convert(ColinearSpinFirst(), ρₙ)
    @test axes(ρ₀) == (axes(ρₙ, Axis{:spin}), AXES...)
    @test !(typeof(ρ₀.data) <: AxisArray)

    ρ = @inferred convert(ColinearSpinLast(), ρ₀)
    @test axes(ρ) == (AXES..., axes(ρₙ, Axis{:spin}))
    @test all(ρ .== ρₙ)
    @test @inferred(convert(ColinearSpinLast(), ρ)) === ρ

    @test_throws InexactError convert(typeof(1u"m^-3"), ColinearSpinFirst(), ρₙ)
    ρₘ = @inferred convert(typeof(1.0u"m^-3"), ColinearSpinFirst(), ρₙ)
    @test axes(ρₘ) == axes(ρ₀)
    @test ρₘ ≈ ρ₀
    @test !(ustrip.(ρₘ) ≈ ustrip.(ρ₀))

    ϵ = similar(ρ₀, DH.Scalars.ϵ, SpinAware())
    @test is_spin_polarized(ϵ)
    @test !(typeof(ϵ.data) <: AxisArray)
end

@testset "Creation via keyword arguments" begin
    axis = 0u"a₀":0.1u"a₀":4u"a₀"
    ρ = zeros(DH.Scalars.ρ{Float64}, false, radius=axis)
    @test eltype(ρ) == DH.Scalars.ρ{Float64}
    @test size(ρ) == (length(axis), )
    @test axes(ρ, 1) == Axis{:radius}(axis)

    ρ = rand(DH.Scalars.ρ{Float64}, false, radius=axis)
    @test eltype(ρ) == DH.Scalars.ρ{Float64}
    @test size(ρ) == (length(axis), )
    @test axes(ρ, 1) == Axis{:radius}(axis)
end

@testset "Wrap a dimensionless array" begin
    @test is_spin_polarized(wrap([1, 2, 3]u"m^-3")) == false
    @test eltype(wrap([1, 2, 3]u"m^-3")) <: DD.Scalars.ρ
    @test_throws ArgumentError wrap(ColinearSpin(), [1, 2, 3]u"m^-3")
    @test is_spin_polarized(wrap(ColinearSpinFirst(), [1 2 3; 4 5 6]u"m^-3"))
    @test is_spin_polarized(wrap(DD.Scalars.ρ, ColinearSpinFirst(), [1 2 3; 4 5 6]))
    @test_throws ArgumentError wrap(DD.Scalars.ρ, ColinearSpin(), [1 2 3; 4 5 6])
end
