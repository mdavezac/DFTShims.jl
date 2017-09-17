module Traits
export has_axis, is_spin_polarized, ColinearSpin, SpinDegenerate, components,
       SpinCategory, FunctionalCategory, LDA, GGA, ColinearSpinFirst, ColinearSpinLast,
       concretize_type, standard_name

using DocStringExtensions
using AxisArrays
using Unitful
using ..UnitfulHartree
using ..Dispatch
const DD = Dispatch.Dimensions
const DH = Dispatch.Hartree

macro lintpragma(s) end
@lintpragma("Ignore use of undeclared variable unitful_dimensions")
@lintpragma("Ignore use of undeclared variable Q")
@lintpragma("Ignore use of undeclared variable T")
@lintpragma("Ignore use of undeclared variable D")
@lintpragma("Ignore use of undeclared variable U")
@lintpragma("Ignore use of undeclared variable H")
@lintpragma("Ignore unused array")

const UH = UnitfulHartree

has_axis(array::Type{<:AxisArray}, axis::Type{<:Axis}) =
    axisnames(axis)[1] ∈ axisnames(array)
has_axis(a::AxisArray, axis_type::Type{<:Axis}) = has_axis(typeof(a), axis_type)
has_axis(array::Type{<:AxisArray}, name::Symbol) = has_axis(array, Axis{name})
has_axis(a::AxisArray, name::Symbol) = has_axis(typeof(a), name)
    

""" Spin axis should have this type """
const SpinAxis = Axis{:spin, Q} where {Q <: Tuple{T, T, Vararg{T}} where T}
"""
True if the array is spin polarized, meaning it has an axis named `:spin` with
two or more values.

$(SIGNATURES)

Note that if the argument is a type (rather than an instance) then the values of the spin
axis should be a tuple, e.g.:

`typeof(Axis(:spin, (:+, :-))) === Axis{:spin, Tuple{Symbol, Symbol}}`

More explicitly, a spin axis is an `AxisArrays.Axis` instance of type
`Axis{:spin, <: Tuple{T, T, Vararg{T}} where T}`.
"""
is_spin_polarized(::Type{<: SpinAxis}) = true
is_spin_polarized(::Type{<: Axis}) = false
is_spin_polarized(ax::Axis) = is_spin_polarized(typeof(ax))
is_spin_polarized(array::AxisArray) = any(is_spin_polarized.(axes(array)))
is_spin_polarized(array::Type{<: AxisArray}) = findfirst(is_spin_polarized, array) ≠ 0
Base.findfirst(::typeof(is_spin_polarized), array::Type{<: AxisArray}) = begin
    index = 1
    for T in array.parameters[end].parameters
        is_spin_polarized(T) && return index
        index += 1
    end
    0
end

""" Union of all polarization traits """
abstract type SpinCategory end
""" Trait for functions accepting input in the colinear spin approximation """
abstract type AbstractColinearSpin <: SpinCategory end
""" Trait for functions accepting input in the colinear spin approximation """
struct ColinearSpin{L} <: AbstractColinearSpin end

""" Spin-polarized input where the spin axis is always the fastest changing """
const ColinearSpinFirst = ColinearSpin{1}
""" Spin-polarized input where the spin axis is always the slowest changing """
const ColinearSpinLast = ColinearSpin{:end}

"""
Trait to differentiate between Base and specialized spin functions

This type is used to select overloads of `Base.zeros` and friends that will change the
spin-axis if needed. More specifically, if `ρ` is a spin-polarized `AxisArray`, then
`Base.zeros(ρ, Dispatch.Scalars.∂³ϵ_∂ρ³, SpinAware)` is also spin-polarized and the
spin-axis has the correct dimension (4 in this case). If `ρ` is not spin-polarized, then
the call is equivalent to `Base.zeros(ρ, Dispatch.Scalars.∂³ϵ_∂ρ³)`.
"""
struct SpinAware <: SpinCategory end

""" Traits for unpolarized inputs """
struct SpinDegenerate <: SpinCategory end

""" Looks at axes to figure out the spin category

- no spin-axis -> SpinDegenerate
- fastest changing dimension -> ColinearSpinFirst
- `n`th dim -> ColinearSpin{n}

"""
(::Type{SpinCategory})(array::AxisArray) = SpinCategory(axes(array))
@generated (::Type{SpinCategory})(ax::Tuple{Axis, Vararg{Axis}}) = begin
    index = findfirst(is_spin_polarized, ax.parameters)
    if index == 0
        :(SpinDegenerate())
    elseif index == length(ax.parameters)
        :(ColinearSpinLast())
    else
        :(ColinearSpin{$index}())
    end
end
(::Type{SpinCategory})(S::SpinCategory) = S
""" ColinearSpin defaults to settings spin axis in last position

However, if the spin axis is already given, functions following this traits should not move
it.
"""
(::Type{ColinearSpin})() = ColinearSpin{:end}()

""" Supertype of all functional categories """
abstract type FunctionalCategory end
""" Trait identifying the LDA functional category """
type LDA <: FunctionalCategory end
""" Trait identifying the LDA functional category """
type GGA <: FunctionalCategory end

"""
Same as Unitful.dimension but for still abstract quantities
"""
unitful_dimensions(::Type{<: Quantity{T, D} where T}) where D = D()

"""
Figures out functional category, whether LDA or GGA

$(SIGNATURES)
"""
(::Type{FunctionalCategory})(::typeof(dimension(UH.ρ)))::LDA = LDA()
(::Type{FunctionalCategory})(::typeof(dimension(UH.∂ϵ_∂ρ))) = LDA()
(::Type{FunctionalCategory})(::typeof(dimension(UH.∂²ϵ_∂ρ²))) = LDA()
(::Type{FunctionalCategory})(::typeof(dimension(UH.∂³ϵ_∂ρ³))) = LDA()
(::Type{FunctionalCategory})(::typeof(dimension(UH.σ))) = GGA()
(::Type{FunctionalCategory})(::typeof(dimension(UH.∂ϵ_∂σ))) = GGA()
(::Type{FunctionalCategory})(::typeof(dimension(UH.∂²ϵ_∂ρ∂σ))) = GGA()
(::Type{FunctionalCategory})(::typeof(dimension(UH.∂²ϵ_∂σ²))) = GGA()
(::Type{FunctionalCategory})(::typeof(dimension(UH.∂³ϵ_∂ρ²∂σ))) = GGA()
(::Type{FunctionalCategory})(::typeof(dimension(UH.∂³ϵ_∂ρ∂σ²))) = GGA()
(::Type{FunctionalCategory})(::typeof(dimension(UH.∂³ϵ_∂σ³))) = GGA()
(::Type{FunctionalCategory})(u::Unitful.FreeUnits) = FunctionalCategory(dimension(u))
(::Type{FunctionalCategory})(u::DD.Scalars.All) = FunctionalCategory(dimension(u))
(::Type{FunctionalCategory})(T::Type{<: DD.Scalars.All}) =
    FunctionalCategory(unitful_dimensions(T))

"""
Labels of the standard components for LDA and GGA inputs

$(SIGNATURES)

ρ refers to the unpolarized density, α and β to the two spin channels of the density, σ to
the contracted gradient density (∇ρ⋅∇ρ), and σαα (∇α⋅∇α), σαβ, σββ  to the polarized
contracted gradient densities.
"""
components(::typeof(dimension(UH.ρ)), ::SpinDegenerate) = (:ρ,)
components(::typeof(dimension(UH.ϵ)), ::SpinDegenerate) = (:ϵ,)
components(::typeof(dimension(UH.∂ϵ_∂ρ)), ::SpinDegenerate) = (:∂ρ,)
components(::typeof(dimension(UH.∂²ϵ_∂ρ²)), ::SpinDegenerate) = (:∂ρ²,)
components(::typeof(dimension(UH.∂³ϵ_∂ρ³)), ::SpinDegenerate) = (:∂ρ³,)
components(::typeof(dimension(UH.ρ)), ::ColinearSpin) = :α, :β
components(::typeof(dimension(UH.ϵ)), ::ColinearSpin) = :α, :β
components(::typeof(dimension(UH.∂ϵ_∂ρ)), ::ColinearSpin) = :∂α, :∂β
components(::typeof(dimension(UH.∂²ϵ_∂ρ²)), ::ColinearSpin) = :∂α², :∂α∂β, :∂²β 
components(::typeof(dimension(UH.∂³ϵ_∂ρ³)), ::ColinearSpin) = :∂α³, :∂α∂β², :∂α∂β², :∂β³
components(::typeof(dimension(UH.σ)), ::ColinearSpin) = :σαα, :σαβ, :σββ
components(::typeof(dimension(UH.∂ϵ_∂σ)), ::ColinearSpin) = :∂σαα, :∂σαβ, :∂σββ
components(::typeof(dimension(UH.∂²ϵ_∂ρ∂σ)), ::ColinearSpin) = 
    :∂α∂σαα, :∂α∂σαβ, :∂α∂σββ, :∂β∂σαα, :∂β∂σαβ, :∂β∂σββ
components(::typeof(dimension(UH.∂²ϵ_∂σ²)), ::ColinearSpin) =
    :∂σαα², :∂σαα∂σαβ, :∂σαα∂σββ, :∂σαβ², :∂σαβσββ, :∂σββ² 
components(::typeof(dimension(UH.∂³ϵ_∂ρ²∂σ)), ::ColinearSpin) = (
    :∂α²∂σαα, :∂α²∂σαβ, :∂α²∂σββ,
    :∂α∂β∂σαα, :∂α∂β∂σαβ, :∂α∂β∂σββ,
    :∂β²∂σαα, :∂β²∂σαβ, :∂β²∂σββ
)
components(::typeof(dimension(UH.∂³ϵ_∂ρ∂σ²)), ::ColinearSpin) = (
    :∂α∂σαα², :∂α∂σαα∂σαβ, :∂α∂σαα∂σββ, :∂α∂σαβ², :∂α∂σαβσββ, :∂α∂σββ²,
    :∂β∂σαα², :∂β∂σαα∂σαβ, :∂β∂σαα∂σββ, :∂β∂σαβ², :∂β∂σαβσββ, :∂β∂σββ² 
)
components(::typeof(dimension(UH.∂³ϵ_∂σ³)), ::ColinearSpin) = (
    :∂σαα³, :∂σαα²∂σαβ, :∂σαα²∂σββ, :∂σαα∂σαβ², :∂σαα∂σαβ∂σββ, :∂σαα∂σββ², 
    :∂σαβ³, :∂σαβ²∂σββ, :∂σαβ∂σββ², :∂σββ³
)
components(::typeof(dimension(UH.σ)), ::SpinDegenerate) = (:σ,)
components(::typeof(dimension(UH.∂ϵ_∂σ)), ::SpinDegenerate) = (:∂σ,)
components(::typeof(dimension(UH.∂²ϵ_∂ρ∂σ)), ::SpinDegenerate) = (:∂ρ∂σ,)
components(::typeof(dimension(UH.∂²ϵ_∂σ²)), ::SpinDegenerate) = (:∂σ²,)
components(::typeof(dimension(UH.∂³ϵ_∂ρ²∂σ)), ::SpinDegenerate) = (:∂ρ²∂σ,)
components(::typeof(dimension(UH.∂³ϵ_∂ρ∂σ²)), ::SpinDegenerate) = (:∂ρ∂σ²,)
components(::typeof(dimension(UH.∂³ϵ_∂σ³)), ::SpinDegenerate) = (:∂σ³,)

components(u::Unitful.FreeUnits, P::SpinCategory) = components(dimension(u), P)
components(u::DD.Scalars.All, P::SpinCategory) = components(dimension(u), P)
components(T::Type{<: DD.Scalars.All}, P::SpinCategory) =
    components(unitful_dimensions(T), P)

""" helper function to concretize some types using Hartree units """
function hartree_concretize_type end
for (_, name) in Dispatch.Dimensioned
    @eval begin
        @inline hartree_concretize_type(T::Type{<:DD.Scalars.$name}) =
            hartree_concretize_type(DH.Scalars.$name, T)
        @inline hartree_concretize_type(::Type{DD.Scalars.$name},
                                        ::Type{<: Quantity{T}}) where T =
            DH.Scalars.$name{T}
        @inline hartree_concretize_type(::Type{DD.Scalars.$name}, T::Type{<: Number}) =
            DH.Scalars.$name{T}
        @inline hartree_concretize_type(::Type{DD.Scalars.$name{TT, U} where TT},
                                        ::Type{T}) where U where T <: Number =
            DD.Scalars.$name{T, U}
        @inline hartree_concretize_type(::Type{DD.Scalars.$name{TT, U} where TT},
                                        ::Type{<: Quantity{T}}) where {U, T} =
            DD.Scalars.$name{T, U}
    end
end

"""
Tries and figures a concrete type given `Q <: Quantity` and `T <: Number`

`Q` may be abstract: the underlying type or the units may be unknown (but the dimensions are
always specified). Whatever is missing is taken from `T`. If the units are unspecified, and
`Q <: Dispatch.Dimensions.Scalars.All`, then Hartree units are used by default.

# Examples

concretize_type(Dispatch.Dimensions.Scalars.ρ, Int16) === Dispatch.Hartree.Scalars.ρ{Int16}
concretize_type(typeof(1.0u"m^-3"), Int16) === typeof(Int16(1)u"m^-3")
"""
@inline concretize_type(::Type{Quantity{T, D, U}}, ::Type{<:Number}) where {T, D, U} =
    Quantity{T, D, U}
@inline concretize_type(::Type{Quantity{T, D}},
                        ::Type{Quantity{TT, D, U} where TT}) where {T, D, U} =
    Quantity{T, D, U}
@inline concretize_type(::Type{Quantity{TT, D} where TT},
                        ::Type{Quantity{T, D, U}}) where {T, D, U} = Quantity{T, D, U}
@inline concretize_type(H::Type{<: Quantity{TT, D} where TT}, T::Type{<: Number}) where D =
    hartree_concretize_type(H, T)
@inline concretize_type(Q::Type{<: Quantity}, x::DenseArray) = concretize_type(Q, eltype(x))
@inline concretize_type(Q::Type{<: Quantity}, x::AxisArray) = concretize_type(Q, eltype(x))
@inline concretize_type(Q::Type{<: Quantity}, x::Number) = concretize_type(Q, typeof(x))

""" Gets standard name for a given physical dimension """
function standard_name end
for (_, name) in Dispatch.Dimensioned
    @eval standard_name(::Type{<: DD.Scalars.$name}) = $(QuoteNode(name))
end
standard_name(n::DD.Scalars.All) = standard_name(typeof(n))
standard_name(array::AbstractArray) = standard_name(eltype(array))
end
