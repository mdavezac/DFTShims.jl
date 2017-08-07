module ArrayInitialization

using Unitful
using AxisArrays
using ..UnitfulHartree
using ..Traits: components, ColinearSpin, SpinDegenerate, SpinCategory, ColinearSpin
using ..Traits: ColinearSpinFirst, ColinearSpinLast, is_spin_polarized
using ..Dispatch

macro lintpragma(s) end
@lintpragma("Ignore unused args")
@lintpragma("Ignore unused T")

const DD = Dispatch.Dimensions

""" Creates the tuple of axes for a spin-polarized quantity """
@generated polarized_axis(comps::Tuple, defaults::Tuple, ax::Tuple) = begin
    if length(ax.parameters) == length(defaults.parameters)
        :(defaults)
    elseif length(ax.parameters) == length(defaults.parameters) - 1
        :(ax..., Axis{:spin}(comps))
    else
        :(Base.front(defaults)..., Axis{:spin}(comps))
    end
end

""" Creates an unpolarized array for the given DFT quantity """
Base.zeros(T::Type{<:DD.Scalars.All}, ::ColinearSpin,  dims::Tuple, ax::Tuple) = begin
    length(ax) > (length(dims) + 1) && throw(ArgumentError("Too many axes"))
    comps = components(T, ColinearSpin())
    data = zeros(T, (dims..., length(comps)))
    defaults = AxisArrays.default_axes(data, ax)
    AxisArray(data, polarized_axis(comps, defaults, ax))
end
Base.zeros(T::Type{<:DD.Scalars.All}, ::ColinearSpin, dims::Tuple) = begin
    comps = components(T, ColinearSpin())
    data = zeros(T, (dims..., length(comps)))
    defaults = AxisArrays.default_axes(data)
    AxisArray(data, Base.front(defaults)..., Axis{:spin}(comps))
end
""" Creates an unpolarized array for the given DFT quantity """
Base.zeros(T::Type{<:DD.Scalars.All}, ::SpinDegenerate, dims::Tuple, ax::Tuple) =
    AxisArray(zeros(T, dims), ax...)
Base.zeros(T::Type{<:DD.Scalars.All}, ::SpinDegenerate, dims::Tuple) =
    AxisArray(zeros(T, dims))
    
    
"""
Creates an array for the given DFT quantity

The spin components, if any, are added as the last dimension.
Note should `dims` does not include the spin components.
"""
Base.zeros(T::Type{<:DD.Scalars.All}, P::SpinCategory,
           args::Vararg{<:Union{Integer, Axis}}) = begin
    @lintpragma("Ignore use of undeclared variable x")
    zeros(T, P, ((x for x in args if typeof(x) <: Integer)...),
            ((x for x in args if typeof(x) <: Axis)...))
end
Base.zeros(T::Type{<:DD.Scalars.All}, polarized::Bool,
           args::Vararg{<:Union{Integer, Axis}}) =
    zeros(T, polarized ? ColinearSpin(): SpinDegenerate(), args...)

_size(T::Type{<: DD.Scalars.All}, ::SpinDegenerate,
      ::SpinDegenerate, object::DD.AxisArrays.All) = size(object)
_size(T::Type{<: DD.Scalars.All}, ::ColinearSpin,
      ::SpinDegenerate, object::DD.AxisArrays.All) =
    (size(object)..., length(components(object)))

add_spin_axis(::ColinearSpin, tup::Tuple, axis::Any) = tup..., axis
add_spin_axis(::ColinearSpinFirst, tup::Tuple, axis::Any) = axis, tup...

@generated replace_spin_axis(T::Type{<: DD.Scalars.All},
                             S::ColinearSpin, sizes::Tuple, ax::Tuple) = begin
    @assert(length(sizes.parameters) == length(ax.parameters))
    index = findfirst(is_spin_polarized, ax.parameters)
    @assert(index ≠ 0)
    left = :(sizes, ax)
    for i in index:length(ax.parameters)
        left = :(Base.front($(left.args[1])), Base.front($(left.args[2])))
    end
    right = :(Base.reverse(sizes), Base.reverse(ax))
    for i in 1:index
        right = :(Base.front($(right.args[1])), Base.front($(right.args[2])))
    end
    right = :(Base.reverse($(right.args[1])), Base.reverse($(right.args[2])))
    quote
        comps = components(T, S)
        (
            ($(left.args[1])..., length(comps), $(right.args[1])...),
            ($(left.args[2])..., Axis{:spin}(comps), $(right.args[2])...)
        )
    end
end

@generated replace_spin_axis(T::Type{<: DD.Scalars.All},
                             S::ColinearSpinFirst, sizes::Tuple, ax::Tuple) = begin
    @assert(length(sizes.parameters) == length(ax.parameters))
    index = findfirst(is_spin_polarized, ax.parameters)
    @assert(index ≠ 0)
    left = :(sizes, ax)
    for i in index:length(ax.parameters)
        left = :(Base.front($(left.args[1])), Base.front($(left.args[2])))
    end
    right = :(Base.reverse(sizes), Base.reverse(ax))
    for i in 1:index
        right = :(Base.front($(right.args[1])), Base.front($(right.args[2])))
    end
    right = :(Base.reverse($(right.args[1])), Base.reverse($(right.args[2])))
    quote
        comps = components(T, S)
        (
            (length(comps), $(left.args[1])..., $(right.args[1])...),
            (Axis{:spin}(comps), $(left.args[2])..., $(right.args[2])...)
        )
    end
end

@generated replace_spin_axis(T::Type{<: DD.Scalars.All},
                             S::ColinearSpinLast, sizes::Tuple, ax::Tuple) = begin
    @assert(length(sizes.parameters) == length(ax.parameters))
    index = findfirst(is_spin_polarized, ax.parameters)
    @assert(index ≠ 0)
    left = :(sizes, ax)
    for i in index:length(ax.parameters)
        left = :(Base.front($(left.args[1])), Base.front($(left.args[2])))
    end
    right = :(Base.reverse(sizes), Base.reverse(ax))
    for i in 1:index
        right = :(Base.front($(right.args[1])), Base.front($(right.args[2])))
    end
    right = :(Base.reverse($(right.args[1])), Base.reverse($(right.args[2])))
    quote
        comps = components(T, S)
        (
            ($(left.args[1])..., $(right.args[1])..., length(comps)),
            ($(left.args[2])..., $(right.args[2])..., Axis{:spin}(comps))
        )
    end
end
# _similar(T::Type{<: DD.Scalars.All},
#          wanted::ColinearSpin, has::ColinearSpin,
#          object::DD.AxisArrays.All) = begin
#     AxisArray(similar(object.data, T, _size(T, has, object)), _axes(T, has, object))
# end
# Base.similar(T::Type{<: DD.Scalars.All}, object::DD.AxisArrays.All) =
#     _similar(T, SpinCategory(object), SpinCategory(object), object)
# Base.similar(T::Type{<: DD.Scalars.All},
#              wanted_pol::SpinCategory, object::DD.AxisArrays.All) =
#     _similar(T, wanted_po, SpinCategory(object), object)
end
