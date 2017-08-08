module ArrayInitialization

using Unitful
using AxisArrays
using ..ConstantArrays: ConstantArray
using ..UnitfulHartree
using ..Traits: components, ColinearSpin, SpinDegenerate, SpinCategory, ColinearSpin
using ..Traits: ColinearSpinFirst, ColinearSpinLast, is_spin_polarized
using ..Dispatch

macro lintpragma(s) end
@lintpragma("Ignore unused args")
@lintpragma("Ignore unused T")

const DD = Dispatch.Dimensions

_size(T::Type{<: DD.Scalars.All}, ::SpinDegenerate,
      ::SpinDegenerate, object::DD.AxisArrays.All) = size(object)
_size(T::Type{<: DD.Scalars.All}, ::ColinearSpin,
      ::SpinDegenerate, object::DD.AxisArrays.All) =
    (size(object)..., length(components(object)))

add_spin_axis(::SpinDegenerate, tup::Tuple, ::Any) = tup
add_spin_axis(::ColinearSpin, tup::Tuple, axis::Any) = tup..., axis
add_spin_axis(::ColinearSpinFirst, tup::Tuple, axis::Any) = axis, tup...

replace_spin_axis(T::Type{<: DD.Scalars.All}, ::SpinDegenerate,
                  sizes::Tuple, ax::Tuple) = begin
    @assert length(sizes.parameters) == length(ax.parameters)
    @assert findfirst(is_spin_polarized, ax.parameters) == 0
    sizes, ax
end

@generated replace_spin_axis(T::Type{<: DD.Scalars.All},
                             S::ColinearSpin, sizes::Tuple, ax::Tuple) = begin
    @assert length(sizes.parameters) == length(ax.parameters)
    index = findfirst(is_spin_polarized, ax.parameters)
    @assert index ≠ 0
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

for extension in [:zeros, :ones, :rand]
    @eval begin
        """
        Creates an array for the given DFT quantity
        
        The spin axis is automatically added if required. Neither `dims` nor `ax` should
        include the spin dimension.
        """
        Base.$extension(T::Type{<:DD.Scalars.All}, C::SpinCategory,
                        dims::Tuple, ax::Tuple) = begin
            length(ax) > (length(dims) + 1) && throw(ArgumentError("Too many axes"))
            comps = components(T, C)
            data = $extension(T, add_spin_axis(C, dims, length(comps)))
            # we can use a dummy array here since the underlying array (and indices) will be the
            # standard one
            defaults = AxisArrays.default_axes(ConstantArray(0, dims), ax)
            AxisArray(data, add_spin_axis(C, defaults, Axis{:spin}(comps)))
        end
        Base.$extension(T::Type{<:DD.Scalars.All}, C::SpinCategory, dims::Tuple) = begin
            comps = components(T, C)
            data = $extension(T, add_spin_axis(C, dims, length(comps)))
            defaults = AxisArrays.default_axes(ConstantArray(0, dims))
            AxisArray(data, add_spin_axis(C, defaults, Axis{:spin}(comps)))
        end

        """
        Creates an array for the given DFT quantity

        The spin components, if any, are added as the last dimension.
        Note should `dims` does not include the spin components.
        """
        Base.$extension(T::Type{<:DD.Scalars.All}, P::SpinCategory,
                   args::Vararg{<:Union{Integer, Axis}}) = begin
            @lintpragma("Ignore use of undeclared variable x")
            $extension(T, P, ((x for x in args if typeof(x) <: Integer)...),
                       ((x for x in args if typeof(x) <: Axis)...))
        end
        Base.$extension(T::Type{<:DD.Scalars.All}, polarized::Bool,
                        args::Vararg{<:Union{Integer, Axis}}) =
            $extension(T, polarized ? ColinearSpin(): SpinDegenerate(), args...)
    end
end
end
