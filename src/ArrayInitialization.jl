module ArrayInitialization

using DocStringExtensions
using Unitful
using AxisArrays
using ..ConstantArrays: ConstantArray
using ..UnitfulHartree
using ..Traits: components, ColinearSpin, SpinDegenerate, SpinCategory, ColinearSpin,
                ColinearSpinFirst, ColinearSpinLast, ColinearSpinPreferLast,
                is_spin_polarized
using ..Dispatch
export add_spin_axis

macro lintpragma(s) end
@lintpragma("Ignore unused args")
@lintpragma("Ignore unused T")

const DD = Dispatch.Dimensions

_size(::Type{<: DD.Scalars.ϵ}, ::SpinCategory, ::SpinCategory, object::DD.AxisArrays.ϵ) =
    size(object)
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

"""
Position for the spin axis for a given spin category

$(SIGNATURES)

`n` is the number of dimensions in the array, and `s` is the location of the current spin
axis, or 0 if there are none.
"""
spin_axis_position(::Type{ColinearSpinPreferLast}, n::Integer, s::Integer) =
    !(1 ≤ s ≤ n) ? n + 1: s
spin_axis_position(::Type{ColinearSpinFirst}, ::Integer, ::Integer) = 1
spin_axis_position(::Type{ColinearSpinLast}, n::Integer, s::Integer) =
    !(1 ≤ s ≤ n) ? n + 1: n
spin_axis_position(S::SpinCategory, n::Integer, s::Integer) =
    spin_axis_position(typeof(S), n, s)
spin_axis_position(array::DD.AxisArrays.All) = findfirst(is_spin_polarized, axes(array))
spin_axis_position(T::Type{DD.AxisArrays.All}) =
    findfirst(is_spin_polarized, T.parameters[end].parameters)

function _move_and_re(name::Symbol, expr::Expr, orig::Integer, final::Integer, n::Integer)
    result = Expr(:tuple)
    for i in 1:n
        i > orig && push!(result.args, :($name[$i]))
        i == final && push!(result.args, expr)
        i < orig && push!(result.args, :($name[$i]))
    end
    final == n + 1 && push!(result.args, expr)
    result
end

@generated replace_spin_axis(T::Type{<: DD.Scalars.All},
                             S::ColinearSpin, sizes::Tuple, ax::Tuple) = begin
    @lintpragma("Ignore use of undeclared variable comps")
    @assert length(sizes.parameters) == length(ax.parameters)
    current = findfirst(is_spin_polarized, ax.parameters)
    final = spin_axis_position(S, length(sizes.parameters), current)
    @assert final ≠ 0
    d = _move_and_re(:sizes, :(length(comps)), current, final, length(sizes.parameters))
    a = _move_and_re(:ax, :(Axis{:spin}(comps)), current, final, length(sizes.parameters))
    quote
        comps = components(T, S)
        ($d, $a)
    end
end

@generated remove_spin_axis(sizes::Tuple, ax::Tuple) = begin
    @lintpragma("Ignore use of undeclared variable comps")
    @assert length(sizes.parameters) == length(ax.parameters)
    current = findfirst(is_spin_polarized, ax.parameters)
    d = _move_and_re(:sizes, :(length(comps)), current, 0, length(sizes.parameters))
    a = _move_and_re(:ax, :(Axis{:spin}(comps)), current, 0, length(sizes.parameters))
    :(($d, $a))
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
        Base.$extension(T::Type{<:DD.Scalars.ϵ}, ::ColinearSpin, dims::Tuple, ax::Tuple) =
            $extension(T, SpinDegenerate(), dims, ax)
        Base.$extension(T::Type{<:DD.Scalars.ϵ}, ::ColinearSpin, dims::Tuple) =
            $extension(T, SpinDegenerate(), dims)

        """
        Creates an array for the given DFT quantity

        The spin components, if any, are added as the last dimension.
        Note should `dims` does not include the spin components.
        """
        Base.$extension(T::Type{<:DD.Scalars.All}, P::SpinCategory,
                        args::Vararg{Union{Integer, Axis}}) = begin
            @lintpragma("Ignore use of undeclared variable x")
            $extension(T, P, ((x for x in args if typeof(x) <: Integer)...),
                       ((x for x in args if typeof(x) <: Axis)...))
        end
        Base.$extension(T::Type{<:DD.Scalars.All}, polarized::Bool,
                        args::Vararg{Union{Integer, Axis}}) =
            $extension(T, polarized ? ColinearSpin(): SpinDegenerate(), args...)
    end
end

for extension in [:zeros, :ones, :similar]
    private = Symbol(:_, extension)
    @eval begin
        $private(T::Type{<:DD.Scalars.All}, ::SpinDegenerate, ::SpinDegenerate,
               array::DD.AxisArrays.All) = $extension(array, T)

        $private(T::Type{<:DD.Scalars.All}, C::ColinearSpin,
                    ::SpinDegenerate, array::DD.AxisArrays.All) = begin
            comps = components(T, C)
            dims = add_spin_axis(C, size(array), length(comps))
            data = $extension(array.data, T, dims)
            defaults = AxisArrays.default_axes(ConstantArray(0, size(array)), axes(array))
            AxisArray(data, add_spin_axis(C, defaults, Axis{:spin}(comps)))
        end

        $private(T::Type{<:DD.Scalars.All}, C::ColinearSpin,
                    ::ColinearSpin, array::DD.AxisArrays.All) = begin
            dims = replace_spin_axis(T, C, size(array), axes(array))
            AxisArray($extension(array.data, T, dims[1]), dims[2])
        end

        $private(T::Type{<:DD.Scalars.All}, ::SpinDegenerate,
                    ::ColinearSpin, array::DD.AxisArrays.All) = begin
            dims = remove_spin_axis(size(array), axes(array))
            AxisArray($extension(array.data, T, dims[1]), dims[2])
        end

        Base.$extension(T::Type{<:DD.Scalars.All}, array::DD.AxisArrays.All) =
            $private(T, SpinCategory(array), SpinCategory(array), array)
        Base.$extension(T::Type{<:DD.Scalars.ϵ}, array::DD.AxisArrays.All) =
            $private(T, SpinDegenerate(), SpinCategory(array), array)

        Base.$extension(T::Type{<:DD.Scalars.All},
                        wanted::SpinCategory, array::DD.AxisArrays.All) =
            $private(T, wanted, SpinCategory(array), array)
        Base.$extension(T::Type{<:DD.Scalars.ϵ}, ::SpinCategory, array::DD.AxisArrays.All) =
            $private(T, SpinDegenerate(), SpinCategory(array), array)
    end
end

Base.reinterpret(T::Type{<: DD.Scalars.All}, ::SpinDegenerate, array::DenseArray) =
    AxisArray(reinterpret(T, array), axes(array))
Base.reinterpret(T::Type{<: DD.Scalars.All}, ::ColinearSpinFirst, array::DenseArray) = 
    AxisArray(reinterpret(T, array),
              Axis{:spin}(components(T, ColinearSpinFirst())), Base.tail(axes(array))...)
Base.reinterpret(T::Type{<: DD.Scalars.All}, C::ColinearSpin, array::DenseArray) =
    AxisArray(reinterpret(T, array),
              Base.front(axes(array))..., Axis{:spin}(components(T, C)))
Base.reinterpret(T::Type{<: DD.Scalars.ϵ}, ::ColinearSpin, array::DenseArray) =
    reinterpret(T, SpinDegenerate(), array)
Base.reinterpret(T::Type{<: DD.Scalars.ϵ}, ::ColinearSpinFirst, array::DenseArray) =
    reinterpret(T, SpinDegenerate(), array)

"""
Converts axis to the requisite spin-axis location

$(SIGNATURES)

If the spin-axis does not move, then a reference to the original array is returned.
Otherwise, a new array is returned.
"""
Base.convert(C::Type{<: ColinearSpin}, array::DD.AxisArrays.All) = begin
    original = spin_axis_position(array)
    final = spin_axis_position(C, ndims(array), original)
    final == original && return array

    iₛ = insert!(deleteat!(collect(1:ndims(array)), original), final, original)
    permutedims(array, iₛ)
end
Base.convert(::Type{ColinearSpinPreferLast}, array::DD.AxisArrays.All) = array
Base.convert(::Type{<: ColinearSpin}, array::DD.AxisArrays.ϵ) = array
Base.convert(S::SpinCategory, array::DD.AxisArrays.All) = convert(typeof(S), array)

Unitful.uconvert(u::Unitful.Units, array::AxisArray) = begin
    data = similar(array.data, typeof(uconvert(u, oneunit(eltype(array)))))
    data .= array.data
    AxisArray(data, axes(array))
end

end
