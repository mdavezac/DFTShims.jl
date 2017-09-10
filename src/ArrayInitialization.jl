module ArrayInitialization

using DocStringExtensions
using Unitful
using AxisArrays
using ArgCheck
using ..ConstantArrays: ConstantArray
using ..UnitfulHartree
using ..Traits: components, ColinearSpin, SpinDegenerate, SpinCategory, ColinearSpin,
                ColinearSpinFirst, ColinearSpinLast, is_spin_polarized, concretize_type,
                SpinAware
using ..Dispatch
export wrap

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

"""
Position for the spin axis for a given spin category

$(SIGNATURES)

`n` is the number of dimensions in the array, and `s` is the location of the current spin
axis, or 0 if there are none.
"""
spin_axis_position(::Type{ColinearSpin{T}}, ::Integer, ::Integer) where T = T
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
            @argcheck length(ax) <= (length(dims) + 1) "Too many axes"
            comps = components(T, C)
            data = $extension(T.parameters[1], add_spin_axis(C, dims, length(comps)))
            # we can use a dummy array here since the underlying array (and indices) will be
            # the standard one
            defaults = AxisArrays.default_axes(ConstantArray(0, dims), ax)
            AxisArray(reinterpret(T, data), add_spin_axis(C, defaults, Axis{:spin}(comps)))
        end
        Base.$extension(T::Type{<:DD.Scalars.All}, C::SpinCategory, dims::Tuple) = begin
            comps = components(T, C)
            data = $extension(T, add_spin_axis(C, dims, length(comps)))
            defaults = AxisArrays.default_axes(ConstantArray(0, dims))
            AxisArray(data, add_spin_axis(C, defaults, Axis{:spin}(comps)))
        end

        """
        Creates an array for the given DFT quantity

        The spin components, if any, are added as the dimension specified by the spin trait.
        Note that the spin trait is the only way polarization should enter the input
        arguments, e.g. the spin-axis should not be given explicitly, neither in the axes
        nor in the dimensions.
        """
        Base.$extension(T::Type{<:DD.Scalars.All}, P::SpinCategory,
                        args::Vararg{Union{Integer, Axis}}; kwargs...) = begin
            @lintpragma("Ignore use of undeclared variable x")
            dims = (
                (x for x in args if typeof(x) <: Integer)...,
                (length(v) for (_, v) in kwargs)...,
            )
            ax = (
                (x for x in args if typeof(x) <: Axis)...,
                (Axis{k}(v) for (k, v) in kwargs)...,
            )
            $extension(T, P, dims, ax)
        end
        Base.$extension(T::Type{<:DD.Scalars.All}, polarized::Bool,
                        args::Vararg{Union{Integer, Axis}}; kwargs...) =
            $extension(T, polarized ? ColinearSpin(): SpinDegenerate(), args...; kwargs...)
    end
end

for extension in [:zeros, :ones, :similar]
    private = Symbol(:_, extension)
    @eval begin
        $private(T::Type{<:DD.Scalars.All}, ::SpinDegenerate, ::SpinDegenerate,
            array::DD.AxisArrays.All) = $extension(array, concretize_type(T, array))

        $private(T::Type{<:DD.Scalars.All}, C::ColinearSpin,
                 ::SpinDegenerate, array::DD.AxisArrays.All) = begin
            comps = components(T, C)
            dims = add_spin_axis(C, size(array), length(comps))
            data = $extension(array.data, concretize_type(T, array), dims)
            defaults = AxisArrays.default_axes(ConstantArray(0, size(array)), axes(array))
            AxisArray(data, add_spin_axis(C, defaults, Axis{:spin}(comps)))
        end

        $private(T::Type{<:DD.Scalars.All}, C::ColinearSpin,
                 ::ColinearSpin, array::DD.AxisArrays.All) = begin
            dims = replace_spin_axis(T, C, size(array), axes(array))
            AxisArray($extension(array.data, concretize_type(T, array), dims[1]), dims[2])
        end

        $private(T::Type{<:DD.Scalars.All}, ::SpinDegenerate,
                    ::ColinearSpin, array::DD.AxisArrays.All) = begin
            dims = remove_spin_axis(size(array), axes(array))
            AxisArray($extension(array.data, concretize_type(T, array), dims[1]), dims[2])
        end

        Base.$extension(array::DD.AxisArrays.All, T::Type{<:DD.Scalars.All}, ::SpinAware) =
            $private(T, SpinCategory(array), SpinCategory(array), array)

        Base.$extension(array::DD.AxisArrays.All, T::Type{<:DD.Scalars.All},
                        wanted::SpinCategory) =
            $private(T, wanted, SpinCategory(array), array)
    end
end

Base.reinterpret(T::Type{<: DD.Scalars.All}, ::SpinDegenerate, array::DenseArray) =
    AxisArray(reinterpret(concretize_type(T, array), array), axes(array))
Base.reinterpret(T::Type{<: DD.Scalars.All}, ::ColinearSpinFirst, array::DenseArray) = 
    AxisArray(reinterpret(concretize_type(T, array), array),
              Axis{:spin}(components(T, ColinearSpinFirst())), Base.tail(axes(array))...)
Base.reinterpret(T::Type{<: DD.Scalars.All}, C::ColinearSpin, array::DenseArray) =
    AxisArray(reinterpret(concretize_type(T, array), array),
              Base.front(axes(array))..., Axis{:spin}(components(T, C)))

"""
    wrap([T::Type{<: Dispatch.Scalars.All}], [P::SpinCategory], array)

Wraps an existing array as a DFT array. The default spin-category is `SpinDegenerate`.
If the type and units are explicitly specified (first argument)), then the array must be
dimensionless. If they are not, then the array must have an element type `T` such that
`T <: Dispatch.Scalars.All`.
"""
wrap(T::Type{<: DD.Scalars.All}, P::SpinCategory, array::DenseArray) = begin
    @argcheck dimension(eltype(array)) == NoDims "Input array should be dimensionless"
    reinterpret(T, P, array)
end
wrap(::SpinDegenerate, array::DenseArray{<: DD.Scalars.All}) = AxisArray(array)
wrap(P::ColinearSpinFirst, array::DenseArray{<: DD.Scalars.All}) = begin
    c = components(eltype(array), P)
    @argcheck length(c) == size(array, 1) "Size of spin-axis does not match expectations"
    AxisArray(array, Axis{:spin}(c))
end
wrap(P::ColinearSpin, array::DenseArray{<: DD.Scalars.All}) = begin
    c = components(eltype(array), P)
    @argcheck(length(c) == size(array, ndims(array)),
              "Size of spin-axis does not match expectations")
    AxisArray(array, Base.front(axes(array))..., Axis{:spin}(c))
end
wrap(array::DenseArray{<: DD.Scalars.All}) = wrap(SpinDegenerate(), array)
wrap(T::Type{<: DD.Scalars.All}, array::DenseArray) =
    wrap(concretize_type(T, array), SpinDegenerate(), array)

""" Helper functions to convert between arrays with different units and memory layout """
_convert(T::Type{<:DD.Scalars.All}, C::ColinearSpin,
         C′::ColinearSpin, array::DD.AxisArrays.All) = begin
    T′ = concretize_type(T, array)
    (T′ == T && C == C′) ? array: copy!(similar(array, T, C), array)
end

"""
Converts axis to the requisite spin-axis location

$(SIGNATURES)

If the spin-axis does not move, then a reference to the original array is returned.
Otherwise, a new array is returned.
"""
Base.convert(::SpinAware, array::DD.AxisArrays.All) = array
Base.convert(T::Type{<: DD.Scalars.All}, ::SpinAware, array::DD.AxisArrays.All) =
    _convert(concretize_type(T, array), SpinCategory(array), SpinCategory(array), array)
Base.convert(C::ColinearSpin, array::DD.AxisArrays.All) =
    _convert(eltype(array), C, SpinCategory(array), array)
Base.convert(T::Type{<: DD.Scalars.All}, C::SpinCategory, array::DD.AxisArrays.All) =
    _convert(T, C, SpinCategory(array), array)
Base.convert(T::Type{<: DD.Scalars.All}, array::DD.AxisArrays.All) =
    _convert(T, SpinCategory(array), SpinCategory(array), array)

Unitful.uconvert(u::Unitful.Units, array::AxisArray) = begin
    data = similar(array.data, typeof(uconvert(u, oneunit(eltype(array)))))
    data .= array.data
    AxisArray(data, axes(array))
end

end
