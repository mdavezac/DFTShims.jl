module ArrayInitialization

using Unitful
using AxisArrays
using ..UnitfulHartree
using ..Traits: components, Polarized, Unpolarized, PolarizationCategory
using ..Dispatch

macro lintpragma(s) end

const DD = Dispatch.Dimensions

""" Creates an unpolarized array for the given DFT quantity """
zeros_impl_polarized(T::Type{<:DD.Scalars.All}, dims::Tuple, ax::Tuple) = begin
    length(ax) > (length(dims) + 1) && throw(ArgumentError("Too many axes"))
    comps = components(T, Polarized)
    data = zeros(T, (dims..., length(comps)))
    if length(ax) == length(dims) + 1
        AxisArray(data, ax...)
    else
        AxisArray(data, AxisArrays.default_axes(data, ax)[1:end - 1]..., Axis{:spin}(comps))
    end
end
""" Creates a polarized array for the given DFT quantity """
zeros_impl_unpolarized(T::Type{<:DD.Scalars.All}, dims::Tuple, ax::Tuple) =
    AxisArray(zeros(T, dims), ax...)
    
    
@lintpragma("Ignore unused args")
"""
Creates an array for the given DFT quantity

The spin components, if any, are added as the last dimension.
Note should `dims` does not include the spin components.
"""
Base.zeros(T::Type{<:DD.Scalars.All}, ::Type{Polarized},
           args::Vararg{<:Union{Integer, Axis}}) = begin
    @lintpragma("Ignore use of undeclared variable x")
    zeros_impl_polarized(T, ((x for x in args if typeof(x) <: Integer)...),
                         ((x for x in args if typeof(x) <: Axis)...))
end
Base.zeros(T::Type{<:DD.Scalars.All}, ::Type{Unpolarized},
           args::Vararg{<:Union{Integer, Axis}}) = begin
    @lintpragma("Ignore use of undeclared variable x")
    zeros_impl_unpolarized(T, ((x for x in args if typeof(x) <: Integer)...),
                           ((x for x in args if typeof(x) <: Axis)...))
end
Base.zeros(T::Type{<:DD.Scalars.All}, polarized::Bool,
           args::Vararg{<:Union{Integer, Axis}}) =
    zeros(T, polarized ? Polarized: Unpolarized, args...)
end
