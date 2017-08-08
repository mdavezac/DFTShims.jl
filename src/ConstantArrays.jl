module ConstantArrays
export ConstantArray

struct ConstantArray{T, N} <: AbstractArray{T, N}
    _value::T
    _size::NTuple{N, Int64}
end

(::Type{ConstantArray})(value::Any, dims::Vararg{<: Integer}) =
    ConstantArray{typeof(value), length(dims)}(value,
                                               convert(NTuple{length(dims), Int64}, dims))

Base.size(q::ConstantArray) = q._size
Base.getindex(q::ConstantArray, ::Integer) = q._value
Base.getindex(q::ConstantArray, ::Vararg{<: Integer}) = q._value

Base.show(io::IO, q::ConstantArray) = show(io, "ConstantArray($(q.value), $(q.dims))")
end
