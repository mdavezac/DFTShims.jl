Base.copy!(dest::AxisArray{TT, N, AA, As},
           source::AxisArray{T, N, A, As}) where {TT, N, AA, T, A, As} = begin
    copy!(dest.data, source.data)
    dest
end

"""
    copy!(dest::AxisArray, source::AxisArray, ignore_axes::Bool=false) -> AxisArray

Copies data from one AxisArray to another. If `ignore_axes` is `false` (default), then all
axis names should match, though they may be permutated. If `ignore_axes` is `true`, then the
function copies the underlying arrays, using the behaviour standard for the types of the
underlying arrays.
"""
Base.copy!(dest::AxisArray{TT, N, AA, AAs}, source::AxisArray{T, N, A, As},
           ignore_axes::Bool=false) where {TT, N, AA, AAs, T, A, As} = begin
    if ignore_axes
        copy!(dest.data, source.data)
    else
        perm = indexin(collect(axisnames(dest)), collect(axisnames(source)))
        any(x==nothing for x in perm) &&
                throw(ArgumentError("""
                    Axes of source and destination do not match.
                    Use copy!(A, B, false) if you want to ignore axis information.
                    """))

        #REVIEW: not sure if this is correct?
        permutedims!(dest.data, convert.(TT, source.data), perm)
    end
    dest
end

"""
    copy!(dest::AxisArray, daxes::Tuple, src::AxisArray, saxes::Tuple, ignore_axes=false)

Copies a block from one array to a block of another. The tuples identifying the blocks can
be anything that `view` accepts. Hence it may include axis information.
"""
Base.copy!(dest::AxisArray, daxes::Tuple{Axis, Vararg{Axis}},
           src::AxisArray, saxes::Tuple{Axis, Vararg{Axis}}) =
    copy!(view(dest, daxes...), view(src, saxes...))
