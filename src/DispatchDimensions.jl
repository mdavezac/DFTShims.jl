module Dimensions
using ...UnitfulHartree
using ...Dispatch: Dimensioned, DFTAxisArray, DFTArray, Dimensioned
using AxisArrays, Unitful

const UH = UnitfulHartree

macro lintpragma(s) end
@lintpragma("Ignore use of undeclared variable T")
@lintpragma("Ignore use of undeclared variable U")
@lintpragma("Ignore use of undeclared variable N")
@lintpragma("Ignore use of undeclared variable A")
@lintpragma("Ignore use of undeclared variable Eₕ")

for (name, abbr) in Dimensioned
    DAA = Symbol("$(name)AxisArray")
    DDA = Symbol("$(name)DenseArray")
    DA = Symbol("$(name)Array")
    @eval begin
        const $name = Quantity{T, typeof(dimension(UH.$abbr)), U} where U where T
        const $DAA = DFTAxisArray{$name{T, U}, N, A} where U where A where N where T
        const $DDA = DenseArray{$name{T, U}, N} where U where N where T
        const $DA = DFTArray{$name{T, U}, N} where U where N where T
    end
end

""" Dimenional dispatch types for scalars """
module Scalars
    using ...Dispatch: Dimensioned, Dimensions
    for (name, abbr) in Dimensioned
        @eval begin
            const $name = Dimensions.$name
            const $abbr = $name
        end
    end
    const All = @eval Union{$((Expr(:curly, n, :T) for (n, a) in Dimensioned)...)} where T
    const ϵ = Eₕ
end

""" Dimensional dispatch types for Arrays """
module Arrays
    using ...Dispatch: Dimensioned, Dimensions
    for (name, abbr) in Dimensioned
        @eval begin
            const $name = Dimensions.$(Symbol("$(name)Array"))
            const $abbr = $name
        end
    end
    const All = @eval Union{
            $((Expr(:curly, n, :T, :N) for (n, a) in Dimensioned)...)} where {T, N}
    const ϵ = Eₕ
end

""" Dimenional dispatch types for AxisArrays """
module AxisArrays
    using ...Dispatch: Dimensioned, Dimensions
    for (name, abbr) in Dimensioned
        @eval begin
            const $name = Dimensions.$(Symbol("$(name)AxisArray"))
            const $abbr = $name
        end
    end
    const All = @eval Union{
            $((Expr(:curly, n, :T, :N, :A) for (n, a) in Dimensioned)...)} where {T, N, A}
    const ϵ = Eₕ
end

""" Dimenional dispatch types for DenseArrays """
module DenseArrays
    using ...Dispatch: Dimensioned, Dimensions
    for (name, abbr) in Dimensioned
        @eval begin
            const $name = Dimensions.$(Symbol("$(name)DenseArray"))
            const $abbr = $name
        end
    end
    const All = @eval Union{
                    $((Expr(:curly, n, :T, :N) for (n, a) in Dimensioned)...)} where {T, N}
    const ϵ = Eₕ
end


end
