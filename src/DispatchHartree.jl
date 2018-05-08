module Hartree
using ...UnitfulHartree
using ...Dispatch: Dimensioned, DFTAxisArray, DFTArray, Dimensions
using AxisArrays, Unitful

const UH = UnitfulHartree

macro lintpragma(s) end
@lintpragma("Ignore use of undeclared variable T")
@lintpragma("Ignore use of undeclared variable N")
@lintpragma("Ignore use of undeclared variable Eₕ")

for (name, abbr) in Dimensioned
    DAA = Symbol("$(name)AxisArray")
    DDA = Symbol("$(name)DenseArray")
    DA = Symbol("$(name)Array")
    if abbr == :σ
       abbr = :σₑ
    end
    @eval begin
        const $name = Dimensions.$name{T, typeof(UH.$abbr)} where T
        const $DAA = Dimensions.$DAA{T, N, typeof(UH.$abbr)} where N where T
        const $DDA = Dimensions.$DDA{T, N, typeof(UH.$abbr)} where N where T
        const $DA = Dimensions.$DA{T, N, typeof(UH.$abbr)} where N where T
    end
end

""" Dimensional dispatch types for scalars """
module Scalars
    using ...Dispatch: Dimensioned, Hartree
    for (name, abbr) in Dimensioned
        @eval begin
            const $name = Hartree.$name
            const $abbr = $name
        end
    end
    const All = @eval Union{$((Expr(:curly, n, :T) for (n, a) in Dimensioned)...)} where T
    const ϵ = Eₕ
end

""" Dimensional dispatch types for Arrays """
module Arrays
    using ...Dispatch: Dimensioned, Hartree
    for (name, abbr) in Dimensioned
        @eval begin
            const $name = Hartree.$(Symbol("$(name)Array"))
            const $abbr = $name
        end
    end
    const All = @eval Union{
            $((Expr(:curly, n, :T, :N) for (n, a) in Dimensioned)...)} where {T, N}
    const ϵ = Eₕ
end

""" Dimenional dispatch types for AxisArrays """
module AxisArrays
    using ...Dispatch: Dimensioned, Hartree
    for (name, abbr) in Dimensioned
        @eval begin
            const $name = Hartree.$(Symbol("$(name)AxisArray"))
            const $abbr = $name
        end
    end
    const All = @eval Union{
            $((Expr(:curly, n, :T, :N, :A) for (n, a) in Dimensioned)...)} where {T, N, A}
    const ϵ = Eₕ
end

""" Dimenional dispatch types for DenseArrays """
module DenseArrays
    using ...Dispatch: Dimensioned, Hartree
    for (name, abbr) in Dimensioned
        @eval begin
            const $name = Hartree.$(Symbol("$(name)DenseArray"))
            const $abbr = $name
        end
    end
    const All = @eval Union{
                    $((Expr(:curly, n, :T, :N) for (n, a) in Dimensioned)...)} where {T, N}
    const ϵ = Eₕ
end

end
