using DFTShims: ConstantArrays
using AxisArrays

@test size(ConstantArray(1.0, (3, 4))) == (3, 4)
array = AxisArray(ConstantArray(1, (2, 2)), Axis{:spin}((:u, :d)))
@test size(array) == (2, 2)
@test axes(array, 1) == Axis{:spin}((:u, :d))
