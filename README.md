[![Build Status](https://travis-ci.org/mdavezac/DFTShims.jl.svg?branch=master)](https://travis-ci.org/mdavezac/DFTShims.jl)
[![Coverage Status](https://coveralls.io/repos/mdavezac/DFTShims.jl/badge.svg)](https://coveralls.io/r/mdavezac/DFTShims.jl)

# DFTShims

This modules provides an interface to manipulate
[Unitful](https://github.com/ajkeller34/Unitful.jl) arrays of interest in DFT.  The basic
underlying type is an [AxisArray](https://github.com/JuliaArrays/AxisArrays.jl) wrapping a
`DenseArray` of a kind or another, with physical dimensions corresponding to the quantitiy
of interest.

# Hartree units

The sub-package `DFTShims.UnitfulHartree` provides the following standard units for Hartree
atomic units, as well as some DFT specific quantities.

| abbreviation | Typename                | abbreviation | Typename                        |
|--------------|-------------------------|--------------|---------------------------------|
| mₑ           | ElectronMass            | ρ            | Density                         |
| e₀           | ElementaryCharge        | σₑ           | ContractedDensityGradient       |
| kₑ           | CoulombForceConstant    | ∂ϵ_∂ρ        | FirstDensityDerivative          |
| ħ            | ReducedPlanckConstant   | ∂ϵ_∂σ        | FirstGradientDerivative         |
| a₀           | BohrRadius              | ∂²ϵ_∂ρ²      | SecondDensityDerivative         |
| Eₕ           | HartreeEnergy           | ∂²ϵ_∂σ²      | SecondGradientDerivative        |
| Ry           | RydbergEnergy           | ∂²ϵ_∂ρ∂σ     | SecondDensityGradientDerivative |
| rₑ           | ClassicalElectronRadius | ∂³ϵ_∂ρ³      | ThirdDensityDerivative          |
|              |                         | ∂³ϵ_∂σ³      | ThirdGradientDerivative         |
|              |                         | ∂³ϵ_∂ρ²∂σ    | ThirdDensity2GradientDerivative |
|              |                         | ∂³ϵ_∂ρ∂σ²    | ThirdDensityGradient2Derivative |
|              |                         | ϵ(=== Eₕ)    | HartreeEnergy                   |

The following constants are also declared:

~~~Julia
const α = 1e₀^2*1kₑ/(1Unitful.c*ħ)
const mₚ = 1836.15mₑ
const μ_b = e₀*ħ/(2mₑ)
const ϵ₀ = 1/(4π*kₑ)
~~~

# Axis-arrays

Axis arrays provide a conveniently flexible data type which can be used to describe most any
data with homogeneous units. The main advantage is the ability to describe each axis
explicitly, e.g. whether the axis relates to spin, positions in Cartesian coordinates,
wavefunction number etc... DFTShims makes it easy to define functions taking axis-arrays
with specific physical dimensions (say ρ, or ∂ϵ/∂ρ). It also provides traits to detect
whether an axis-array is spin-polarized. Finally, it overloads a number of functions to make
it easier to create and manipulate such arrays.

## Dispatch over physical dimensions

The objective is to easily specify functions that take axis-arrays with physical dimensions
and physical units. In Unitful, dimension relates to the physical meaning of the quantity,
e.g. length, whereas units implies both a physical dimension and a specific measurement,
e.g. _meters_. `DFTShims` allows dispatch over both separately, for both scalars and axis
arrays. The dispatch types are separated into modules:

- `DFTShims.Dispatch`
  * `Dimensions`
    + `Scalars`
    + `AxisArrays`
  * `Hartree`
    + `Scalars`
    + `AxisArrays`

All `Scalars` and `AxisArrays` contain parametric types for `ρ`, `σ`, `ϵ`, and the derivate
of the latter versus the two former up to degree three. The parametric types are first
parameterized of the underlying scalar (`Float64`, 'Int16', etc). `Dimensions` are
also parameterized over the actual units. And `AxisArrays` are further parameterized over
the exact underlying array and axis.

Hence we can create the following:

~~~Julia
using DFTShims: Dispatch
f(x::Dispatch.Dimensions.Scalars.ρ) =  "any ρ scalar"
f(x::Dispatch.Dimensions.Scalars.ρ{<: Integer}) = "integers with dimension ρ"
f(x::Dispatch.Hartree.Scalars.ρ) =  "any ρ scalar with coorrect hartree units"
f(x::Dispatch.Hartree.Scalars.ρ{<: Integer}) = "integers with hartree units ρ"
~~~

The reader is invited to play with `f(1.0u"ρ")`, `f(1u"ρ")` (both in Hartree),
`f(1u"m^-3")`, and so on. After reading the section below on easily creating `AxisArrays`,
the reader may want to define methods such ash `f(x::Dispatch.Hartree.AxisArrays.ρ)`.

To be more specific, each of `ρ` and friends in `Scalars` is an alias for
`Quantity{T, D, U}` where `D` - the physical dimension - is specified, `U` - the units - is
left unspecified in `Dispatch` but pertains to the atomic units in `Hartree`. In
`AxisArrays`, similar aliases are defined for
`AxisArrays{Quantity{T, D, U}, N, <: DenseArray{Quantity{T, D, U}, N}, AXES}`. In
other words the aliases in `Scalars` allow multiple dispatch over scalars, whereas the
aliases in `AxisArrays` allow multiple dispatch of axis-arrays of the scalars.

For convenience, `DFTShims.Dispatch` provides
`const Scalars === Dispatch.Dimensions.Scalars`.

## Spin

The main advantage of using axis-arrays is that it allows us to define whether a quantity is
spin-polarized from its type, as well as figure out how the polarization is managed in
memory.

An axis-array is polarized if it sports an axis with then name `:spin` and containing more
than one component. For instance:

~~~Julia
using AxisArrays
using DFTShims
@assert is_spin_polarized(AxisArrays(zeros(2, 3), Axis{:spin}((:α, :β))))
~~~

By default, the components for a given quantity are obtained from the function `components`.

In general, the spin-axis can be the fastest changing (first) axis, or the slowest (last),
or anything in between. A specialized trait hierarchy deriving from `SpinCategory` is
available to specify the preferred option:

- `struct SpinDenegenerate <: SpinCategory end`
- `abstract type ColinearSpin <: SpinCategory end`: all spin-polarized options in the
colinear spin approximation
- `abstract type ColinearSpinFirst <: ColinearSpin end`: spin-axis is always the first axis
(fastest changing)
- `abstract type ColinearSpinLast <: ColinearSpin end`: spin-axis is always the last axis
(slowest changing)
- `abstract type ColinearSpinPreferLast <: ColinearSpin end`: spin-axis is set to last
unless it is already set. This is convenient for functions that can handle any spin-axis
location and would want to avoid copying data where possible. However, it still specifies
a preferred axis location in cases where it is not known from the input. In practice, this
trait is useful only when creating a new array from another.

The function `SpinCategory` can be used on an array to guess the relevant trait, whether
`SpinDegenerate` or some flavor of `ColinearSpin`. 

## Array creation

The array creation functions `zeros`, `ones`, `rand`, and `similar`, have been overloaded to
easily create spin-polarized and spin-degenerate arrays with the correct units. The location
of the spin-axis depends on the input trait:

```Julia
a = zeros(Dispatch.Hartree.Scalars.ρ{Int64}, SpinDegenerate(), (8, 9))
@test size(a) == (8, 9) && unit(a) === UnitfulHartree.ρ

b = similar(Dispatch.Hartree.Scalars.∂³ϵ_∂ρ³{Int64}, ColinearSpinLast(), (8, 9))
@test size(b) == (8, 9, 4)

c = ones(Dispatch.Hartree.Scalars.∂³ϵ_∂ρ³{Int64}, ColinearSpinFirst(), b)
@test size(c) == (4, 8, 9)

d = ones(Dispatch.Dimensions.Scalars.ρ, c)
@test size(d) == (4, 8, 9)
```

Note in the last example, the first argument is an abstract type which specifies only the
physical dimension (but not the underlying type, nor the units). The underlying type will
guessed from the type of the second argument. And the units, will be either taken from `c`
or default to Hartrees (Here, `dimension(1u"ρ") ≠ dimension(c)`, hence the units are
defaulted to Hartree).

When specifying the dimensions of the array, the spin axis should be omitted, it is fully
specified by the spin-trait and the physical units. This ensures that the name of the
components of the spin-axis are correct (try `axes(c, 1)`).

It is also possible to convert between different spin-axis locations with:

```Julia
convert(ColinearSpinLast(), c)
convert(Dispatch.Hartree.ρ{Float64}, ColinearSpinFirst(), c)
```
