module UnitfulHartree
using Unitful
import Unitful

macro lintpragma(s) end
@lintpragma("Ignore use of undeclared variable ElectronMass")
@lintpragma("Ignore use of undeclared variable mₑ")
@lintpragma("Ignore use of undeclared variable ElementaryCharge")
@lintpragma("Ignore use of undeclared variable e₀")
@lintpragma("Ignore use of undeclared variable CoulombForceConstant")
@lintpragma("Ignore use of undeclared variable kₑ")
@lintpragma("Ignore use of undeclared variable ReducedPlanckConstant")
@lintpragma("Ignore use of undeclared variable ħ")
@lintpragma("Ignore use of undeclared variable BohrRadius")
@lintpragma("Ignore use of undeclared variable a₀")
@lintpragma("Ignore use of undeclared variable HartreeEnergy")
@lintpragma("Ignore use of undeclared variable Eₕ")
@lintpragma("Ignore use of undeclared variable RydbergEnergy")
@lintpragma("Ignore use of undeclared variable Ry")
@lintpragma("Ignore use of undeclared variable ClassicalElectronRadius")
@lintpragma("Ignore use of undeclared variable rₑ")
@lintpragma("Ignore use of undeclared variable Density")
@lintpragma("Ignore use of undeclared variable ρ")
@lintpragma("Ignore use of undeclared variable ContractedDensityGradient")
@lintpragma("Ignore use of undeclared variable σ")
@lintpragma("Ignore use of undeclared variable FirstDensityDerivative")
@lintpragma("Ignore use of undeclared variable ∂ϵ_∂ρ")
@lintpragma("Ignore use of undeclared variable FirstGradientDerivative")
@lintpragma("Ignore use of undeclared variable ∂ϵ_∂σ")
@lintpragma("Ignore use of undeclared variable SecondDensityDerivative")
@lintpragma("Ignore use of undeclared variable ∂²ϵ_∂ρ²")
@lintpragma("Ignore use of undeclared variable SecondGradientDerivative")
@lintpragma("Ignore use of undeclared variable ∂²ϵ_∂σ²")
@lintpragma("Ignore use of undeclared variable SecondDensityGradientDerivative")
@lintpragma("Ignore use of undeclared variable ∂²ϵ_∂ρ∂σ")
@lintpragma("Ignore use of undeclared variable ThirdDensityDerivative")
@lintpragma("Ignore use of undeclared variable ∂³ϵ_∂ρ³")
@lintpragma("Ignore use of undeclared variable ThirdGradientDerivative")
@lintpragma("Ignore use of undeclared variable ∂³ϵ_∂σ³")
@lintpragma("Ignore use of undeclared variable ThirdDensity2GradientDerivative")
@lintpragma("Ignore use of undeclared variable ∂³ϵ_∂ρ²∂σ")
@lintpragma("Ignore use of undeclared variable ThirdDensityGradient2Derivative")
@lintpragma("Ignore use of undeclared variable ∂³ϵ_∂ρ∂σ²")

@unit mₑ "mₑ" ElectronMass 9.1093835611e-31*u"kg" false
@unit e₀ "eₒ" ElementaryCharge 1.602176620898e-19*u"C" false
@unit kₑ "kₑ" CoulombForceConstant 1/(4π)u"ϵ0" false
@unit ħ "ħ" ReducedPlanckConstant Unitful.ħ false
@unit a₀ "a₀" BohrRadius 1ħ^2/(1kₑ*mₑ*e₀^2) false
@unit Eₕ "Eₕ" HartreeEnergy mₑ*e₀^4*kₑ^2/(1ħ^2) true
@unit Ry "Ry" RydbergEnergy 0.5Eₕ true
@unit rₑ "rₑ" ClassicalElectronRadius 1e₀^2*kₑ/(1mₑ*Unitful.c^2) false
const α = 1e₀^2*1kₑ/(1Unitful.c*ħ)
const mₚ = 1836.15mₑ
const μ_b = e₀*ħ/(2mₑ)
const ϵ₀ = 1/(4π*kₑ)

@unit ρ         "ρ"         Density                         (1a₀)^-3       false
@unit σ         "σ"         ContractedDensityGradient       (1a₀)^-8       false
@unit ∂ϵ_∂ρ     "∂ϵ_∂ρ"     FirstDensityDerivative          (1Eₕ)*(1a₀)^3  false
@unit ∂ϵ_∂σ     "∂ϵ_∂σ"     FirstGradientDerivative         (1Eₕ)*(1a₀)^8  false
@unit ∂²ϵ_∂ρ²   "∂²ϵ_∂ρ²"   SecondDensityDerivative         (1Eₕ)*(1a₀)^6  false
@unit ∂²ϵ_∂σ²   "∂²ϵ_∂σ²"   SecondGradientDerivative        (1Eₕ)*(1a₀)^16 false
@unit ∂²ϵ_∂ρ∂σ  "∂²ϵ_∂ρ∂σ"  SecondDensityGradientDerivative (1Eₕ)*(1a₀)^11 false
@unit ∂³ϵ_∂ρ³   "∂³ϵ_∂ρ³"   ThirdDensityDerivative          (1Eₕ)*(1a₀)^9  false
@unit ∂³ϵ_∂σ³   "∂³ϵ_∂σ³"   ThirdGradientDerivative         (1Eₕ)*(1a₀)^24 false
@unit ∂³ϵ_∂ρ²∂σ "∂³ϵ_∂ρ²∂σ" ThirdDensity2GradientDerivative (1Eₕ)*(1a₀)^14 false
@unit ∂³ϵ_∂ρ∂σ² "∂³ϵ_∂ρ∂σ²" ThirdDensityGradient2Derivative (1Eₕ)*(1a₀)^19 false
const ϵ = Eₕ


# Some gymnastics required here because if we precompile, we cannot add to
# Unitful.basefactors at compile time and expect the changes to persist to runtime.
const localunits = Unitful.basefactors
function __init__()
    merge!(Unitful.basefactors, localunits)
    Unitful.register(UnitfulHartree)
end


end # module
