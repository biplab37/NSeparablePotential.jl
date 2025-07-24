module NSeparablePotential

using LinearAlgebra
using NLsolve
using RecipesBase
using QuadGK
using UsefulFunctions

include("custom_types.jl")
include("helper_functions.jl")
include("separablerepresentation.jl")
include("gapeqn.jl")
include("LSequation.jl")
include("plot_recipes.jl")
include("one_separable.jl")
include("nseparable_phase_shift.jl")

end # module NSeparablePotential
