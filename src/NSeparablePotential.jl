module NSeparablePotential

using LinearAlgebra
using NLsolve
using RecipesBase
using QuadGK

include("custom_types.jl")
include("helper_funcitons.jl")
include("separablerepresentation.jl")
include("gapeqn.jl")
include("LSequation.jl")
include("plot_recipes.jl")

end # module NSeparablePotential
