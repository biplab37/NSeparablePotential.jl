module NSeparablePotential

using LinearAlgebra
using NLsolve
using RecipesBase

include("custom_types.jl")
include("separablerepresentation.jl")
include("gapeqn.jl")
include("LDequation.jl")
include("plot_recipes.jl")

end # module NSeparablePotential