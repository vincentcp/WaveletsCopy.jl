# disable precompilation during development
#__precompile__()


module Wavelets
using RecipesBase
using Reexport
# Code from early Wavelet Toolbox
include("util.jl")
@reexport using .Util

include("sequences/sequence.jl")
include("filterbank.jl")
include("dwt/discretewavelets.jl")
include("util/recipes.jl")

@reexport using .DWT, .Sequences, .Filterbanks

end
