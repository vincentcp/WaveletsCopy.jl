# disable precompilation during development
#__precompile__()


module Wavelets
using RecipesBase
using Reexport
using CardinalBSplines

# Load module Sequences
include("sequences/sequence.jl")
# Load module Filterbank
include("filterbanks/filterbank.jl")

include("dwt/discretewavelets.jl")

@reexport using .DWT

end
