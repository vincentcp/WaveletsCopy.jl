# disable precompilation during development
__precompile__(true)


module WaveletsCopy
using RecipesBase
using Reexport
using CardinalBSplines

include("sequences/sequence.jl")
include("filterbank.jl")
include("dwt/discretewavelets.jl")

@reexport using .DWT

end
