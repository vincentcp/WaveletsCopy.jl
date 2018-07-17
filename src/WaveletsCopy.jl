# disable precompilation during development
__precompile__(true)


module WaveletsCopy

using RecipesBase, Reexport, CardinalBSplines

if VERSION < v"0.7-"
    ComplexF64 = Complex128
else
    import Base: lastindex, firstindex
    # linspace(a,b,c) = range(a, stop=b, length=c)
end

include("sequences/sequence.jl")
include("filterbanks/filterbank.jl")
include("dwt/discretewavelets.jl")

@reexport using .DWT

end
