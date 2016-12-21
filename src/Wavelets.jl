# disable precompilation during development
#__precompile__()

module Wavelets

function name
end

using Reexport
# Code from early Wavelet Toolbox
include("util.jl")
include("wt.jl")
include("transforms.jl")
include("threshold.jl")
include("plot.jl")
@reexport using .Util, .WT, .Transforms, .Threshold, .Plot


include("filtering.jl")
include("sequences/sequence.jl")
include("firfilter.jl")
include("filterbank.jl")
include("dwt/discretewavelets.jl")

@reexport using .DWT, .Sequences, .Filterbanks

end
