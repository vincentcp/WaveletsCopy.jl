# disable precompilation during development
#__precompile__()

module Wavelets

using Reexport
include("util.jl")

include("filtering.jl")
include("sequence.jl")

include("wt.jl")
include("transforms.jl")


include("firfilter.jl")
include("filterbank.jl")
include("discretewavelets.jl")

include("threshold.jl")
include("plot.jl")

@reexport using .Util, .WT, .Transforms, .Threshold, .Plot, .Sampling, .Filters, .DWT, .Sequences
#@reexport using .Util, .WT, .Transforms, .Threshold, .Plot

end
