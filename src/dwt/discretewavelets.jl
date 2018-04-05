# discretewavelets.jl
module DWT

using ..Sequences
using ..Filterbanks

import ..Filterbanks: Filterbank
import Base: eltype, filter
import ..Sequences: support
export filter


# from wvlt.jl
export is_symmetric, is_orthogonal, is_biorthogonal
export Primal, Dual, scaling, wavelet, coefficient
export support, vanishingmoments, support_length, filter
export evaluate, evaluate_periodic, evaluate_periodic_in_dyadic_points, evaluate_in_dyadic_points
# from DiscreteWavelets.jl
export print_implemented_wavelets, IMPLEMENTED_WAVELETS

include("util/util_functions.jl")
include("boundaries.jl")
include("transform_types.jl")
include("wvlt.jl")

# Transformation of function evaluations on dyadic grid to scaling coefficients
include("scaling_coefficients.jl")

# Transformation of scaling coefficients to wavelet coefficients
include("dwtstep.jl")
include("dwttransform.jl")

# Convenience function
name{T}(::T) = name(T)

include("wvlt_daubechies.jl")
include("wvlt_cdf.jl")

IMPLEMENTED_WAVELETS = (IMPLEMENTED_DB_WAVELETS..., IMPLEMENTED_CDF_WAVELETS...)

print_implemented_wavelets() = println(map(name, IMPLEMENTED_WAVELETS))

include("util/recipes.jl")
end
