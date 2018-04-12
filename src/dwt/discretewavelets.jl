# discretewavelets.jl
module DWT

using ..Sequences
using ..Filterbanks

import ..Filterbanks: Filterbank
import Base: eltype, filter
import ..Sequences: support
export filter

# from scaling_coefficients.jl

# from boudaries
export perbound, symbound, zerobound


# from wvlt.jl
export is_symmetric, is_orthogonal, is_biorthogonal
export Primal, Dual, scaling, wavelet, coefficient
export support, vanishingmoments, support_length, filter
# from evaluation.jl
export evaluate, evaluate_periodic, evaluate_periodic_in_dyadic_points, evaluate_in_dyadic_points
# from DiscreteWavelets.jl
export print_implemented_wavelets, IMPLEMENTED_WAVELETS, print_all_implemented_wavelets, ALL_IMPLEMENTED_WAVELETS

# from wvlt_cdf.jl
export CDFWavelet
# wvlt_daubechies.jl
export HaarWavelet, DaubechiesWavelet

# from util/waveletindex.jl
export WaveletIndex, wavelet_indices, kind, level, offset, value

# from util/scratchspace.jl
export EvalScratchSpace, EvalPeriodicScratchSpace



include("util/util_functions.jl")

include("boundaries.jl")
include("transform_types.jl")
include("wvlt.jl")


Filterbank(w::DiscreteWavelet) =
    Filterbank( FilterPair(filter(Prl(), Scl(), w), filter(Prl(), Wvl(), w)),
                FilterPair(filter(Dul(), Scl(), w), filter(Dul(), Wvl(), w)) )
DualFilterbank(w::DiscreteWavelet) =
    Filterbank( FilterPair(filter(Dul(), Scl(), w), filter(Dul(), Wvl(), w)),
                FilterPair(filter(Prl(), Scl(), w), filter(Prl(), Wvl(), w)) )

"DWT groups the data that fully characterize a discrete wavelet transform."
struct DWT_Data
    "The wavelet to use in the discrete wavelet transform."
    wavelet ::  DiscreteWavelet
    "The treatment of boundaries."
    bnd     ::  WaveletBoundary
    "The transform type."
    transformtype   ::  TransformType
end


# Transformation of function evaluations on dyadic grid to scaling coefficients
include("scaling_coefficients.jl")

# Transformation of scaling coefficients to wavelet coefficients
include("dwtstep.jl")
include("dwttransform.jl")

# Convenience function
name{T}(::T) = name(T)

include("wvlt_daubechies.jl")
include("wvlt_cdf.jl")

IMPLEMENTED_WAVELETS = (IMPLEMENTED_DB_WAVELETS..., IMPLEMENTED_CDF_WAVELETS_Float64...)
ALL_IMPLEMENTED_WAVELETS = (IMPLEMENTED_DB_WAVELETS...,IMPLEMENTED_CDF_WAVELETS...)

print_implemented_wavelets() = println(map(name, IMPLEMENTED_WAVELETS))
print_all_implemented_wavelets() = println(map(name, ALL_IMPLEMENTED_WAVELETS))

include("util/recipes.jl")
end
