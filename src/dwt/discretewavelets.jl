# discretewavelets.jl
module DWT

using ..Sequences
using ..Filterbanks

typealias False Val{false}
typealias True Val{true}


import ..Filterbanks: Filterbank
import Base: eltype
# Convenience method
eltype(x, y) = promote_type(eltype(x), eltype(y))
import Wavelets: name

export filter

# Traits
export is_symmetric, is_orthogonal, is_biorthogonal

export name, class

export primal_scalingfilter, dual_scalingfilter, primal_support, dual_support
export primal_vanishingmoments, dual_vanishingmoments, primal_waveletfilter, dual_waveletfilter
export primal_coefficientfilter, dual_coefficientfilter

# BOUNDARY TYPES

abstract WaveletBoundary

# Periodic boundary condition
immutable PeriodicBoundary <: WaveletBoundary
end

# Symmetric extension
immutable SymmetricBoundary <: WaveletBoundary
end

# zero padding
immutable ZeropaddingBoundary <: WaveletBoundary
end

# constant padding
immutable CPBoundary{T} <: WaveletBoundary
    constant    ::  T
end

perbound = DWT.PeriodicBoundary()
symbound = DWT.SymmetricBoundary()
zerobound = DWT.ZeropaddingBoundary()

# TRANSFORM TYPES

abstract TransformType

immutable FilterbankTransform
end

immutable LiftingTransform
end



abstract DiscreteWavelet{T}

# Symmetry trait
is_symmetric{W <: DiscreteWavelet}(::Type{W}) = False

# Orthogonality traits
is_orthogonal{W <: DiscreteWavelet}(::Type{W}) = False
is_biorthogonal{W <: DiscreteWavelet}(::Type{W}) = False
is_semiorthogonal{W <: DiscreteWavelet}(::Type{W}) = False
# Not sure yet whether this one makes sense:
#is_semiorthogonal{W <: DiscreteWavelet}(::Type{W}) = False
eltype{T}(::Type{DiscreteWavelet{T}}) = T
eltype{W <: DiscreteWavelet}(::Type{W}) = eltype(supertype(W))
eltype(w::DiscreteWavelet) = eltype(typeof(w))

for op in (:is_symmetric, :is_orthogonal, :is_biorthogonal, :is_semiorthogonal)
    @eval $op(w::DiscreteWavelet) = $op(typeof(w))()
end

# Vanishing Moments
primal_vanishingmoments{T}(w::DiscreteWavelet{T}) = primal_vanishingmoments(typeof(w))
dual_vanishingmoments{T}(w::DiscreteWavelet{T}) = dual_vanishingmoments(typeof(w))

primal_vanishingmoments{WT<:DiscreteWavelet}(::Type{WT}) = error("primal_vanishingmoments not implemented for wavelet ", WT)
dual_vanishingmoments{WT<:DiscreteWavelet}(W::Type{WT}) = _primal_vanishingmoments(W, is_orthogonal(W))
_primal_vanishingmoments(W, is_orthogonal::Type{True}) = primal_vanishingmoments(W)

# Support
primal_support{WT<:DiscreteWavelet}(::Type{WT}) = error("primal_support not implemented for wavelet ", WT)
dual_support{WT<:DiscreteWavelet}(W::Type{WT}) = _primal_support(W, is_orthogonal(W))
_primal_support(W, is_orthogonal::Type{True}) = primal_support(W)

primal_support{T}(w::DiscreteWavelet{T}) = primal_support(typeof(w))
dual_support{T}(w::DiscreteWavelet{T}) = dual_support(typeof(w))

# Filters
analysis_lowpassfilter(w::DiscreteWavelet) = primal_scalingfilter(w)
analysis_highpassfilter(w::DiscreteWavelet) = primal_waveletfilter(w)

synthesis_lowpassfilter(w::DiscreteWavelet) = dual_scalingfilter(w)
synthesis_highpassfilter(w::DiscreteWavelet) = dual_waveletfilter(w)

# By default, the wavelet filters are associated with the dual scaling filters via the alternating flip relation
dual_waveletfilter(w) = alternating_flip(primal_scalingfilter(w))
primal_waveletfilter(w) = alternating_flip(dual_scalingfilter(w))


dual_scalingfilter(w::DiscreteWavelet) = _dual_scalingfilter(w, is_orthogonal(w))
_dual_scalingfilter(w, is_orthogonal::True) = primal_scalingfilter(w)

dual_coefficientfilter(w::DiscreteWavelet) = _dual_coefficientfilter(w, is_orthogonal(w))
_dual_coefficientfilter(w, is_orthogonal::True) = primal_coefficientfilter(w)


"""
  Evaluate the shifted and dilated scaling function of a wavelet in a point x.

  ϕ_jk = 2^(k/2) ϕ(2^k-j)
"""
evaluate_transformed_primal_scalingfunction{W<:DiscreteWavelet}(wlt::W, x::Number, k::Int=0, j::Int=0) =
    2.0^(k/2)*evaluate_primal_scalingfunction(wlt, 2.0^k-j)

"""
  Evaluate the primal scaling function of a wavelet in a point x.

  ϕ(x)
"""
function evaluate_primal_scalingfunction end
evaluate_primal_scalingfunction{W<:DiscreteWavelet}(wlt::W, x::Number) = error("No explicit formula of scaling function provided for wavelet: ", wlt)
evaluate_dual_scalingfunction{W<:DiscreteWavelet}(wlt::W, x::Number) =
    _evaluate_dual_scaling_function(wlt, x, is_orthogonal(wlt))
_evaluate_dual_scaling_function{W<:DiscreteWavelet}(wlt::W, x::Number, is_orthogonal::True) = evaluate_primal_scalingfunction(wlt, x)

Filterbank(w::DiscreteWavelet) =
    Filterbank( FilterPair(primal_scalingfilter(w), primal_waveletfilter(w)),
                FilterPair(  dual_scalingfilter(w),   dual_waveletfilter(w)) )

"DWT groups the data that fully characterize a discrete wavelet transform."
immutable DWT_Data
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
name{T}(::Type{T}) = "Not defined"
include("util/cardinal_b_splines.jl")
include("util/cascade.jl")

include("wvlt_daubechies.jl")
include("wvlt_cdf.jl")

IMPLEMENTED_WAVELETS = (IMPLEMENTED_DB_WAVELETS..., IMPLEMENTED_CDF_WAVELETS...)
end
