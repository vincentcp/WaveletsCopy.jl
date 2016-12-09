# discretewavelets.jl
module DWT

using ..Sequences
using ..Filterbanks

typealias False Val{false}
typealias True Val{true}

import ..Filterbanks: Filterbank
import Base: eltype
import Wavelets: name

export filter

# Traits
export is_symmetric, is_orthogonal, is_biorthogonal

export name, class

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
# Not sure yet whether this one makes sense:
#is_semiorthogonal{W <: DiscreteWavelet}(::Type{W}) = False
eltype{T}(::Type{DiscreteWavelet{T}}) = T
eltype{W <: DiscreteWavelet}(::Type{W}) = eltype(supertype(W))
eltype(w::DiscreteWavelet) = eltype(typeof(w))

for op in (:is_symmetric, :is_orthogonal, :is_biorthogonal)
    @eval $op(w::DiscreteWavelet) = $op(typeof(w))()
end


analysis_lowpassfilter(w::DiscreteWavelet) = primal_scalingfilter(w)
analysis_highpassfilter(w::DiscreteWavelet) = primal_waveletfilter(w)

synthesis_lowpassfilter(w::DiscreteWavelet) = dual_scalingfilter(w)
synthesis_highpassfilter(w::DiscreteWavelet) = dual_waveletfilter(w)

# By default, the wavelet filters are associated with the dual scaling filters via the alternating flip relation
dual_waveletfilter(w) = alternating_flip(primal_scalingfilter(w))
primal_waveletfilter(w) = alternating_flip(dual_scalingfilter(w))


dual_scalingfilter(w::DiscreteWavelet) = _dual_scalingfilter(w, is_orthogonal(w))
_dual_scalingfilter(w, is_orthogonal::True) = primal_scalingfilter(w)

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


include("dwtstep.jl")
include("dwttransform.jl")
# Convenience function
name{T}(::T) = name(T)
name{T}(::Type{T}) = "Not defined"
include("wvlt_daubechies.jl")
include("wvlt_cdf.jl")



end # module
