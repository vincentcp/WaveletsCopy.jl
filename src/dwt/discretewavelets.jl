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
export filter

# Traits
export is_symmetric, is_orthogonal, is_biorthogonal

export primal_scalingfilter, dual_scalingfilter, primal_support, dual_support
export primal_support_length, dual_support_length
export primal_vanishingmoments, dual_vanishingmoments, primal_waveletfilter, dual_waveletfilter
export primal_coefficientfilter, dual_coefficientfilter

export primal, dual, scaling, wavelet, coefficient
export support, vanishingmoments, support_length, filter

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
immutable TestWavelet{T} <: DWT.DiscreteWavelet{T} end

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
abstract Kind
abstract Side
type Prl <: Side end
type Dul <: Side end
type Scl <: Kind end
type Wvl <: Kind end
primal = Prl()
dual = Dul()
scaling = Scl()
wavelet = Wvl()
Base.inv(::Prl) = Dul()
Base.inv(::Dul) = Prl()
# Vanishing Moments
vanishingmoments{WT<:DiscreteWavelet}(side::Side, kind::Kind, ::Type{WT}) = vanishingmoments(s, WT)
vanishingmoments{WT<:DiscreteWavelet}(::Prl, ::Type{WT}) = throw("unimplemented")
vanishingmoments{WT<:DiscreteWavelet}(::Dul, W::Type{WT}) = _vanishingmoments(Prl(), W, is_orthogonal(W))
_vanishingmoments(::Prl, W, is_orthogonal::Type{True}) = vanishingmoments(Prl(), W)

# Support
DWT.support{WT<:DiscreteWavelet}(side::Side, kind::Scl, ::Type{WT}, j::Int=0, k::Int=0) = (j == 0 && k == 0) ?
    Sequences.support(filter(side, Scl(), WT)) :
    Sequences.support(filter(side, Scl(), WT), j, k)

function DWT.support{WT<:DiscreteWavelet}(side::Side, kind::Wvl, ::Type{WT}, j::Int=0, k::Int=0)
  supp1 = support(side, Scl(), WT)
  supp2 = support(inv(side), Scl(), WT)
  S1 = Int(1/2*(supp1[1]-supp2[2]+1))
  S2 = Int(1/2*(supp1[2]-supp2[1]+1))
  (j == 0 && k == 0) ? (S1,S2) : (1/(1<<j)*(S1[1]+k), 1/(1<<j)*(S2+k))
end

support_length{WT<:DiscreteWavelet}(side::Side, kind::Kind,  ::Type{WT}) = support(side, kind, WT)[2] - support(side, kind, WT)[1]

# Filters

# By default, the wavelet filters are associated with the dual scaling filters via the alternating flip relation
Base.filter{W<:DiscreteWavelet}(side::Side, kind::Wvl, ::Type{W}) = alternating_flip(filter(inv(side), Scl(), W))

# If orthogonal, dual and primal scaling functions are equal
Base.filter{W<:DiscreteWavelet}(side::Dul, kind::Scl, ::Type{W}) = _filter(side, kind, W, is_orthogonal(W))
_filter{W<:DiscreteWavelet}(::Dul, ::Scl, ::Type{W}, is_orthogonal::Type{True}) = filter(Prl(), Scl(), W)

# coefficient filter is just, √2 times the scaling filter, overwrite if it can have nice (rational) values
type Cof <: Kind end
coefficient = Cof()

Base.filter{W<:DiscreteWavelet}(side::Side, ::Cof, ::Type{W}) = sqrt(eltype(W)(2))*filter(side, Scl(), W)


evaluate_function{W<:DiscreteWavelet}(::Side, ::Kind, ::Type{W}) = throw("unimplemented")
# TODO implement default evaluate_function(side, wavelet, w) using wavelet coefficients and evaluate_function(side, scaling, w)

# All previous functions where applicable on types, now make methods for instances.
for op in (:support, :support_length, :filter, :evaluate_function)
  @eval DWT.$op(side::Side, kind::Kind, w::DiscreteWavelet, args...) = DWT.$op(side, kind, typeof(w), args...)
end

for op in (:vanishingmoments,)
  @eval $op(side::Side, w::DiscreteWavelet, args...) = $op(side, typeof(w), args...)
end

Filterbank(w::DiscreteWavelet) =
    Filterbank( FilterPair(filter(Prl(), Scl(), w), filter(Prl(), Wvl(), w)),
                FilterPair(filter(Dul(), Scl(), w), filter(Dul(), Wvl(), w)) )

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

include("util/cardinal_b_splines.jl")
include("util/cascade.jl")

include("wvlt_daubechies.jl")
include("wvlt_cdf.jl")

IMPLEMENTED_WAVELETS = (IMPLEMENTED_DB_WAVELETS..., IMPLEMENTED_CDF_WAVELETS...)
end
