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
immutable Prl <: Side end
immutable Dul <: Side end
immutable Scl <: Kind end
immutable Wvl <: Kind end
function name end
primal = Prl(); DWT.name(::Prl) = "primal"
dual = Dul(); DWT.name(::Dul) = "dual"
scaling = Scl(); DWT.name(::Scl) = "scaling"
wavelet = Wvl(); DWT.name(::Wvl) = "wavelet"
Base.inv(::Prl) = Dul()
Base.inv(::Dul) = Prl()
###############################################################################
# vanishingmoments
###############################################################################
vanishingmoments{WT<:DiscreteWavelet}(side::Side, kind::Kind, ::Type{WT}) = vanishingmoments(s, WT)
vanishingmoments{WT<:DiscreteWavelet}(::Prl, ::Type{WT}) = throw("unimplemented")
vanishingmoments{WT<:DiscreteWavelet}(::Dul, W::Type{WT}) = _vanishingmoments(Prl(), W, is_orthogonal(W))
_vanishingmoments(::Prl, W, is_orthogonal::Type{True}) = vanishingmoments(Prl(), W)
###############################################################################
# support/support_length
###############################################################################
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
###############################################################################
# filter
###############################################################################
# By default, the wavelet filters are associated with the dual scaling filters via the alternating flip relation
Base.filter{W<:DiscreteWavelet}(side::Side, kind::Wvl, ::Type{W}) = alternating_flip(filter(inv(side), Scl(), W))

# If orthogonal, dual and primal scaling functions are equal
Base.filter{W<:DiscreteWavelet}(side::Dul, kind::Scl, ::Type{W}) = _filter(side, kind, W, is_orthogonal(W))
_filter{W<:DiscreteWavelet}(::Dul, ::Scl, ::Type{W}, is_orthogonal::Type{True}) = filter(Prl(), Scl(), W)

# coefficient filter is just, √2 times the scaling filter, overwrite if it can have nice (rational) values
type Cof <: Kind end
coefficient = Cof()

Base.filter{W<:DiscreteWavelet}(side::Side, ::Cof, ::Type{W}) = sqrt(eltype(W)(2))*filter(side, Scl(), W)
###############################################################################
# All previous functions where applicable on types, now make methods for instances.
###############################################################################
for op in (:support, :support_length, :filter)
  @eval DWT.$op(side::Side, kind::Kind, w::DiscreteWavelet, args...) = DWT.$op(side, kind, typeof(w), args...)
end

for op in (:vanishingmoments,)
  @eval $op(side::Side, w::DiscreteWavelet, args...) = $op(side, typeof(w), args...)
end
###############################################################################
# Evaluation of scaling and wavelet function
###############################################################################
export eval_periodic, eval_periodic_in_dyadic_points, eval_in_dyadic_points
include("util/cascade.jl")
include("util/periodize.jl")
function eval_periodic{T, S<:Real}(side::Side, kind::Kind, w::DiscreteWavelet{T}, j::Int, k::Int, x::S; xtol::S=1e-5, options...)
  a = T(0); b = T(1)
  offset = floor(x)
  x -= offset
  d, kpoint = closest_dyadic_point(x, xtol; options...)
  (!in_periodic_support(kpoint/2^d, periodic_support(side, kind, w, j, k, a, b)...)) && return T(0)
  f = eval_periodic_in_dyadic_points(side, kind, w, j, k, d)
  f[kpoint+1]
end

#TODO test this eval method
function DWT.eval{T, S<:Real}(side::Side, kind::Scl, w::DiscreteWavelet{T}, j::Int, k::Int, x::S; xtol::S=1e-5, options...)
  d, kpoint = closest_dyadic_point(x, xtol; options...)
  s = support(side, kind, w, j, k)
  (!in_support(kpoint/2^d, s)) && return T(0)
  f = eval_periodic_in_dyadic_points(side, kind, w, j, k, d)
  f[-Int(s[1]*(1<<d))+kpoint+1]
end

function eval_periodic_in_dyadic_points{T}(side::Side, kind::Kind, w::DiscreteWavelet{T}, j=0, k=0, d=10; points=false, options...)
  a = T(0); b = T(1)
  s = eval_in_dyadic_points(side, kind, w, j ,k ,d)
  L = round(Int, (1<<d)*(b-a))
  @assert abs(L-(1<<d)*(b-a)) < eps(real(T))
  f = _periodize(L, s, Int((1<<d)*DWT.support(side, kind, w, j, k)[1])+1)
  if points
    f, linspace(a,b,L+1)[1:end-1]
  else
    f
  end
end

function eval_in_dyadic_points{T}(side::Side, kind::Kind, w::DiscreteWavelet{T}, j=0, k=0, d=10; points=false, options...)
  @assert (d-j) >= 0
  f = eval_in_dyadic_points(filter(side, kind, w), j, k, d; options...)
  if points
    f, dyadicpointsofcascade(side, kind, w, j, k, d; options...)
  else
    f
  end
end

function eval_in_dyadic_points{T}(side::Side, kind::Wvl, w::DiscreteWavelet{T}, j=0, k=0, d=10; points=false, options...)
  s = eval_in_dyadic_points(side, Scl(), w, j+1, k, d; options...)
  filter = DWT.filter(side, Wvl(), w)
  L = DWT.support_length(side, Wvl(), w)
  f = zeros(T, (1<<(d-j))*L+1)
  for (i,l) in enumerate(firstindex(filter):lastindex(filter))
    offset = (1<<(d-1-j))*(i-1)
    f[offset+1:offset+length(s)] += filter[l]*s
  end
  if points
    f, dyadicpointsofcascade(side, kind, w, j, k, d; options...)
  else
    f
  end
end

function eval_in_dyadic_points{T}(s::CompactSequence{T}, j=0, k=0, d=0; options...)
  T(2)^(T(j)/2)*cascade_algorithm(s, (d-j); options...)
end

###############################################################################

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

include("wvlt_daubechies.jl")
include("wvlt_cdf.jl")

IMPLEMENTED_WAVELETS = (IMPLEMENTED_DB_WAVELETS..., IMPLEMENTED_CDF_WAVELETS...)
end
