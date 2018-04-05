abstract type DiscreteWavelet{T} end
struct TestWavelet{T} <: DWT.DiscreteWavelet{T} end

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
abstract type Kind end
abstract type Side end
struct Prl <: Side end
struct Dul <: Side end
struct Scl <: Kind end
struct Wvl <: Kind end
function name end
Primal = Prl(); DWT.name(::Prl) = "primal"
Dual = Dul(); DWT.name(::Dul) = "dual"
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
support{WT<:DiscreteWavelet}(side::Side, kind::Scl, ::Type{WT}, j::Int=0, k::Int=0) = (j == 0 && k == 0) ?
    Sequences.support(filter(side, Scl(), WT)) :
    Sequences.support(filter(side, Scl(), WT), j, k)

function support{WT<:DiscreteWavelet}(side::Side, kind::Wvl, ::Type{WT}, j::Int=0, k::Int=0)
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
filter{W<:DiscreteWavelet}(side::Side, kind::Wvl, ::Type{W}) = alternating_flip(filter(inv(side), Scl(), W))

# If orthogonal, dual and primal scaling functions are equal
filter{W<:DiscreteWavelet}(side::Dul, kind::Scl, ::Type{W}) = _filter(side, kind, W, is_orthogonal(W))
_filter{W<:DiscreteWavelet}(::Dul, ::Scl, ::Type{W}, is_orthogonal::Type{True}) = filter(Prl(), Scl(), W)

# coefficient filter is just, √2 times the scaling filter, overwrite if it can have nice (rational) values
type Cof <: Kind end
coefficient = Cof()

filter{W<:DiscreteWavelet}(side::Side, ::Cof, ::Type{W}) = sqrt(eltype(W)(2))*filter(side, Scl(), W)
###############################################################################
# All previous functions where applicable on types, now make methods for instances.
###############################################################################
for op in (:support, :support_length, :filter)
  @eval $op(side::Side, kind::Kind, w::DiscreteWavelet, args...) = $op(side, kind, typeof(w), args...)
end

for op in (:vanishingmoments,)
  @eval $op(side::Side, w::DiscreteWavelet, args...) = $op(side, typeof(w), args...)
end
###############################################################################
# Evaluation of scaling and wavelet function
###############################################################################
include("util/cascade.jl")
include("util/periodize.jl")
function evaluate_periodic{T, S<:Real}(side::Side, kind::Kind, w::DiscreteWavelet{T}, j::Int, k::Int, x::S; xtol::S=1e-5, options...)
  a = T(0); b = T(1)
  offset = floor(x)
  x -= offset
  d, kpoint = closest_dyadic_point(x, xtol; options...)
  (!in_periodic_support(kpoint/2^d, periodic_support(side, kind, w, j, k, a, b)...)) && return T(0)
  f = evaluate_periodic_in_dyadic_points(side, kind, w, j, k, d)
  f[kpoint+1]
end

function evaluate{T, S<:Real}(side::Side, kind::Scl, w::DiscreteWavelet{T}, j::Int, k::Int, x::S; xtol::S=1e-5, options...)
  d, kpoint = closest_dyadic_point(x, xtol; options...)
  s = support(side, kind, w, j, k)
  (!in_support(kpoint/2^d, s)) && return T(0)
  f = evaluate_periodic_in_dyadic_points(side, kind, w, j, k, d)
  f[-Int(s[1]*(1<<d))+kpoint+1]
end

function evaluate_periodic_in_dyadic_points{T}(side::DWT.Side, kind::DWT.Scl, w::DWT.DiscreteWavelet{T}, j=0, k=0, d=10, scratch=nothing, scratch2=nothing, scratch3=nothing; points=false, options...)
  f = zeros(T, DWT.evaluate_periodic_in_dyadic_points_length(side, kind, w, j, k, d))
  scratch = zeros(T, DWT.evaluate_periodic_in_dyadic_points_scratch_length(side, kind, w, j, k, d))
  scratch2 = zeros(T, DWT.evaluate_periodic_in_dyadic_points_scratch2_length(side, kind, w, j, k, d))
  DWT.evaluate_periodic_in_dyadic_points!(f, side, kind, w, j, k, d, scratch, scratch2; options...)
  if points
    f, linspace(T(0),T(1),(1<<d)+1)[1:end-1]
  else
    f
  end
end

function evaluate_periodic_in_dyadic_points{T}(side::DWT.Side, kind::DWT.Wvl, w::DWT.DiscreteWavelet{T}, j=0, k=0, d=10, scratch=nothing, scratch2=nothing, scratch3=nothing; points=false, options...)
  f = zeros(T, DWT.evaluate_periodic_in_dyadic_points_length(side, kind, w, j, k, d))
  scratch = zeros(T, DWT.evaluate_periodic_in_dyadic_points_scratch_length(side, kind, w, j, k, d))
  scratch2 = zeros(T, DWT.evaluate_periodic_in_dyadic_points_scratch2_length(side, kind, w, j, k, d))
  scratch3 = zeros(T, DWT.evaluate_periodic_in_dyadic_points_scratch3_length(side, kind, w, j, k, d))
  DWT.evaluate_periodic_in_dyadic_points!(f, side, kind, w, j, k, d, scratch, scratch2, scratch3; options...)
  if points
    f, linspace(T(0),T(1),(1<<d)+1)[1:end-1]
  else
    f
  end
end

function evaluate_in_dyadic_points{T}(side::DWT.Side, kind::DWT.Kind, w::DWT.DiscreteWavelet{T}, j=0, k=0, d=10; points=false, options...)
  f = zeros(T, DWT.evaluate_in_dyadic_points_length(side, kind, w, j, k, d))
  scratch = zeros(T, DWT.evaluate_in_dyadic_points_scratch_length(side, kind, w, j, k, d))
  DWT.evaluate_in_dyadic_points!(f, side, kind, w, j, k, d, scratch; options...)
  if points
    f, DWT.dyadicpointsofcascade(side, kind, w, j, k, d; options...)
  else
    f
  end
end

function evaluate_in_dyadic_points{T}(side::DWT.Side, kind::DWT.Wvl, w::DWT.DiscreteWavelet{T}, j=0, k=0, d=10; points=false, options...)
  f = zeros(T, DWT.evaluate_in_dyadic_points_length(side, kind, w, j, k, d))
  scratch = zeros(T, DWT.evaluate_in_dyadic_points_scratch_length(side, kind, w, j, k, d))
  scratch2 = zeros(T, DWT.evaluate_in_dyadic_points_scratch2_length(side, kind, w, j, k, d))
  DWT.evaluate_in_dyadic_points!(f, side, kind, w, j, k, d, scratch, scratch2; options...)
  if points
    f, DWT.dyadicpointsofcascade(side, kind, w, j, k, d; options...)
  else
    f
  end
end

# In place methods
function evaluate_periodic_in_dyadic_points!{T}(f::AbstractArray{T}, side::DWT.Side, kind::DWT.Wvl, w::DWT.DiscreteWavelet{T}, j=0, k=0, d=10, scratch=nothing, scratch2=nothing, scratch3=nothing; options...)
  DWT.evaluate_in_dyadic_points!(scratch, side, kind, w, j ,k ,d, scratch2, scratch3; options...)
  DWT._periodize!(f, scratch, Int((1<<d)*DWT.support(side, kind, w, j, k)[1])+1)
end

function evaluate_periodic_in_dyadic_points!{T}(f::AbstractArray{T}, side::DWT.Side, kind::DWT.Scl, w::DWT.DiscreteWavelet{T}, j=0, k=0, d=10, scratch=nothing, scratch2=nothing; options...)
  DWT.evaluate_in_dyadic_points!(scratch, side, kind, w, j ,k ,d, scratch2; options...)
  DWT._periodize!(f, scratch, Int((1<<d)*DWT.support(side, kind, w, j, k)[1])+1)
end

function evaluate_in_dyadic_points!{T}(f::AbstractArray{T,1}, side::DWT.Side, kind::DWT.Wvl, w::DWT.DiscreteWavelet{T}, j=0, k=0, d=10, scratch = nothing, scratch2 = nothing; options...)
  @assert length(f) == DWT.evaluate_in_dyadic_points_length(side, kind, w, j, k, d)
  if (d-j) >= 1
    DWT.evaluate_in_dyadic_points!(scratch, side, DWT.Scl(), w, j+1, k, d, scratch2; options...)
    scratchlength = length(scratch)
    filter = DWT.filter(side, DWT.Wvl(), w)
    L = DWT.support_length(side, DWT.Wvl(), w)
    m = 1<<(d-1-j)
    for (i,l) in enumerate(firstindex(filter):lastindex(filter))
      offset = m*(i-1)
      f[offset+1:offset+scratchlength] += filter[l]*scratch
    end
  else
    DWT.evaluate_in_dyadic_points!(scratch2, side, kind, w, j, k, j+1, scratch, nothing; options...)
    copy!(f, scratch2[1:1<<(j+1-d):end])
    nothing
  end
end

function evaluate_in_dyadic_points!{T}(f::AbstractArray{T,1}, side::Side, kind::Kind, w::DiscreteWavelet{T}, j=0, k=0, d=10, scratch=nothing; verbose=true, options...)
  @assert length(f) == evaluate_in_dyadic_points_length(side, kind, w, j, k, d)
  if (d-j) >= 0
    evaluate_in_dyadic_points!(f, filter(side, kind, w), j, k, d; points=false, options...)
  else
    @assert length(scratch) == DWT.cascade_length(filter(side, kind, w), 0)
    evaluate_in_dyadic_points!(scratch, filter(side, kind, w), j, k, j; options...)
    copy!(f, scratch[1:1<<(j-d):end])
  end
end

function evaluate_in_dyadic_points!{T}(f::AbstractArray{T,1}, s::CompactSequence{T}, j=0, k=0, d=0; options...)
  @assert length(f) == DWT.cascade_length(s, (d-j))
  cascade_algorithm!(f, s, (d-j); options...)
  scale!(f,T(2)^(T(j)/2))
  nothing
end

evaluate_periodic_in_dyadic_points_length(side::DWT.Side, kind::DWT.Kind, w, j, k, d) =
    evaluate_periodic_in_dyadic_points_length(d)
evaluate_periodic_in_dyadic_points_length(d) = 1<<d

evaluate_periodic_in_dyadic_points_scratch_length(side::DWT.Side, kind::DWT.Kind, w, j, k, d) =
    evaluate_in_dyadic_points_length(side, kind, w, j, k, d)
evaluate_periodic_in_dyadic_points_scratch2_length(side::DWT.Side, kind::DWT.Kind, w, j, k, d) =
    evaluate_in_dyadic_points_scratch_length(side, kind, w, j, k, d)
evaluate_periodic_in_dyadic_points_scratch3_length(side::DWT.Side, kind::DWT.Kind, w, j, k, d) =
    evaluate_in_dyadic_points_scratch2_length(side, kind, w, j, k, d)

evaluate_in_dyadic_points_length(side::DWT.Side, kind::DWT.Kind, w, j, k, d) = (d-j) >= 0 ?
    (1<<(d-j))*DWT.support_length(side, kind, w)+1 :
    cld(DWT.support_length(side, kind, w)+1, 1<<(j-d))

evaluate_in_dyadic_points_scratch_length(side::DWT.Side, kind::DWT.Kind, w, j, k, d) = (d-j) >= 0 ?
    0 :
    DWT.support_length(side, kind, w)+1
evaluate_in_dyadic_points_scratch_length(side::DWT.Side, kind::DWT.Wvl, w, j, k, d) = (d-j) >= 1 ?
    evaluate_in_dyadic_points_length(side, DWT.Scl(), w, j+1, k, d):
    evaluate_in_dyadic_points_scratch_length(side, kind, w, j, k, j+1)

evaluate_in_dyadic_points_scratch2_length(side::DWT.Side, kind::DWT.Wvl, w, j, k, d) = (d-j) >= 1 ?
    evaluate_in_dyadic_points_scratch_length(side, DWT.Scl(), w, j+1, k, d) :
    evaluate_in_dyadic_points_length(side, kind, w, j, k, j+1)

###############################################################################

Filterbank(w::DiscreteWavelet) =
    Filterbank( FilterPair(filter(Prl(), Scl(), w), filter(Prl(), Wvl(), w)),
                FilterPair(filter(Dul(), Scl(), w), filter(Dul(), Wvl(), w)) )

"DWT groups the data that fully characterize a discrete wavelet transform."
struct DWT_Data
    "The wavelet to use in the discrete wavelet transform."
    wavelet ::  DiscreteWavelet
    "The treatment of boundaries."
    bnd     ::  WaveletBoundary
    "The transform type."
    transformtype   ::  TransformType
end
