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
Base.inv(::Type{Prl}) = Dul
Base.inv(::Type{Dul}) = Prl


include("util/waveletindex.jl")

###############################################################################
# vanishingmoments
###############################################################################
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

include("util/recursion.jl")
include("util/cascade.jl")
include("util/periodize.jl")

include("evaluation.jl")
