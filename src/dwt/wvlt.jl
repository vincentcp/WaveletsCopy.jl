abstract type DiscreteWavelet{T} end
struct TestWavelet{T} <: DWT.DiscreteWavelet{T} end

# Symmetry trait
is_symmetric(::Type{W}) where {W <: DiscreteWavelet} = False

# Orthogonality traits
is_orthogonal(::Type{W}) where {W <: DiscreteWavelet} = False
is_biorthogonal(::Type{W}) where {W <: DiscreteWavelet} = False
is_semiorthogonal(::Type{W}) where {W <: DiscreteWavelet} = False
# Not sure yet whether this one makes sense:
#is_semiorthogonal{W <: DiscreteWavelet}(::Type{W}) = False
eltype(::Type{DiscreteWavelet{T}}) where {T} = T
eltype(::Type{W}) where {W <: DiscreteWavelet} = eltype(supertype(W))
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
if !(VERSION < v"0.7-")
    Base.broadcastable(s::Side) = Ref(s)
    Base.broadcastable(k::Kind) = Ref(k)
    Base.broadcastable(w::DiscreteWavelet{T}) where T = Ref(w)
end


include("util/waveletindex.jl")

###############################################################################
# vanishingmoments
###############################################################################
vanishingmoments(::Prl, ::Type{WT})  where {WT <: DiscreteWavelet}= throw("unimplemented")
vanishingmoments(::Dul, W::Type{WT})  where {WT <: DiscreteWavelet}= _vanishingmoments(Prl(), W, is_orthogonal(W))
_vanishingmoments(::Prl, W, is_orthogonal::Type{True}) = vanishingmoments(Prl(), W)
###############################################################################
# support/support_length
###############################################################################
support(side::Side, kind::Scl, ::Type{WT})  where {WT <: DiscreteWavelet}=
    Sequences.support(filter(side, Scl(), WT))
function support(side::Side, kind::Wvl, ::Type{WT}) where {WT <: DiscreteWavelet}
    l1, r1 = support(side, Scl(), WT)
    l2, r2 = support(inv(side), Scl(), WT)
    ((l1-r2+1)>>1, (r1-l2+1)>>1)
end
function support(side::Side, kind::Kind, ::Type{WT}, j::Int, k::Int) where {WT <: DiscreteWavelet}
    l, r = support(side, kind, WT)
    (1/(1<<j)*(l+k), 1/(1<<j)*(r+k))
end

support_length(side::Side, kind::Kind,  ::Type{WT})  where {WT <: DiscreteWavelet}= support(side, kind, WT)[2] - support(side, kind, WT)[1]
###############################################################################
# filter
###############################################################################
# By default, the wavelet filters are associated with the dual scaling filters via the alternating flip relation
filter(side::Side, kind::Wvl, ::Type{W}) where {W<:DiscreteWavelet} = alternating_flip(filter(inv(side), Scl(), W))

# If orthogonal, dual and primal scaling functions are equal
filter(side::Dul, kind::Scl, ::Type{W}) where {W<:DiscreteWavelet} = _filter(side, kind, W, is_orthogonal(W))
_filter(::Dul, ::Scl, ::Type{W}, is_orthogonal::Type{True}) where {W<:DiscreteWavelet} = filter(Prl(), Scl(), W)

# coefficient filter is just, √2 times the scaling filter, overwrite if it can have nice (rational) values
struct Cof <: Kind end
coefficient = Cof()

filter(side::Side, ::Cof, ::Type{W}) where {W<:DiscreteWavelet} = sqrt(eltype(W)(2))*filter(side, Scl(), W)
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
