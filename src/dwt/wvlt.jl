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
###############################################################################
# Evaluation of scaling and wavelet function
###############################################################################
include("util/recursion.jl")
include("util/periodize.jl")


"""
evaluation of 2^{j/2}ϕ(2^jx-k) where f is the primal/dual scaling/wavelet function of type `w` in `x`. The result is periodized with period 1.

See also `evaluate`
"""
function evaluate_periodic{T, S<:Real}(side::Side, kind::Kind, w::DiscreteWavelet{T}, j::Int, k::Int, x::S; xtol::S=1e-5, options...)
    # Map x to [0,1]
    a = T(0); b = T(1)
    offset = floor(x)
    x -= offset

    # Look for kpoint*2^{-d} close to x
    d, kpoint = closest_dyadic_point(x, xtol; options...)

    # If x is outside of support return 0
    (!in_periodic_support(kpoint/2^d, periodic_support(side, kind, w, j, k, a, b)...)) && return T(0)

    # Evaluate wavelet in dyadic points
    f = evaluate_periodic_in_dyadic_points(side, kind, w, j, k, d)
    # and select the correct points
    f[kpoint+1]
end

"""
Evaluation of 2^{j/2}ϕ(2^jx-k) where f is the primal/dual scaling/wavelet function of type `w` in `x`.

This is only approximate we look for a dyadic point `l2^{-d}` l∈Z, d∈N close to x. (this point can not always be found)
"""
function evaluate{T, S<:Real}(side::Side, kind::Kind, w::DiscreteWavelet{T}, j::Int, k::Int, x::S; xtol::S=1e-5, options...)
    # Look for kpoint*2^{-d} close to x
    d, kpoint = closest_dyadic_point(x, xtol; options...)

    # If x is outside of support return 0
    s = support(side, kind, w, j, k)
    (!in_support(kpoint/2^d, s)) && return T(0)

    # Evaluate wavelet in dyadic points
    f = evaluate_in_dyadic_points(side, kind, w, j, k, d)
    # and select the correct points knowing that zero is always present
    f[-ceil(Int,s[1]*(1<<d))+kpoint+1]
end

"""
Periodic evaluation of 2^{j/2}ϕ(2^jx-k) where f is the primal/dual scaling/wavelet function of type `w` in an equispaced grid with separation `2^-d`.

The result is periodized with period 1.
See also `evaluate_in_dyadic_points`
"""
function evaluate_periodic_in_dyadic_points{T}(side::DWT.Side, kind::DWT.Scl, w::DWT.DiscreteWavelet{T}, j=0, k=0, d=10,
        scratch=nothing, scratch2=nothing, scratch3=nothing; points=false, options...)
    # Wrapper for allocating memory and using evaluate_periodic_in_dyadic_points!

    # allocate the right amount of scratch space
    f = zeros(T, DWT.evaluate_periodic_in_dyadic_points_length(side, kind, w, j, k, d))
    scratch = zeros(T, DWT.evaluate_periodic_in_dyadic_points_scratch_length(side, kind, w, j, k, d))
    scratch2 = zeros(T, DWT.evaluate_periodic_in_dyadic_points_scratch2_length(side, kind, w, j, k, d))

    # Evaluate wavelet in dyadic points
    evaluate_periodic_in_dyadic_points!(f, side, kind, w, j, k, d, scratch, scratch2; options...)

    # Return points also if points if true
    if points
        f, linspace(T(0),T(1),(1<<d)+1)[1:end-1]
    else
        f
    end
end

function evaluate_periodic_in_dyadic_points{T}(side::DWT.Side, kind::DWT.Wvl, w::DWT.DiscreteWavelet{T}, j=0, k=0, d=10,
        scratch=nothing, scratch2=nothing, scratch3=nothing; points=false, options...)
    # Wrapper for allocating memory and using evaluate_periodic_in_dyadic_points!,
    # wavelet evaluation needs more scratch than scaling evaluation.

    # allocate the right amount of scratch space
    f = zeros(T, DWT.evaluate_periodic_in_dyadic_points_length(side, kind, w, j, k, d))
    scratch = zeros(T, DWT.evaluate_periodic_in_dyadic_points_scratch_length(side, kind, w, j, k, d))
    scratch2 = zeros(T, DWT.evaluate_periodic_in_dyadic_points_scratch2_length(side, kind, w, j, k, d))
    scratch3 = zeros(T, DWT.evaluate_periodic_in_dyadic_points_scratch3_length(side, kind, w, j, k, d))

    # Evaluate wavelet in dyadic points
    evaluate_periodic_in_dyadic_points!(f, side, kind, w, j, k, d, scratch, scratch2, scratch3; options...)

    # Return points also if points if true
    if points
        f, linspace(T(0),T(1),(1<<d)+1)[1:end-1]
    else
        f
    end
end


"""
Evaluation of 2^{j/2}ϕ(2^jx-k) where f is the primal/dual scaling/wavelet function of type `w` in an equispaced grid with separation `2^-d`.

If the options `points=true` is added, the equispaced grid is returned after the evaluations
"""
function evaluate_in_dyadic_points{T}(side::DWT.Side, kind::DWT.Kind, w::DWT.DiscreteWavelet{T}, j=0, k=0, d=10;
        points=false, options...)
    # Wrapper for allocating memory and using evaluate_in_dyadic_points!

    # allocate the right amount of scratch space
    f = zeros(T, DWT.evaluate_in_dyadic_points_length(side, kind, w, j, k, d))
    scratch = zeros(T, DWT.evaluate_in_dyadic_points_scratch_length(side, kind, w, j, k, d))

    # Evaluate wavelet in dyadic points
    evaluate_in_dyadic_points!(f, side, kind, w, j, k, d, scratch; options...)

    # Return points also if points if true
    if points
        f, dyadicpointsofrecursion(side, kind, w, j, k, d; options...)
    else
        f
    end
end

function evaluate_in_dyadic_points{T}(side::DWT.Side, kind::DWT.Wvl, w::DWT.DiscreteWavelet{T}, j=0, k=0, d=10;
        points=false, options...)
    # Wrapper for allocating memory and using evaluate_in_dyadic_points!
    # wavelet evaluation needs more scratch than scaling evaluation.

    # allocate the right amount of scratch space
    f = zeros(T, DWT.evaluate_in_dyadic_points_length(side, kind, w, j, k, d))
    scratch = zeros(T, DWT.evaluate_in_dyadic_points_scratch_length(side, kind, w, j, k, d))
    scratch2 = zeros(T, DWT.evaluate_in_dyadic_points_scratch2_length(side, kind, w, j, k, d))

    # Evaluate wavelet in dyadic points
    evaluate_in_dyadic_points!(f, side, kind, w, j, k, d, scratch, scratch2; options...)

    # Return points also if points if true
    if points
        f, dyadicpointsofrecursion(side, kind, w, j, k, d; options...)
    else
        f
    end
end

# In place methods
function evaluate_periodic_in_dyadic_points!{T}(f::AbstractArray{T}, side::DWT.Side, kind::DWT.Wvl, w::DWT.DiscreteWavelet{T}, j=0, k=0, d=10,
        scratch=nothing, scratch2=nothing, scratch3=nothing; options...)
    # Periodic evaluation consist of a usual evaluation and a periodization
    # Wavelet evaluation uses more schratch than scaling evaluation

    DWT.evaluate_in_dyadic_points!(scratch, side, kind, w, j ,k ,d, scratch2, scratch3; options...)
    DWT._periodize!(f, scratch, -ceil(Int,(1<<d)*DWT.support(side, kind, w, j, k)[1])+1)
    # try
    #     DWT._periodize!(f, scratch, -Int((1<<d)*DWT.support(side, kind, w, j, k)[1])+1)
    # catch InexactError
    #     DWT._periodize!(f, scratch, -Int(cld(DWT.support(side, kind, w, 0, k)[1],1<<(j-d))/(1<<d))+1)
    # end
end

function evaluate_periodic_in_dyadic_points!{T}(f::AbstractArray{T}, side::DWT.Side, kind::DWT.Scl, w::DWT.DiscreteWavelet{T}, j=0, k=0, d=10,
        scratch=nothing, scratch2=nothing; options...)
    # Periodic evaluation consist of a usual evaluation and a periodization

    DWT.evaluate_in_dyadic_points!(scratch, side, kind, w, j ,k ,d, scratch2; options...)
    # try
        DWT._periodize!(f, scratch, -ceil(Int,(1<<d)*DWT.support(side, kind, w, j, k)[1])+1)
        # DWT._periodize!(f, scratch, -Int((1<<d)*DWT.support(side, kind, w, j, k)[1])+1)
    # catch InexactError
        # DWT._periodize!(f, scratch, -Int(cld(DWT.support(side, kind, w, 0, k)[1],1<<(j-d))/(1<<d))+1)
    # end
end

function evaluate_in_dyadic_points!{T}(f::AbstractArray{T,1}, side::DWT.Side, kind::DWT.Wvl, w::DWT.DiscreteWavelet{T}, j=0, k=0, d=10,
        scratch = nothing, scratch2 = nothing; options...)
    # Wavelet evaluation uses scaling evaluation on a finer level and a the linear relation between wavelet and scaling functions

    @assert length(f) == DWT.evaluate_in_dyadic_points_length(side, kind, w, j, k, d)

    # if (d-j) >= 1, then a finer level exist
    if (d-j) >= 1
        # evaluate scaling functions in finer level
        DWT.evaluate_in_dyadic_points!(scratch, side, DWT.Scl(), w, j+1, k, d, scratch2; options...)

        # get the wavelet filter
        filter = DWT.filter(side, DWT.Wvl(), w)

        # use the wavelet filter d_n to calculate
        # ψ(t) = ∑_n d_nϕ(2t-n)
        # ψ_jk(t) = 2^{j/2}ψ(2^jt-k)
        # so scratch = ϕ_{j+1,k}(t_d)
        # and ψ_{j,k}(t_d) = ∑_n d_n ϕ_{j+1,k}(t_d-n)
        scratchlength = length(scratch)
        m = 1<<(d-1-j)
        for (i,l) in enumerate(firstindex(filter):lastindex(filter))
            offset = m*(i-1)
            f[offset+1:offset+scratchlength] += filter[l]*scratch
        end
    else
        # finer level does not exist calculate wavelet on finer level
        DWT.evaluate_in_dyadic_points!(scratch2, side, kind, w, j, k, j+1, scratch, nothing; options...)
        # and select the interesting values including zero

        # index of value related with x=0
        istart = -Int((1<<(j+1))*DWT.support(side, kind, w, j, k)[1])+1
        # Look for the first dyadic point in the array
        istart = mod(istart-1, 1<<(j+1-d))+1
        # Copy the interesting values
        copy!(f, scratch2[istart:1<<(j+1-d):end])
        nothing
    end
end

function evaluate_in_dyadic_points!{T}(f::AbstractArray{T,1}, side::Side, kind::Kind, w::DiscreteWavelet{T}, j=0, k=0, d=10,
        scratch=nothing; verbose=true, options...)
    # Scaling evaluation uses scaling evaluation on a finer level and a the linear relation between wavelet and scaling functions

    @assert length(f) == evaluate_in_dyadic_points_length(side, kind, w, j, k, d)

    # If grid is fine enough for the level
    if (d-j) >= 0
        # Do scaling evalutaion using the scaling filter
        evaluate_in_dyadic_points!(f, filter(side, kind, w), j, k, d; points=false, options...)

    # If grid is not fine enough
    else
        @assert length(scratch) == DWT.recursion_length(filter(side, kind, w), 0)

        # do evaluation on a finer grid
        evaluate_in_dyadic_points!(scratch, filter(side, kind, w), j, k, j; options...)
        # index of value related with x=0
        istart = -Int((1<<(j+1))*DWT.support(side, kind, w, j, k)[1])+1
        # Look for the first dyadic point in the array
        istart = mod(istart-1, 1<<(j+1-d))+1
        # Copy the interesting values
        copy!(f, scratch[istart:1<<(j+1-d):end])
    end
end

function evaluate_in_dyadic_points!{T}(f::AbstractArray{T,1}, s::CompactSequence{T}, j=0, k=0, d=0; options...)
    # Evaluation is done through the recursion_algorithm

    @assert length(f) == DWT.recursion_length(s, (d-j))
    recursion_algorithm!(f, s, (d-j); options...)
    scale!(f,T(2)^(T(j)/2))
    nothing
end

## scratch util function #######################################################
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
