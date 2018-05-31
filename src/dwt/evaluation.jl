
include("util/scratchspace.jl")

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
    # and select the correct point
    f[mod(kpoint,length(f))+1]
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

periodic_dyadic_points(d::Int, ::Type{T}=Float64) where {T} = linspace(T(0),T(1),(1<<d)+1)[1:end-1]

"""
Periodic evaluation of 2^{j/2}ϕ(2^jx-k) where f is the primal/dual scaling/wavelet function of type `w` in an equispaced grid with separation `2^-d`.

The result is periodized with period 1.
See also `evaluate_in_dyadic_points`
"""
function evaluate_periodic_in_dyadic_points{T}(side::DWT.Side, kind::DWT.Scl, w::DWT.DiscreteWavelet{T}, j=0, k=0, d=10,
        scratch=nothing, scratch2=nothing, scratch3=nothing; points=false, options...)
    # Wrapper for allocating memory and using evaluate_periodic_in_dyadic_points!

    # allocate the right amount of scratch space
    f, scratch, scratch2 = get_evaluate_periodic_scratch_space(side, kind, w, j, k, d)

    # Evaluate wavelet in dyadic points
    evaluate_periodic_in_dyadic_points!(f, side, kind, w, j, k, d, scratch, scratch2; options...)

    # Return points also if points if true
    if points
        f, periodic_dyadic_points(d, T)
    else
        f
    end
end

function evaluate_periodic_in_dyadic_points{T}(side::DWT.Side, kind::DWT.Wvl, w::DWT.DiscreteWavelet{T}, j=0, k=0, d=10,
        scratch=nothing, scratch2=nothing, scratch3=nothing; points=false, options...)
    # Wrapper for allocating memory and using evaluate_periodic_in_dyadic_points!,
    # wavelet evaluation needs more scratch than scaling evaluation.

    # allocate the right amount of scratch space
    f, scratch, scratch2, scratch3 = get_evaluate_periodic_scratch_space(side, kind, w, j, k, d)

    # Evaluate wavelet in dyadic points
    evaluate_periodic_in_dyadic_points!(f, side, kind, w, j, k, d, scratch, scratch2, scratch3; options...)

    # Return points also if points if true
    if points
        f, periodic_dyadic_points(d, T)
    else
        f
    end
end

"""
Evaluation of 2^{j/2}ϕ(2^jx-k) where f is the primal/dual scaling/wavelet function of type `w` in an equispaced grid with separation `2^-d`.

If the options `points=true` is added, the equispaced grid is returned after the evaluations
"""
function evaluate_in_dyadic_points{T}(side::DWT.Side, kind::DWT.Scl, w::DWT.DiscreteWavelet{T}, j=0, k=0, d=10;
        points=false, options...)
    # Wrapper for allocating memory and using evaluate_in_dyadic_points!

    # allocate the right amount of scratch space
    f, scratch = get_evaluate_scratch_space(side, kind, w, j, k, d)

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
    f, scratch, scratch2 = get_evaluate_scratch_space(side, kind, w, j, k, d)

    # Evaluate wavelet in dyadic points
    evaluate_in_dyadic_points!(f, side, kind, w, j, k, d, scratch, scratch2; options...)

    # Return points also if points if true
    if points
        f, dyadicpointsofrecursion(side, kind, w, j, k, d; options...)
    else
        f
    end
end

evaluate_periodic_wavelet_basis_in_dyadic_points(s::DWT.Side, w::DWT.DiscreteWavelet{T}, coeffs, d::Int=10; options...) where {T} =
    evaluate_periodic_scaling_basis_in_dyadic_points(s, w, DWT.idwt(coeffs, w, perbound), d; options...)

function evaluate_periodic_scaling_basis_in_dyadic_points(s::DWT.Side, w::DWT.DiscreteWavelet{T}, coeffs, d::Int=10; options...) where {T}
    f = zeros(T, evaluate_periodic_in_dyadic_points_scratch_length(s, scaling, w, Int(log2(length(coeffs))), 0, d))
    scratch = zeros(T, evaluate_periodic_in_dyadic_points_scratch2_length(s, scaling, w, Int(log2(length(coeffs))), 0, d))
    f_scaled = similar(f)
    evaluate_periodic_wavelet_basis_in_dyadic_points!(zeros(T, 1<<d), s, w, coeffs, d, f, f_scaled, scratch; options...)
end

evaluate_periodic_wavelet_basis_in_dyadic_points!(y::Vector{T}, s::DWT.Side, w::DWT.DiscreteWavelet{T}, coeffs, d, f, f_scaled, scratch; options...) where T=
    evaluate_periodic_scaling_basis_in_dyadic_points!(y, s, w, coeffs, d, f, f_scaled, scratch; options...)

function evaluate_periodic_scaling_basis_in_dyadic_points!(y::Vector{T}, s::DWT.Side, w::DWT.DiscreteWavelet{T}, coeffs, d::Int, f, f_scaled, scratch;
    points = false,
    options...) where {T}
    j = Int(log2(length(coeffs)))
    DWT.evaluate_in_dyadic_points!(f, s, scaling, w, j, 0, d, scratch)
    f_scaled = similar(f)
    for k in 0:(1<<j)-1
        f_scaled .= f .* coeffs[k+1]
        DWT._periodize_add!(y, f_scaled, -ceil(Int,(1<<d)*DWT.support(s, scaling, w, j, k)[1])+1)
    end
    if points
        y, DWT.periodic_dyadic_points(d, T)
    else
        y
    end
end




# Not correctly implented
# function evaluate_periodic_in_dyadic_points(side::DWT.Side, w::DWT.DiscreteWavelet{T}, coeffs, d::Int=10;
#         points=false, options...) where {T}
#     f = cascade_algorithm_linear_combination(side, w, coeffs, d, perbound)
#     if points
#         f, periodic_dyadic_points(d, T)
#     else
#         f
#     end
# end
#
# inv_evaluate_periodic_in_dyadic_points(side::DWT.Side, w::DWT.DiscreteWavelet{T}, coeffs, d::Int=10) where {T} =
#     inv_cascade_algorithm_linear_combination(side, w, coeffs, d, perbound)



###############################################################################
# In place methods ############################################################
###############################################################################
evaluate_in_dyadic_points!(f::AbstractArray{T,1}, side::Side, kind::Kind, w::DiscreteWavelet{T}, j::Int, k::Int, d::Int, SS::ScratchSpace; options...) where {T} =
    evaluate_in_dyadic_points!(f, side, kind, w, j, k, d, get_evaluate_scratch_space(SS, side, kind, w, j, k, d)[2:end]...; options...)

evaluate_periodic_in_dyadic_points!(f::AbstractArray{T,1}, side::Side, kind::Kind, w::DiscreteWavelet{T}, j::Int, k::Int, d::Int, SS::ScratchSpace; options...) where {T} =
    evaluate_periodic_in_dyadic_points!(f, side, kind, w, j, k, d, get_evaluate_periodic_scratch_space(SS, side, kind, w, j, k, d)[2:end]...; options...)


function evaluate_periodic_in_dyadic_points!{T}(f::AbstractArray{T}, side::DWT.Side, kind::DWT.Wvl, w::DWT.DiscreteWavelet{T}, j::Int, k::Int, d::Int,
        scratch, scratch2, scratch3; options...)
    # Periodic evaluation consist of a usual evaluation and a periodization
    # Wavelet evaluation uses more schratch than scaling evaluation

    DWT.evaluate_in_dyadic_points!(scratch, side, kind, w, j ,k ,d, scratch2, scratch3; options...)
    DWT._periodize!(f, scratch, -ceil(Int,(1<<d)*DWT.support(side, kind, w, j, k)[1])+1)
end

function evaluate_periodic_in_dyadic_points!{T}(f::AbstractArray{T}, side::DWT.Side, kind::DWT.Scl, w::DWT.DiscreteWavelet{T}, j::Int, k::Int, d::Int,
        scratch, scratch2; options...)
    # Periodic evaluation consist of a usual evaluation and a periodization

    DWT.evaluate_in_dyadic_points!(scratch, side, kind, w, j ,k ,d, scratch2; options...)
    DWT._periodize!(f, scratch, -ceil(Int,(1<<d)*DWT.support(side, kind, w, j, k)[1])+1)
end

function evaluate_in_dyadic_points!{T}(f::AbstractArray{T,1}, side::DWT.Side, kind::DWT.Wvl, w::DWT.DiscreteWavelet{T}, j::Int, k::Int, d::Int,
        scratch, scratch2 ; options...)
    # Wavelet evaluation uses scaling evaluation on a finer level and a the linear relation between wavelet and scaling functions

    @assert length(f) == DWT.evaluate_in_dyadic_points_length(side, kind, w, j, k, d)
    f[:] = 0
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

function evaluate_in_dyadic_points!{T}(f::AbstractArray{T,1}, side::Side, kind::Kind, w::DiscreteWavelet{T}, j::Int, k::Int, d::Int,
        scratch; verbose=true, options...)
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

function evaluate_in_dyadic_points!{T}(f::AbstractArray{T,1}, s::CompactSequence{T}, j::Int, k::Int, d::Int; options...)
    # Evaluation is done through the recursion_algorithm

    @assert length(f) == DWT.recursion_length(s, (d-j))
    recursion_algorithm!(f, s, (d-j); options...)
    scale!(f,T(2)^(T(j)/2))
    nothing
end
