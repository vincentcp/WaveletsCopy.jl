# scaling_coefficients.jl

"""
    The scaling coeffients function deals with the wavelt crime.
    The input of the DWT are inner products between dual functions and the function `f`
    ```
    v_Nk = <ϕ̃_Nk, f>
    ```
    or also the scaling coefficients of the projection of `f`, `P_Nf`, on `V_N`
    ```
    P_Nf(t) = \sum_k v_Nkϕ_Nk(t)
    ```
    not function samples
    ```
    f_Nk = f(kΔt), Δt = 2^{-N}, N = 2^L
    ```
    Replacing `v_Nk` with `f_k` is the wavelet crime.

    However, assuming that `L` is large and hence that `ϕ` and `ϕ̃` have a narrow support we can think of the
    `v_Nk`s as scaled function values `f_k`
    ```
    v_Nk ≈\sqrt{Δt/θ}f_Nk, θ=\sqrt{2π}
    ```

    Then
    ```
    P_Nf(t) = 1/θ\sum_l f(t-lΔt)ϕ(l) + O(Δt)
    ```
    but convergence is slow...

    Therefore, we use approximants ṽ_Nk for v_Nk using the Riemann sum

    ```
    ṽ_Nk = 1/\sqrt{2π} \sum_n f_n2^{N/2}\overline{ϕ̃(2^NnΔt-k)}Δt
         = 2^{-N/2}/\sqrt{2π}\sum_n f_n\overline{ϕ̃(n-k)}
    ```
    which is a convolution between function values and \{\overline{ϕ̃(-n)}\}
    which is in turn a low pass filtering.
"""
function scaling_coefficients end

scaling_coefficients(f::Function, w::DWT.DiscreteWavelet, L::Int, fembedding, a::Real=0, b::Real=1; options...) =
    scaling_coefficients(f, Dual, w, L, fembedding, a, b; options...)


# Function on the interval [a,b] to function on [0,1]
function scaling_coefficients(f::Function, side::DWT.Side, w::DWT.DiscreteWavelet, L::Int, fembedding, a::Real=0, b::Real=1; options...)
      T = promote_type(eltype(w), eltype(a), eltype(b))
      a = T(a); b = T(b)
      flt = filter(side, Scl(), w)
      (b-a)*scaling_coefficients(x->f((b-a)*x+a), flt, L::Int, fembedding; options...)
end

# Function on the interval [a,b] to samples
function scaling_coefficients{T}(f::Function, s::CompactSequence{T}, L::Int, fembedding; options...)
      x = linspace(T(0), T(1), 1<<L + 1)[1:end-1]
      fcoefs = map(f, x)
      @assert eltype(fcoefs)==T
      filter = _scalingcoefficient_filter(s)
      fembedding == nothing && (fembedding = FunctionEmbedding(k -> f(k*T(2)^(-T(L)))))
      scaling_coefficients(fcoefs, filter, fembedding; options...)
end

# Convenience function: wavelet to filter
scaling_coefficients{T}(f::AbstractArray, s::Side, w::DiscreteWavelet{T}, bnd::PeriodicBoundary; options...) =
    scaling_coefficients(f, _scalingcoefficient_filter(filter(inv(s), Scl(), w)), PeriodicEmbedding(); options...)

# function samples to scaling coeffients
function scaling_coefficients{T}(f::AbstractArray, filter::CompactSequence{T}, fembedding; n::Int=length(f), options...)
      @assert isdyadic(f)
      c = Array{T}(n)
      scaling_coefficients!(c, f, filter, fembedding; options...)
      c
end

"In place method of scaling_coefficients"
function scaling_coefficients!{T}(c, f, filter::CompactSequence{T}, fembedding; offset::Int=0, options...)
    # TODO write a convolution function
    # convolution between low pass filter and function values gives approximation of scaling coefficients
    for j in offset:offset+length(c)-1
        ci = zero(T)
        for l in firstindex(filter):lastindex(filter)
            ci += filter[l]*fembedding[f, j-l]
        end
        c[j+1-offset] = T(1)/T(sqrt(length(f)))*ci
    end
end

_scalingcoefficient_filter(f::CompactSequence) =
    reverse(CompactSequence(recursion_algorithm(f, 0), f.offset))

# Remove
# "Transforms scaling coeffients back to function evaluations on the dyadic grid."
# scaling_coefficients_to_dyadic_grid{T}(scaling_coefficients::AbstractArray, w::DWT.DiscreteWavelet{T}, bnd::WaveletBoundary, d=ndyadicscales(scaling_coefficients); options...) =
    # scaling_coefficients_to_dyadic_grid(scaling_coefficients, Primal, w, bmd, d; options...)

"Transforms scaling coeffients back to function evaluations on the dyadic grid."
function scaling_coefficients_to_dyadic_grid{T}(scaling_coefficients::AbstractArray, s::Side, w::DWT.DiscreteWavelet{T}, bnd::WaveletBoundary, d=ndyadicscales(scaling_coefficients); grid=false, options...)
    function_evals = zeros(T,DWT.scaling_coefficients_to_dyadic_grid_length(s,w,d))
    scratch = zeros(DWT.scaling_coefficients_to_dyadic_grid_scratch_length(s,w,d))
    scratch2 = zeros(DWT.scaling_coefficients_to_dyadic_grid_scratch2_length(s,w,d))

    scaling_coefficients_to_dyadic_grid!(function_evals, scaling_coefficients, w, bnd, scratch, scratch2; options...)
    grid ?
        (return function_evals, linspace(T(0), T(1), length(function_evals)+1)[1:end-1]) :
        (return function_evals)
end

# Scaling coefficients to function evaluations on dyadic grid (assumes periodic extension)
"In place method of scaling_coefficients_to_dyadic_grid."
scaling_coefficients_to_dyadic_grid!{T}(function_evals::AbstractArray, scaling_coefficients::AbstractArray, w::DWT.DiscreteWavelet{T}, bnd::DWT.PeriodicBoundary, scratch=nothing, scratch2=nothing; options...) =
    scaling_coefficients_to_dyadic_grid!(function_evals, scaling_coefficients, Primal, w, bnd, scratch, scratch2; options...)


"In place method of scaling_coefficients_to_dyadic_grid."
function scaling_coefficients_to_dyadic_grid!{T}(function_evals::AbstractArray, scaling_coefficients::AbstractArray, s::Side, w::DWT.DiscreteWavelet{T}, ::DWT.PeriodicBoundary, scratch=nothing, scratch2=nothing; grid=false, options...)
    d = Int(log2(length(function_evals)))
    @assert length(function_evals) == length(scratch)
    @assert DWT.scaling_coefficients_to_dyadic_grid_scratch2_length(s, w, d) == length(scratch2)
    @assert isdyadic(scaling_coefficients)
    j = ndyadicscales(scaling_coefficients)
    function_evals[:] = T(0)
    for (c_i,c) in enumerate(scaling_coefficients)
        k = c_i - 1
        DWT.evaluate_periodic_in_dyadic_points!(scratch, s, DWT.Scl(), w, j, k, d, scratch2, nothing)
        scale!(scratch, c)
        for i in 1:length(function_evals)
            function_evals[i] += scratch[i]
        end
        # function_evals[:] += c*DWT.evaluate_periodic_in_dyadic_points(Prl(), Scl(), w, j, k, d)
    end
    nothing
end

scaling_coefficients_to_dyadic_grid_length(s::Side, w::DWT.DiscreteWavelet, d::Int) = 1<<d

scaling_coefficients_to_dyadic_grid_scratch_length(s::Side, w::DWT.DiscreteWavelet, d::Int) =
    scaling_coefficients_to_dyadic_grid_length(s, w, d)
scaling_coefficients_to_dyadic_grid_scratch2_length(s::Side, w::DWT.DiscreteWavelet, d::Int) =
    DWT.evaluate_periodic_in_dyadic_points_scratch_length(s, DWT.Scl(), w, d, 0, d)

function DWT.support(side::Side, n::Int, i::Int, l::Int, w::DiscreteWavelet)
    kind, j, k = wavelet_index(n,i,l)
    support(side, kind, w, j, k)
end

" All wavelet indices on a certain level"
function wavelet_indices(l::Int)
    n = 1<<l
    [wavelet_index(n, i, l) for i in 1:n]
end

"""
  The index ([:scaling/:wavelet], j, k) in the (scaling+wavelet) sequence for coefficient i in a sequence of length n after l dwt synthesis_lowpassfilter

  For example, the indices of a sequence with 4 elements after
  0 dwt steps
    (Scl(), 2, 0),    (Scl(), 2, 1),    (Scl(), 2, 2),     (Scl(), 2, 3)
  1 dwt step
    (Scl(), 1, 0),    (Scl(), 1, 1),    (Wvl(), 1, 0),     (Wvl(), 1, 1)
  2 dwt steps
    (Scl(), 0, 0),    (Wvl(), 0, 0),    (Wvl(), 1, 0),     (Wvl(), 1, 1)
"""
function wavelet_index(n::Int, i::Int, l::Int)
    if i > n/(1<<l)
        j = level(n,i)
        k = mod(i-1,1<<j)
        Wvl(), j, k
    else
        Scl(), Int(log2(n))-l, i-1
    end
end

"""
  The coefficient corresponding to the wavelet index (kind, j, k). The integer
  is the index in the array of the DWT.
"""
coefficient_index(kind::Scl, j::Int, k::Int)::Int = k+1
coefficient_index(kind::Wvl, j::Int, k::Int)::Int = 1<<j+k+1

"Return the wavelet level of a coefficient in an array of length `n` with index `i`"
function level(n::Int, i::Int)
    (i == 1 || i == 2) && (return 0)
    for l in 1:round(Int,log2(n))
        if i <= (1<<(l+1))
            return l
        end
    end
end
