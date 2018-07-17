# scaling_coefficients.jl

"""
    The scaling coeffients function deals with the wavelt crime.
    The input of the DWT are inner products between dual functions and the function `f`
    ```
    v_Nk = <ϕ̃_Nk, f>
    ```
    or also the scaling coefficients of the projection of `f`, `P_Nf`, on `V_N`
    ```
    P_Nf(t) = \\sum_k v_Nkϕ_Nk(t)
    ```
    not function samples
    ```
    f_Nk = f(kΔt), Δt = 2^{-N}, N = 2^L
    ```
    Replacing `v_Nk` with `f_k` is the wavelet crime.

    However, assuming that `L` is large and hence that `ϕ` and `ϕ̃` have a narrow support we can think of the
    `v_Nk`s as scaled function values `f_k`
    ```
    v_Nk ≈\\sqrt{Δt/θ}f_Nk, θ=\\sqrt{2π}
    ```

    Then
    ```
    P_Nf(t) = 1/θ\\sum_l f(t-lΔt)ϕ(l) + O(Δt)
    ```
    but convergence is slow...

    Therefore, we use approximants ṽ_Nk for v_Nk using the Riemann sum

    ```
    ṽ_Nk = 1/\\sqrt{2π} \\sum_n f_n2^{N/2}\\overline{ϕ̃(2^NnΔt-k)}Δt
         = 2^{-N/2}/\\sqrt{2π}\\sum_n f_n\\overline{ϕ̃(n-k)}
    ```
    which is a convolution between function values and \\{\\overline{ϕ̃(-n)}\\}
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
function scaling_coefficients(f::Function, s::CompactSequence{T}, L::Int, fembedding; options...) where {T}
      x = linspace(T(0), T(1), 1<<L + 1)[1:end-1]
      fcoefs = map(f, x)
      @assert eltype(fcoefs)==T
      filter = _scalingcoefficient_filter(s)
      fembedding == nothing && (fembedding = FunctionEmbedding(k -> f(k*T(2)^(-T(L)))))
      scaling_coefficients(fcoefs, filter, fembedding; options...)
end

# Convenience function: wavelet to filter
scaling_coefficients(f::AbstractArray, s::Side, w::DiscreteWavelet{T}, bnd::PeriodicBoundary; options...)  where {T}=
    scaling_coefficients(f, _scalingcoefficient_filter(filter(inv(s), Scl(), w)), PeriodicEmbedding(); options...)

# function samples to scaling coeffients
function scaling_coefficients(f::AbstractArray, filter::CompactSequence{T}, fembedding; n::Int=length(f), options...) where {T}
      @assert isdyadic(f)
      c = (VERSION<v"0.7-") ? Array{T}(n) : Array{T}(undef, n)
      scaling_coefficients!(c, f, filter, fembedding; options...)
      c
end

"In place method of scaling_coefficients"
function scaling_coefficients!(c, f, filter::CompactSequence{T}, fembedding; offset::Int=0, options...) where {T}
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

"Transforms scaling coeffients back to function evaluations on the dyadic grid."
function scaling_coefficients_to_dyadic_grid(scaling_coefficients::AbstractArray{T,1}, s::Side, w::DWT.DiscreteWavelet{T}, bnd::DWT.PeriodicBoundary, d=ndyadicscales(scaling_coefficients); grid=false, options...) where {T}
    function_evals = zeros(T,DWT.scaling_coefficients_to_dyadic_grid_length(s,w,d))
    j = Int(log2(length(scaling_coefficients)))
    f = zeros(T, evaluate_periodic_in_dyadic_points_scratch_length(s, scaling, w, j, 0, d))
    scratch = zeros(T, evaluate_periodic_in_dyadic_points_scratch2_length(s, scaling, w, j, 0, d))
    f_scaled = similar(f)

    scaling_coefficients_to_dyadic_grid!(function_evals, scaling_coefficients, w, bnd, f, f_scaled, scratch; options...)
    grid ?
        (return function_evals, linspace(T(0), T(1), length(function_evals)+1)[1:end-1]) :
        (return function_evals)
end

# Scaling coefficients to function evaluations on dyadic grid (assumes periodic extension)
"In place method of scaling_coefficients_to_dyadic_grid."
scaling_coefficients_to_dyadic_grid!(function_evals::AbstractArray{T,1}, scaling_coefficients::AbstractArray{T,1}, w::DWT.DiscreteWavelet{T}, bnd::DWT.PeriodicBoundary, f, f_scaled, scratch; options...)  where {T} =
    evaluate_periodic_scaling_basis_in_dyadic_points!(function_evals, Primal, w, scaling_coefficients, Int(log2(length(function_evals))), f, f_scaled, scratch; options...)

scaling_coefficients_to_dyadic_grid_length(s::Side, w::DWT.DiscreteWavelet, d::Int) = 1<<d

scaling_coefficients_to_dyadic_grid_scratch_length(s::Side, w::DWT.DiscreteWavelet, d::Int) =
    scaling_coefficients_to_dyadic_grid_length(s, w, d)
scaling_coefficients_to_dyadic_grid_scratch2_length(s::Side, w::DWT.DiscreteWavelet, d::Int) =
    DWT.evaluate_periodic_in_dyadic_points_scratch_length(s, DWT.Scl(), w, d, 0, d)

function DWT.support(side::Side, n::Int, i::Int, l::Int, w::DiscreteWavelet)
    kind, j, k = wavelet_index(n,i,l)
    support(side, kind, w, j, k)
end
