# scaling_coefficients.jl
import ..Util: ndyadicscales, isdyadic

# Function on the interval (a, b) to scaling coefficients
"""
  Transformation of a function on the interval [a, b] to `2^L` scaling coefficients.

  The function is assumed to be periodically extended.
"""
scaling_coefficients(f::Function, w::DWT.DiscreteWavelet, L::Int, boundary::PeriodicBoundary, a::Number=0, b::Number=1; options...) =
    scaling_coefficients(f, w, L, PeriodicEmbedding(), a, b; options...)
"""
  Transformation of a function evaluation on a dyadic grid to `2^L` scaling coefficients
"""
scaling_coefficients{T}(f::AbstractArray{T,1}, w::DiscreteWavelet{T}, boundary::PeriodicBoundary; dual=true, options...) = dual?
      scaling_coefficients(f, _scalingcoefficient_filter(dual_scalingfilter(w)), PeriodicEmbedding(); options...) :
      scaling_coefficients(f, _scalingcoefficient_filter(primal_scalingfilter(w)), PeriodicEmbedding(); options...)


function scaling_coefficients(f::Function, w::DWT.DiscreteWavelet, L::Int, fembedding, a::Number=0, b::Number=1; dual=true, options...)
  T = promote_type(eltype(w), eltype(a), eltype(b))
  a = T(a); b = T(b)
  if dual
    filter = dual_scalingfilter
  else
    filter = primal_scalingfilter
  end
  (b-a)*scaling_coefficients(x->f((b-a)*x+a), filter(w), L::Int, fembedding; options...)
end

# Function on the interval (0 1) to scaling coefficients
function scaling_coefficients{T}(f::Function, s::CompactSequence{T}, L::Int, fembedding; options...)
  x = linspace(T(0), T(1), 1<<L + 1)[1:end-1]
  fcoefs = map(f, x)
  @assert eltype(fcoefs)==T
  filter = _scalingcoefficient_filter(s)
  fembedding == nothing && (fembedding = FunctionEmbedding(k -> f(k*T(2)^(-T(L)))))
  scaling_coefficients(fcoefs, filter, fembedding; options...)
end


# Function evaluations on a dyadic grid to scaling coefficients
function scaling_coefficients{T}(f::AbstractArray{T,1}, filter::CompactSequence{T}, fembedding; n::Int=length(f), options...)
  @assert isdyadic(f)
  c = Array(T,n)
  scaling_coefficients!(c, f, filter, fembedding; options...)
end

function scaling_coefficients!{T}(c, f, filter::CompactSequence{T}, fembedding; offset::Int=0, options...)
  # TODO write a convolution function
  # convolution between low pass filter and function values gives approximation of scaling coefficients
  for j in offset:offset+length(c)-1
    ci = zero(T)
    for l in firstindex(filter):lastindex(filter)
      ci += filter[l]*fembedding[f, j-l]
    end
    c[j+1-offset] = ci
  end
  T(1)/T(sqrt(length(f)))*c
end

_scalingcoefficient_filter(f::CompactSequence) =
    reverse(CompactSequence(cascade_algorithm(f, 0), f.offset))

# Scaling coefficients to function evaluations on dyadic grid (assumes periodic extension)
function scaling_coefficients_to_dyadic_grid{T}(scaling_coefficients::Array{T,1}, w::DWT.DiscreteWavelet{T}, d=ndyadicscales(scaling_coefficients); grid=false, options...)
  @assert isdyadic(scaling_coefficients)
  j = ndyadicscales(scaling_coefficients)

  s = primal_support(w)
  # length of equidistand grid associated to function on [0, 1)
  f_l = 1<<d
  # length of grid (same resolution as before) associated to scaling function in [s[1],s[2]]
  scaling_l = 1+(s[2]-s[1])*(1<<(d-j))

  function_evals = zeros(T,f_l)
  for (c_i,c) in enumerate(scaling_coefficients)
    k = c_i - 1
    istart = 1 + max(0,-(1<<(d-j))*(s[1]+k))
    iend =  -(1<<(d-j))*s[1] + min((1<<d)-(1<<(d-j))*k, 1+(1<<(d-j))*s[2])
    offset = (s[1]+k)*(1<<(d-j))

    # Periodic extension, so wrap scalingfunction outside of [0,1] periodically into the interval
    f = c*DWT.scaling_function_in_dyadic_points(w, j, k, d)
    # Part of the scaling function that lays in [0,1]
    ids = istart:iend
    function_evals[ids+offset] += f[ids]
    # Part of scaling function right of [0,1]
    ids = iend+1:f_l:scaling_l
    for (i,index) in enumerate(ids)
      i == length(ids) ?
        function_evals[1:scaling_l-index+1] += f[index:scaling_l] :
        function_evals += f[index:index+f_l-1]
    end
    # Part of scaling function left of [0,1]
    ids = istart-1:-f_l:1
    for (i,index) in enumerate(ids)
      i == length(ids) ?
        function_evals[end+1-index:end] += f[1:index] :
        function_evals += f[index-f_l+1:index]
    end
  end
  grid ?
    (return function_evals, linspace(T(0), T(1), length(function_evals)+1)[1:end-1]) :
    (return function_evals)
end

function primal_support(n::Int, i::Int, l::Int, w::DiscreteWavelet)
  kind, j, k = wavelet_index(n,i,l)
  if      kind  ==  :scaling
    primal_scalingsupport(w, j, k)
  elseif  kind  ==  :wavelet
    primal_waveletsupport(w, j, k)
  end
end

"""
  The index ([:scaling/:wavelet], j, k) in the (scaling+wavelet) sequence for coefficient i in a sequence of length n after l dwt synthesis_lowpassfilter

  For example, the indices of a sequence with 4 elements after
  0 dwt steps
    (:scaling, 2, 0),    (:scaling, 2, 1),    (:scaling, 2, 2),     (:scaling, 2, 3)
  1 dwt step
    (:scaling, 1, 0),    (:scaling, 1, 1),    (:wavelet, 1, 0),     (:wavelet, 1, 1)
  2 dwt steps
    (:scaling, 0, 0),    (:wavelet, 0, 0),    (:wavelet, 1, 0),     (:wavelet, 1, 1)
"""
function wavelet_index(n::Int, i::Int, l::Int)
  if i > n/(1<<l)
    j = level(n,i)
    k = mod(i-1,1<<j)
    :wavelet, j, k
  else
    :scaling, Int(log2(n))-l, i-1
  end
end

function coefficient_index(kind, j::Int, k::Int)::Int
  if kind == :scaling
    k+1
  elseif kind == :wavelet
    1<<j+k+1
  end
end

function level(n::Int, i::Int)
  (i == 1 || i == 2) && (return 0)
  for l in 1:round(Int,log2(n))
    if i <= (1<<(l+1))
      return l
    end
  end
end
