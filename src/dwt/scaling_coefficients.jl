# scaling_coefficients.jl

# Function on the interval (a, b) to scaling coefficients
function scaling_coefficients(f::Function, w::DWT.DiscreteWavelet, fembedding, a=0, b=1; options...)
  T = promote_type(eltype(w), eltype(a), eltype(b))
  a = T(a); b = T(b)
  (b-a)*scaling_coefficients(x->f((b-a)*x+a), dual_scalingfilter(w), fembedding; options...)
end

# Function on the interval (0 1) to scaling coefficients
function scaling_coefficients{T}(f::Function, s::CompactSequence{T}, fembedding; L=4, options...)
  x = linspace(T(0), T(1), 1<<L + 1)[1:end-1]
  fcoefs = map(f, x)
  @assert eltype(fcoefs)==T
  filter = _scalingcoefficient_filter(s)
  fembedding == nothing && (fembedding = FunctionEmbedding(k -> f(k*T(2)^(-T(L)))))
  scaling_coefficients(fcoefs, filter, fembedding; options...)
end

scaling_coefficients{T}(f::AbstractArray{T,1}, w::DiscreteWavelet{T}, fembedding; options...) =
      scaling_coefficients(f, dual_scalingfilter(w), fembedding; options...)
# Function evaluations on a dyadic grid to scaling coefficients
function scaling_coefficients{T}(f::AbstractArray{T,1}, filter::CompactSequence{T}, fembedding; n::Int=length(f), options...)
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
    reverse(CompactSequence(cascade_algorithm(f), f.offset))

# # Scaling coefficients to function evaluations on dyadic grid
# function scaling_coefficients_to_dyadic_grid{T}(scaling_coefficients::Array{T,1}, w::DWT.DiscreteWavelet{T}, a=0, b=1)
#   function_evals = similar(scaling_coefficients)
# end
