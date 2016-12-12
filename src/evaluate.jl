# evaluate.jl

using Polynomials
using .Sequences
using .Filters
import Base: eltype

typealias WLT Union{WT.WaveletClass,DWT.DiscreteWavelet}

# Convenience method
eltype(x, y) = promote_type(eltype(x), eltype(y))
"""
  Evaluate the primal scaling function of a wavelet in a point x.

  ϕ(x)
"""
function evaluate_primal_scalingfunction end
"""
  Evaluate the shifted and dilated scaling function of a wavelet in a point x.

  ϕ_jk = 2^(k/2) ϕ(2^k-j)
"""
evaluate_transformed_primal_scalingfunction{W<:WLT}(wlt::W, x::Number, k::Int=0, j::Int=0) =
    2.0^(k/2)*evaluate_primal_scalingfunction(wlt, 2.0^k-j)

function evaluate_primal_scalingfunction{W<:WLT}(wlt::W, x::Number)
  error("not implemented")
end

for W in (:(WT.Daubechies{1}), :(WT.Haar))
  @eval evaluate_primal_scalingfunction(::$W, x::Number) = evaluate_Bspline(Degree{0}, x, eltype(x, $W))
end

evaluate_primal_scalingfunction{T}(w::DWT.HaarWavelet{T}, x::Number) =
      evaluate_Bspline(Degree{0}, x, eltype(x, w))

evaluate_primal_scalingfunction{N1,N2}(w::WT.CDF{N1,N2}, x::Number) =
      evaluate_Bspline(Degree{N1-1}, x-T(DWT.symmetric_offset(N1+1)), eltype(w, x))

evaluate_primal_scalingfunction{N1,N2,T}(w::DWT.CDFWavelet{N1,N2,T}, x::Number) =
      evaluate_Bspline(Degree{N1-1}, x-T(DWT.symmetric_offset(N1+1)), eltype(w, x))

# Implementation of cardinal B splines of degree N
typealias Degree{N} Val{N}

evaluate_Bspline(N::Int, x, T::Type) = evaluate_Bspline(Degree{N}, x, T)
function evaluate_Bspline{N}(::Type{Degree{N}}, x, T::Type)
  T(x)/T(N)*evaluate_Bspline(Degree{N-1}, x, T) +
      (T(N+1)-T(x))/T(N)*evaluate_Bspline(Degree{N-1}, x-1, T)
end

evaluate_Bspline(::Type{Degree{0}}, x, T::Type) = (0 <= x < 1) ? T(1) : T(0)

function evaluate_Bspline(::Type{Degree{1}}, x, T::Type)
  if (0 <= x < 1)
    s = T[0, 1]
  elseif (1 <= x < 2)
    s = T[2, -1]
  else
    s = T[0, 0]
  end
  @eval @evalpoly $x $(s...)
end
function evaluate_Bspline(::Type{Degree{2}}, x, T::Type)
  if (0 <= x < 1)
    s = T[0, 0, 1/2]
  elseif (1 <= x < 2)
    s = T[-3/2, 3, -1]
  elseif (2 <= x < 3)
    s = T[9/2, -3, 1/2]
  else
    s = T[0, 0]
  end
  @eval @evalpoly $x $(s...)
end

function evaluate_Bspline(::Type{Degree{3}}, x, T::Type)
  if (0 <= x < 1)
    s = T[0, 0, 0, 1/6]
  elseif (1 <= x < 2)
    s = T[2/3, -2, 2, -1/2]
  elseif (2 <= x < 3)
    s = T[-22/3, 10, -4, 1/2]
  elseif (3 <= x < 4)
    s = T[32/3, -8, 2, -1/6]
  else
    s = T[0, 0]
  end
  @eval @evalpoly $x $(s...)
end

function evaluate_Bspline(::Type{Degree{4}}, x, T::Type)
  if (0 <= x < 1)
    s = T[0, 0, 0, 0, 1/24]
  elseif (1 <= x < 2)
    s = T[-5/24, 5/6, -5/4, 5/6, -1/6]
  elseif (2 <= x < 3)
    s = T[155/24, -25/2, 35/4, -5/2, 1/4]
  elseif (3 <= x < 4)
    s = T[-655/24, 65/2, -55/4, 5/2, -1/6]
  elseif (4 <= x < 5)
    s = T[625/24, -125/6, 25/4, -5/6, 1/24]
  else
    s = T[0, 0]
  end
  @eval @evalpoly $x $(s...)
end

"""
The scaling function of a wavelet with filtercoeficients h in diadic points

`x[k] = 2^(-L)k for k=0..K`

where K is the length of h. Thus L=0, gives the evaluation of the
scaling function in the points [0,1,...,K], and L=1, the points [0,.5,1,...,K].
"""
cascade_algorithm(w::DWT.DiscreteWavelet; dual=true, options...) = dual?
    cascade_algorithm(DWT.dual_scalingfilter(w); options...) :
    cascade_algorithm(DWT.primal_scalingfilter(w); options...)
cascade_algorithm(s::CompactSequence; options...) = cascade_algorithm(s.a; options...)
function cascade_algorithm(h; L = 0, tol = 1e-8, options...)
  T = eltype(h)
  @assert sum(h)≈sqrt(T(2))
  N = length(h)
  # find ϕ(0), .. , ϕ(N-1) by solving a eigenvalue problem
  # The eigenvector with eigenvalue equal to 1/√2 is the vector containting ϕ(0), .. , ϕ(N-1)
  # see http://cnx.org/contents/0nnvPGYf@7/Computing-the-Scaling-Function for more mathematics

  # Create matrix from filter coefficients
  H = _get_H(h, N, T)
  # Find eigenvector eigv
  E = eigfact(H)
  index = find(abs(E[:values]-1/sqrt(T(2))).<tol)
  @assert length(index) > 0
  for i in index
    V = E[:vectors][:,i]
    @assert norm(imag(V)) < tol
    eigv = real(V)
    abs(sum(eigv)) > tol*100 ? break : nothing
  end
  eigv /= sum(eigv)

  # Find intermediate values ϕ(1/2), .. ,ϕ(N-1 -1/2)
  K = 2
  for m in 1:L
    eigv = vcat(eigv, zeros(T,length(eigv)-1))
    eigv_length = length(eigv)
    for n in eigv_length:-2:3
      eigv[n] = eigv[ceil(Int,n/2)]
    end
    eigv[2:2:end] = 0
    for n in 2:2:eigv_length
      I = max(1,ceil(Int,(2n-1-eigv_length)/K+1)):min(N,floor(Int, (2n-1)/K+1))
      eigv[n] = sqrt(T(2))*sum(h[i].*eigv[2*n-1-K*(i-1)] for i in I)
    end
    K = K << 1
  end
  eigv
end

dyadicpointsofcascade(s::DWT.DiscreteWavelet, L::Int=10) = dyadicpointsofcascade(DWT.primal_scalingfilter(s), L)
dyadicpointsofcascade(s::CompactSequence, L::Int=10) = dyadicpointsofcascade(s.a, L, s.offset)
dyadicpointsofcascade(h, L::Int=10, offset::Int=0) = dyadicpointsofcascade(eltype(h), length(h), L, offset)
dyadicpointsofcascade(T::Type, H::Int, L::Int, offset::Int=0) = L==0 ?
      T(offset)+linspace(T(0),T(H-1),H) : T(offset)+linspace(T(0),T(H-1),2^L*(H-1)+1)

function _get_H(h, N, T)
  M = zeros(T,N,N)
  for l in 1:2:N
    hh = h[l]
    for m in 1:ceil(Int,N/2)
      M[l>>1+m,1+(m-1)*2] = hh
    end
  end
  for l in 2:2:N
    hh = h[l]
    for m in 1:N>>1
      M[l>>1+m,2+(m-1)*2] = hh
    end
  end
  M
end

# On the interval 0 1

scaling_coefficients(f::Function, w::DWT.DiscreteWavelet, fembedding; options...) =
      scaling_coefficients(f, DWT.dual_scalingfilter(w), fembedding; options...)
function scaling_coefficients(f::Function, s::CompactSequence, fembedding; L=4, options...)
  filter = _scalingcoefficient_filter(s)
  x = linspace(0, 1, 1<<L + 1)
  x = x[1:end-1]; f2 = k -> f(k*2.0^(-L))
  # iseven(filter.n) ? (x = midpoints(x); f2 = k-> f((k+.5)*2.0^(-L)) ) :
  #   (x = x[1:end-1]; f2 = k -> f(k*2.0^(-L)))



  fcoefs = map(f, x)
  T = promote_type(eltype(fcoefs), eltype(s))
  fcoefs = convert(Array{T}, fcoefs)
  fembedding == nothing && (fembedding = FunctionEmbedding(f2))
  scaling_coefficients(fcoefs, filter, fembedding)
end

function scaling_coefficients{T}(f, filter::CompactSequence{T}, fembedding)
  @assert eltype(f)==T
  c = similar(f)
  scaling_coefficients!(c, f, filter, fembedding)
end
function scaling_coefficients!{T}(c, f, filter::CompactSequence{T}, fembedding)
  # TODO write a convolution function
  # convolution between low pass filter and function values gives approximation of scaling coefficients
  for j in 0:length(c)-1
    ci = zero(T)
    for l in firstindex(filter):lastindex(filter)
      ci += filter[l]*fembedding[f, j-l]
    end
    c[j+1] = ci
  end
  T(1)/T(sqrt(length(c)))*c
end

function _scalingcoefficient_filter(f::CompactSequence)
  #TODO remove zeros
  # phi_eval = reverse(cascade_algorithm(f))
  # println(phi_eval)
  # F = length(phi_eval)
  # # remove zeros and shift the resulting CompactSequence accordingly
  # @assert phi_eval[1] < 1e-14
  # @assert F == 2 || phi_eval[end] < 1e-14
  # F == 2 ? phi_eval = phi_eval[2:2] : phi_eval = phi_eval[2:end-1]
  # filter = CompactSequence(phi_eval, f.offset+1)
  reverse(CompactSequence(cascade_algorithm(f), f.offset))
end

_scalingcoefficient_filter(DWT.primal_scalingfilter(DWT.cdf22))
#
# DWT.dual_scalingfilter(DWT.db1)
# cascade_algorithm(DWT.dual_scalingfilter(DWT.db1))
# _scalingcoefficient_filter(DWT.db1)
