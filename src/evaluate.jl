# evaluate.jl

using Polynomials
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
evaluate_primal_scalingfunction{W<:WLT}(wlt::W, k::Int, j::Int, x::Number) =
    2.0^(k/2)*evaluate_primal_scalingfunction(wlt, 2.0^k-j)

function evaluate_primal_scalingfunction{W<:WLT}(wlt::W, x::Number)
  error("not implemented")
end

for W in (:(WT.Daubechies{1}), :(WT.Haar))
  @eval evaluate_primal_scalingfunction(::$W, x::Number) = _evaluate_constant_Bspline(x, eltype(x, $W))
end

evaluate_primal_scalingfunction{T}(w::DWT.HaarWavelet{T}, x::Number) =
    _evaluate_constant_Bspline(x, eltype(x, w))

evaluate_primal_scalingfunction{T}(w::DWT.CDFWavelet{1,1,T}, x::Number) =
    _evaluate_constant_Bspline(x, eltype(x, w))

evaluate_primal_scalingfunction{N1,N2}(w::WT.CDF{N1,N2}, x::Number) =
      evaluate_cardinalBSpline(N1, x, eltype(w, x))

evaluate_primal_scalingfunction{N1,N2,T}(w::DWT.CDFWavelet{N1,N2,T}, x::Number) =
      evaluate_cardinalBSpline(N1, x, eltype(w, x))

# Implementation of cardinal B splines of order m
function evaluate_cardinalBSpline(m::Int, x::Number, T::Type)
  if m <= 0
    error("Order of cardinal B spline should be stricly positive")
  elseif m == 1
    _evaluate_constant_Bspline(x, T)
  elseif m == 2
    _evaluate_linear_Bspline(x, T)
  elseif m == 3
    _evaluate_quadratic_Bspline(x, T)
  elseif m == 4
    _evaluate_cubic_Bspline(x, T)
  elseif m == 5
    _evaluate_biquadratic_Bspline(x, T)
  else
    T(x)/T(m-1)*evaluate_cardinalBSpline(m-1, x, T) +
        (T(m)-T(x))/T(m-1)*evaluate_cardinalBSpline(m-1, x-1, T)
  end
end

_evaluate_constant_Bspline(x, T::Type) = (0 <= x < 1) ? T(1) : T(0)

function _evaluate_linear_Bspline(x, T::Type)
  if (0 <= x < 1)
    s = T[0, 1]
  elseif (1 <= x < 2)
    s = T[2, -1]
  else
    s = T[0, 0]
  end
  @eval @evalpoly $x $(s...)
end

function _evaluate_quadratic_Bspline(x, T::Type)
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

function _evaluate_cubic_Bspline(x, T::Type)
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

function _evaluate_biquadratic_Bspline(x, T::Type)
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
function cascade_algorithm(h, L = 10; tol = 1e-8, options...)
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
  index = index[1]
  V = E[:vectors][:,index]
  @assert norm(imag(V)) < tol
  eigv = real(V)
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

dyadicpointsofcascade(T::Type,H::Int,L::Int) = L==0 ?
      linspace(T(0),T(H-1),H) : linspace(T(0),T(H-1),2^L*(H-1)+1)

dyadicpointsofcascade(h, L::Int=10) = dyadicpointsofcascade(eltype(h), length(h), L)

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
