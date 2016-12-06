# evaluate.jl

using Polynomials
import Base: eltype
typealias WLT Wavelets.WT

# TODO Find a nice place in wt.jl to fit this in
eltype{W <: WLT.WaveletClass}(w::W) = eltype(typeof(w))
eltype(::Type{WLT.WaveletClass}) = Float64
# Convenience method
eltype(x, y) = promote_type(eltype(x), eltype(y))
"""
  Evaluate the primal scaling function of a wavelet in a point x.
"""
function evaluate_primal_scalingfunction end

function evaluate_primal_scalingfunction{W<:WLT.DiscreteWavelet}(wlt::W, x::Number)
  warn("not implenented")
end

for W in (:(WLT.Daubechies{1}), :(WLT.Haar))
  @eval evaluate_primal_scalingfunction(::$W, x::Number) = _evaluate_constant_Bspline(x, eltype(x, $W))
end

evaluate_primal_scalingfunction{N1,N2}(w::WLT.CDF{N1,N2}, x::Number) =
      evaluate_primal_scalingfunction(N1, x, eltype(w, x))

function evaluate_primal_scalingfunction(N1::Int, x::Number, T::Type)
  if N1 == 1
    _evaluate_constant_Bspline(x, T)
  elseif N1 == 2
    _evaluate_linear_Bspline(x, T)
  elseif N1 == 3
    _evaluate_quadratic_Bspline(x, T)
  elseif N1 == 4
    _evaluate_cubic_Bspline(x, T)
  elseif N1 == 5
    _evaluate_biquadratic_Bspline(x, T)
  else
    T(x)/T(N1-1)*evaluate_primal_scalingfunction(N1-1, x, T) +
        (T(N1)-T(x))/T(N1-1)*evaluate_primal_scalingfunction(N1-1, x-1, T)
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
