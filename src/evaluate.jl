# evaluate.jl

using Polynomials

typealias WLT Wavelets.WT
"""
  Evaluate the primal scaling function of a wavelet in a point x.
"""
function evaluate_primal_scalingfunction end

# function evaluate_primal_scalingfunction{W<:WLT.DiscreteWavelet}(wlt::W, x::Number)
#   warn("not implenemted")
# end

for W in (:(WLT.Daubechies{1}), :(WLT.Haar))
  @eval evaluate_primal_scalingfunction(::$W, x::Number) = _evaluate_constant_Bspline(x, eltype(x))
end
_evaluate_constant_Bspline(x, T::Type) = (0 <= x < 1) ? T(1) : T(0)
