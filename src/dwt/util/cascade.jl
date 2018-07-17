
# function cascade_algorithm(side::Side, w::DiscreteWavelet{T}, wind::WaveletIndex, d::Int, bnd::WaveletBoundary) where {T}
#     coeffs = zeros{T}(1<<d)
#     index = coefficient_index(wind)
#     # If the wind does not correspond with an index in `wavelet_indices(1<<d)`: error
#     @assert index < (1<<d)
#     coeffs[index] = 1
#     cascade_algorithm_linear_combination(side, w, coeffs, bnd)
# end

"""
    Evaluation of the wavelet/scaling function with wavelet index `wind` in the
    grid `[k2^-d for k in 0:2^d-1]`. An extension of type `bnd` is assumed. For
    example, a boundary type equal to `perbound` will lead to a periodic
    function with period `1`.
"""
function cascade_algorithm_linear_combination(side::Side, w::DiscreteWavelet{T}, coeffs, d::Int, bnd::WaveletBoundary) where {T}
    # TODO implement for different lengths.
    @assert d ≈ log2(length(coeffs))
    input = zeros(T,1<<d)
    copyto!(input, coeffs)
    T(2)^(d//2)*idwt(input, side, w, bnd)
end


function inv_cascade_algorithm_linear_combination(side::Side, w::DiscreteWavelet{T}, coeffs, d::Int, bnd::WaveletBoundary) where {T}
    # TODO implement for different lengths.
    @assert d ≈ log2(length(coeffs))
    input = zeros(T,1<<d)
    copyto!(input, coeffs)
    T(2)^(d//2)*dwt(input, side, w, bnd)
end
