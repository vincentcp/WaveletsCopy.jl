# dwtstep.jl

using ..Sequences
using ..Filterbanks

export dwtstep, idwtstep

"""
Perform one step of the wavelet transform using the given wavelet.
"""
function dwtstep!
end

"""
Perform one step of the inverse wavelet transform using the given wavelet.
"""
function idwtstep_filter!
end


"""
Compute the size of the vectors of scaling and wavelet coefficients when applying
the given wavelet to the given vector x using the given boundary condition.
"""
function dwtstep_size
end



function dwtstep_size(n::Int, fb::Filterbank, bnd::WaveletBoundary)
    # TODO implement this
    # For odd length, we include one more scaling coefficient than we have wavelet coefficients
    # dc_len = n >> 1  # = ceil(Int, (N-1)/2)
    # sc_len = n - dc_len

    # For odd length, the number of scaling and wavelet coefficients is equal
    sc_len = (n + 1) >> 1
    dc_len = sc_len
    (sc_len, dc_len)
end

dwtstep_size(x, fb::Filterbank, bnd) = dwtstep_size(length(x), fb, bnd)

idwtstep_size(n1::Int, n2::Int, fb::Filterbank, bnd::WaveletBoundary) = n1+n2

idwtstep_size(sc, dc, fb::Filterbank, bnd::WaveletBoundary) = idwtstep_size(length(sc), length(dc), fb, bnd)


function dwtstep(x, fb::Filterbank, bnd::WaveletBoundary)
    sc_len, dc_len = dwtstep_size(x, fb, bnd)
    T = promote_type(eltype(x), eltype(fb))
    sc = zeros(T, sc_len)
    dc = zeros(T, dc_len)
    dwtstep!(sc, dc, x, fb, bnd)
    (sc, dc)
end

function idwtstep(sc, dc, fb::Filterbank, bnd::WaveletBoundary)
    x_len = idwtstep_size(sc, dc, fb, bnd)
    T = promote_type(eltype(sc), eltype(fb))
    x = zeros(T, x_len)
    idwtstep!(x, sc, dc, fb, bnd)
    x
end

## Periodic boundaries

dwtstep!(sc, dc, x, fb::Filterbank, bnd::PeriodicBoundary) =
    polyphase_analysis!(sc, dc, x, fb.pm_analysis, PeriodicEmbedding())

idwtstep!(x, sc, dc, fb::Filterbank, bnd::PeriodicBoundary) =
    polyphase_synthesis!(x, sc, dc, fb.pm_synthesis, PeriodicEmbedding())


## Symmetric boundaries
# TODO implement
#
#
# function dwtstep!(sc, dc, x, fb::Filterbank, bnd::SymmetricBoundary)
#     if iseven(length(x))
#         if isodd(sublength(dual_lowpassfilter(fb)))
#             polyphase_analysis!(sc, dc, x, fb.pm_analysis, SymmetricEmbedding{:wp,:wp,:even,:even}())
#         end
#     end
# end
#
# function idwtstep!(x, sc, dc, fb::Filterbank, bnd::SymmetricBoundary)
#     if iseven(length(x))
#         if isodd(sublength(primal_lowpassfilter(fb)))
#             polyphase_synthesis!(x, sc, dc, fb.pm_synthesis, SymmetricEmbedding{:hp,:hp,:even,:even}())
#         end
#     end
# end
#
#
# "Compute a matrix representation of one DWT step for a signal of size n."
# function dwtstep_matrix(n, fb::Filterbank, bnd::WaveletBoundary)
#     T = eltype(fb)
#     A = zeros(n,n)
#     sc_len, dc_len = dwtstep_size(n, fb, bnd)
#     sc = zeros(T, sc_len)
#     dc = zeros(T, dc_len)
#     x = zeros(T, n)
#     for i = 1:n
#         if i > 1
#             x[i-1] = 0
#         end
#         x[i] = 1
#         dwtstep!(sc, dc, x, fb, bnd)
#         A[1:sc_len,i] = sc
#         A[sc_len+1:n,i] = dc
#     end
#     A
# end
