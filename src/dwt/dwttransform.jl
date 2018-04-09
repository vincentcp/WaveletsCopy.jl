# dwttransform.jl

export dwt, idwt, full_dwt, full_idwt
"Perform `L` steps of the wavelet transform"
function dwt!
end

"Perform `L` steps of the inverse wavelet transform"
function idwt!
end

dwt_size(x, w::DiscreteWavelet, bnd::WaveletBoundary) = dwt_size(x, Filterbank(w), bnd)
dwt_size(x, fb::Filterbank, bnd::WaveletBoundary) = length(x)

dwt(x, w::DiscreteWavelet, bnd::WaveletBoundary, L::Int=maxtransformlevels(x)) =
    dwt(x, Primal, w, bnd, L)

idwt(x, w::DiscreteWavelet, bnd::WaveletBoundary, L::Int=maxtransformlevels(x)) =
    idwt(x, Primal, w, bnd, L)

dwt(x, s::Prl, w::DiscreteWavelet, bnd::WaveletBoundary, L::Int=maxtransformlevels(x)) =
    dwt(x, Filterbank(w), bnd, L)

idwt(x, s::Prl, w::DiscreteWavelet, bnd::WaveletBoundary, L::Int=maxtransformlevels(x)) =
    idwt(x, Filterbank(w), bnd, L)

dwt(x, s::Dul, w::DiscreteWavelet, bnd::WaveletBoundary, L::Int=maxtransformlevels(x)) =
    dwt(x, DualFilterbank(w), bnd, L)

idwt(x, s::Dul, w::DiscreteWavelet, bnd::WaveletBoundary, L::Int=maxtransformlevels(x)) =
    idwt(x, DualFilterbank(w), bnd, L)

function dwt!(y, x, fb::Filterbank, bnd::WaveletBoundary, L::Int=maxtransformlevels(x))
    lx = length(x)
    @assert 2^L <= lx
    copy!(y,x)
    for l in 1:L
        lx2 = lx >> 1
        y[1:lx2], y[lx2+1:lx] = DWT.dwtstep(y[1:lx], fb, bnd)
        lx = lx2
    end
end

function dwt(x, fb::Filterbank, bnd::WaveletBoundary,  L::Int=maxtransformlevels(x))
    @assert isdyadic(x)
    T = promote_type(eltype(x), eltype(fb))
    y = zeros(T, dwt_size(x, fb, bnd))
    dwt!(y, x, fb, bnd, L)
    y
end

function idwt!(y, x, fb::Filterbank, bnd::WaveletBoundary, L::Int=maxtransformlevels(x))
    lx = length(x)
    @assert 2^L <= lx
    copy!(y,x)
    lx = lx >> L
    for l in L:-1:1
        lx2 = lx << 1
        t = DWT.idwtstep(y[1:lx], y[lx+1:lx2], fb, bnd)
        y[1:lx2] = t
        lx = lx2
    end
end

function idwt(x, fb::Filterbank, bnd::WaveletBoundary, L::Int=maxtransformlevels(x))
    @assert isdyadic(x)
    T = promote_type(eltype(x), eltype(fb))
    y = zeros(T, dwt_size(x, fb, bnd))
    idwt!(y, x, fb, bnd, L)
    y
end

full_dwt{T}(x, w::DiscreteWavelet{T}, bnd::WaveletBoundary) =
    full_dwt(x, Primal, w, bnd)

full_idwt{T}(x, w::DiscreteWavelet{T}, bnd::WaveletBoundary) =
    full_idwt(x, Primal, w, bnd)

function full_dwt{T}(x, s::Side, w::DiscreteWavelet{T}, bnd::WaveletBoundary)
    coefs = scaling_coefficients(x, s, w, bnd)
    dwt(coefs, s, w, bnd)
end

function full_idwt{T}(x, s::Side, w::DiscreteWavelet{T}, bnd::WaveletBoundary)
    scalingcoefs = idwt(x, s, w, bnd)
    scaling_coefficients_to_dyadic_grid(scalingcoefs, s, w, bnd)
end
