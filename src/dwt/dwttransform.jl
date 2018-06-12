# dwttransform.jl

export dwt, idwt, full_dwt, full_idwt
"Perform `L` steps of the wavelet transform"
function dwt!
end

"Perform `L` steps of the inverse wavelet transform"
function idwt!
end

dwt_size(x, fb::Filterbank, bnd::WaveletBoundary) = length(x)

dwt(m::Matrix, w::DiscreteWavelet, bnd::WaveletBoundary, d = 2, L::Int=maxtransformlevels(size(m, 3-d))) =
    dwt(m, Primal, w, bnd, d, L)

dwt(m::Matrix{T}, s::Side, w::DiscreteWavelet{T}, bnd::WaveletBoundary, d = 2, L::Int=maxtransformlevels(size(m, 3-d))) where {T}=
    dwt!(copy(m), s, w, bnd, d, L)

function dwt!(m::Matrix{T}, s::Side, w::DiscreteWavelet{T}, bnd::WaveletBoundary, d = 2, L::Int=maxtransformlevels(size(m, 3-d))) where {T}
    if d ==2
        ld, ls = size(m)
        u = Array{T}(ls)
        v = Array{T}(ls)
        t = Array{T}(ls)
        fb = SFilterBank(s, w)
        os = 1
        @inbounds for i in 1:ld
            # @time u .= m[:,i]
            # @time dwt!(v, u, fb, bnd, L, t)
            # m[:,i] .= v
            copy!(u, 1, m, os, ls)
            dwt!(v, u, fb, bnd, L, t)
            copy!(m, os, v, 1, ls)
            os += ls
        end
    end
    m
end

function idwt!(m::Matrix{T}, s::Side, w::DiscreteWavelet{T}, bnd::WaveletBoundary, d = 2, L::Int=maxtransformlevels(size(m, 3-d))) where {T}
    if d ==2
        ld, ls = size(m)
        u = Array{T}(ls)
        v = Array{T}(ls)
        t = Array{T}(ls)
        fb = SFilterBank(s, w)
        os = 1
        @inbounds for i in 1:ld
            copy!(u, 1, m, os, ls)
            idwt!(v, u, fb, bnd, L, t)
            copy!(m, os, v, 1, ls)
            os += ls
        end
    end
    m
end

SFilterBank(s::Side, w::DiscreteWavelet) = s==Primal? Filterbank(w) : DualFilterbank(w)

dwt(x::AbstractVector, w::DiscreteWavelet, bnd::WaveletBoundary, L::Int=maxtransformlevels(x)) =
    dwt(x, Primal, w, bnd, L)

idwt(x::AbstractVector, w::DiscreteWavelet, bnd::WaveletBoundary, L::Int=maxtransformlevels(x)) =
    idwt(x, Primal, w, bnd, L)

dwt(x::AbstractVector, s::Side, w::DiscreteWavelet, bnd::WaveletBoundary, L::Int=maxtransformlevels(x)) =
    dwt(x, SFilterBank(s, w), bnd, L)

idwt(x, s::Side, w::DiscreteWavelet, bnd::WaveletBoundary, L::Int=maxtransformlevels(x)) =
    idwt(x, SFilterBank(s, w), bnd, L)

# function dwt!(y, x, fb::Filterbank, bnd::WaveletBoundary, L::Int=maxtransformlevels(x))
#     lx = length(x)
#     @assert 2^L <= lx
#     copy!(y,x)
#     @inbounds for l in 1:L
#         lx2 = lx >> 1
#         y[1:lx2], y[lx2+1:lx] = DWT.dwtstep(y[1:lx], fb, bnd)
#         lx = lx2
#     end
# end

function dwt!(y, x, fb::Filterbank, bnd::WaveletBoundary, L::Int=maxtransformlevels(x), s=similar(y))
    lx = length(x)
    @assert 2^L <= lx
    copy!(y,x)
    @inbounds for l in 1:L
        lx2 = lx >> 1
        l1, l2 = dwtstep_size(lx, fb, bnd)
        DWT.dwtstep!(s, l1, l2, y, lx, fb, bnd)
        copy!(y, 1, s, 1, lx)
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

# function idwt!(y, x, fb::Filterbank, bnd::WaveletBoundary, L::Int=maxtransformlevels(x))
#     lx = length(x)
#     @assert 2^L <= lx
#     copy!(y,x)
#     lx = lx >> L
#     @inbounds for l in L:-1:1
#         lx2 = lx << 1
#         t = DWT.idwtstep(y[1:lx], y[lx+1:lx2], fb, bnd)
#         copy!(y, 1, t, 1, lx2)
#         lx = lx2
#     end
# end

function idwt!(y, x, fb::Filterbank, bnd::WaveletBoundary, L::Int=maxtransformlevels(x),s=similar(x))
    lx = length(x)
    @assert 2^L <= lx
    copy!(y,x)
    lx = lx >> L
    @inbounds for l in L:-1:1
        lx2 = lx << 1
        # t = DWT.idwtstep(y[1:lx], y[lx+1:lx2], fb, bnd)
        DWT.idwtstep!(s, lx2, y, lx, lx, fb, bnd)
        copy!(y, 1, s, 1, lx2)
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
