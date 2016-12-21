# dwttransform.jl
using ..Util: isdyadic, maxtransformlevels

import ..Transforms: dwt, idwt
export dwt, idwt
"Perform l steps of the wavelet transform"
function dwt!
end

"Perform l steps of the inverse wavelet transform"
function idwt!
end

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
  y = zeros(T, length(x))
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
  y = zeros(T, length(x))
  idwt!(y, x, fb, bnd, L)
  y
end
