# Convenience method
eltype(x, y) = promote_type(eltype(x), eltype(y))
# WAVELET INDEXING AND SIZES

maxtransformlevels(x::AbstractArray) = maxtransformlevels(minimum(size(x)))
function maxtransformlevels(arraysize::Integer)
    arraysize > 1 || return 0
    tl = 0
    while sufficientpoweroftwo(arraysize, tl)
        tl += 1
    end
    return tl - 1
end

ndyadicscales(n::Integer) = round(Int, log2(n))
ndyadicscales(x::AbstractArray) = ndyadicscales(size(x,1))

function isdyadic(x::AbstractArray)
    for i = 1:ndims(x)
        isdyadic(size(x,i)) || return false
    end
    return true
end
isdyadic(n::Integer) = (n == 2^(ndyadicscales(n)))

sufficientpoweroftwo(n::Integer, L::Integer) = (n%(2^L) == 0)

"""
  A dyadic point given by `k/2^d` at most xtol separated from x.
"""
function closest_dyadic_point{T}(x::T, xtol::T; dmax = 20)
    offset = floor(x)
    x -= offset
    up = 1.;  low = 0.
    lowk = 0; upk = 1
    d = 0
    resk = 0
    while true
        mid   = (up+low)/2
        lowk  = lowk  << 1
        upk   = upk   << 1
        d += 1
        midk = Int((lowk+upk)/2)
        if x > mid
            low = mid; lowk = midk
        else
            up = mid ; upk  = midk
        end
        (abs(x-low) <= xtol)  && (resk = lowk; break)
        (abs(x-up)  <= xtol)  && (resk  = upk; break)
        (d == dmax) && (warn("dmax is reached in closest_dyadic_point"); resk = lowk)
    end
    k = (1<<d)*Int(offset) + resk
    d, k
end
