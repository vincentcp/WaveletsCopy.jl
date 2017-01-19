module Util
export   maxtransformlevels,
        isdyadic
using Compat

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

ndyadicscales(n::Integer) = @compat round(Int, log2(n))
ndyadicscales(x::AbstractArray) = ndyadicscales(size(x,1))

function isdyadic(x::AbstractArray)
    for i = 1:ndims(x)
        isdyadic(size(x,i)) || return false
    end
    return true
end
isdyadic(n::Integer) = (n == 2^(ndyadicscales(n)))

sufficientpoweroftwo(n::Integer, L::Integer) = (n%(2^L) == 0)

end # module
