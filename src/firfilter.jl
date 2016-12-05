# FIRFilter.jl

module Filters

using Wavelets.Sequences

using FixedSizeArrays

import Base: eltype, length, eachindex, getindex

import Base: -, ctranspose

import Wavelets.Sequences: firstindex, lastindex

import Base: filter!

typealias FixedFilter FixedSequence

export filter!
export FixedFilter

"The smallest index for y that is unaffected by the boundary region."
lowerboundary(filter, n) = max(0, lastindex(filter))

"The largest index for y that is unaffected by the boundary region."
upperboundary(filter, n) = min(n-1, n-1+firstindex(filter))

"Return a temporary array with elements of the given type and of length L. Memory is allocated only
on the first call of this function and reused afterwards."
@generated temparray{T,L}(::Type{T}, ::Type{Val{L}}) = quote
    a = $(Array(T,L-1))
    return a
end


function filter!{L,OFS,T}(y, x, f::FixedFilter{L,OFS,T}, ext_x, y_i1, y_i2)
    x_vals = temparray(T, Val{L})

    j = y_i1 - OFS - L + 1
    for l in 1:L-1
        if j < 0
            x_vals[l] = ext_x[j]
        else
            x_vals[l] = x[j+1]
        end
        j += 1
    end

    l1 = lowerboundary(f, sublength(ext_x))
    l2 = upperboundary(f, sublength(ext_x))

    # In the boundary region: use extension for x-values
    for i in y_i1:l1-1
        xj = ext_x[j]
        y[i+1] = filter_mainloop!(xj, x_vals, f)
        j += 1
    end

    # Away from the boundary, use the vector x itself
    for i in l1:min(l2, y_i2)
        xj = x[j+1]
        y[i+1] = filter_mainloop!(xj, x_vals, f)
        j += 1
    end

    # Boundary region on the right
    for i in l2+1:y_i2
        xj = ext_x[j]
        y[i+1] = filter_mainloop!(xj, x_vals, f)
        j += 1
    end
end

@inline function filter_mainloop!{L,OFS,T}(xj, x_vals, f::FixedFilter{L,OFS,T})
    z = f[OFS]*xj
    @inbounds for k = 0:L-3
        z += x_vals[k+1] * f[OFS+L-1-k]
        x_vals[k+1] = x_vals[k+2]
    end
    z += x_vals[L-1] * f[OFS+1]
    x_vals[L-1] = xj
    z
end



"""
A FIRFilter is a finite impulse response filter, i.e. a linear and time-invariant filter with
a finite impulse response.

The FIR filter is respresented by its impulse response `h`, which is non-zero only starting at
a certain offset and has finite length L:

`h_i = 0 ` for `i < offset`
`h_i = 0 ` for `i >= offset + L`

A causal filter has offset 0.
"""
immutable FIRFilter{T} <: Sequences.Sequence
    h       ::  Array{T,1}
    offset  ::  Int
end

# By default the filter is causal, i.e. offset = 0.
FIRFilter{T}(filter::Array{T,1}, offset = 0) = FIRFilter{T}(filter, offset)

eltype{T}(::Type{FIRFilter{T}}) = T

length(f::FIRFilter) = length(f.h)

offset(f::FIRFilter) = f.offset

filter(f::FIRFilter) = f.h

"Shift the filter by `pos` positions to the right (forward in time)."
shift(f::FIRFilter, pos) = FIRFilter(f.h, f.offset+pos)

"Is the given filter a causal filter?"
iscausal(f::FIRFilter) = f.offset == 0

(-)(f::FIRFilter) = FIRFilter(-f.h, f.offset)

"The first non-zero index of the FIR filter."
firstindex(f::FIRFilter) = f.offset

"The last non-zero index of the FIR filter."
lastindex(f::FIRFilter) = f.offset+length(f)-1

getindex(f::FIRFilter, i) =
    (i >= firstindex(f)) && (i <= lastindex(f)) ? f.h[i-f.offset+1] : zero(eltype(f))

eachindex(f::FIRFilter) = firstindex(f):lastindex(f)


"""
From a given filter h_i, compute a new filter satisfying the alternating flip relation,
centered around the given pivot:

`g_k = (-1)^k h_{pivot-k}`

The default pivot is 1.
"""
function alternating_flip(f::FIRFilter, pivot = 1)
    hflip = similar(f.h)
    # The first element of hflip has index pivot-lastindex(f).
    # Whether or not we need to flip its sign depend on the parity of this index:
    isodd(pivot-lastindex(f)) ? s = -1 : s = 1
    hflip[1:2:end] =  s * f.h[end:-2:1]
    hflip[2:2:end] = -s * f.h[end-1:-2:1]
    FIRFilter(hflip, pivot - lastindex(f))
end

"""
From a given filter h_i, compute its hermitian conjugate as defined by:
`g_k = conj(h_{-k})`.
This is the complex conjugate of the time-reversed filter.
"""
function ctranspose(f::FIRFilter)
    FIRFilter(eltype(f)[conj(f[k]) for k in lastindex(f):-1:firstindex(f)], -lastindex(f))
end

"""
The j-th discrete moment of a FIR filter is defined as:

`\sum_k h_k k^j`
"""
moment(f::FIRFilter, j) = sum([f[k]*k^j for k in eachindex(f)])


"""
The Z transform of a FIR filter is defined as

`H(z) = \sum_k h_k z^k`
"""
ztransform(f::FIRFilter, z) = sum([f[k]*(1/z)^k for k in eachindex(f)])

"""
The Fourier transform of a FIR filter is defined as

`H(ω) = \sum_k h_k exp(i ω k)`

It is periodic with period `2π`.
"""
fouriertransform(f::FIRFilter, ω) = sum([f[k]*exp(-im*ω*k) for k in eachindex(f)])

function modulation_matrix(h::FIRFilter, g::FIRFilter, z)
    A = zeros(promote_type(eltype(h),eltype(g),eltype(z)), 2, 2)
    A[1,1] = ztransform(h, z)
    A[1,2] = ztransform(h, -z)
    A[2,1] = ztransform(g, z)
    A[2,2] = ztransform(g, -z)
    A
end

end # module
