# periodize.jl
in_support(x, support::Tuple{T,T}) where {T} = support[1] <= x < support[2]

function in_periodic_support(x, support::Tuple{T,T}...) where {T}
  for s in support
    (s[1] <= x <= s[2]) && (return true)
  end
  false
end

function periodic_support(side::Side, kind::Kind, w::DiscreteWavelet{T}, j, k, a=0, b=1) where {T}
  s = DWT.support(side, kind, w, j, k)
  _periodize(s, a, b)
end

function _periodize(s::Tuple{T,T}, a=0, b=1) where {T}
  a = T(a)
  b = T(b)
  p = b-a
  ((s[2]-s[1]) >= p) && (return ((a,b),))
  offset = -p*fld(s[1]-a, p)
  s = (s[1]+offset, s[2]+offset)
  (s[2] <= b) && (return (s,))
  (s[2] > b) && (return (s[1], b), (a, s[2]-(b-a)))
end

# Periodize src into an array of length n as follows (n=3)
#     src1, scr2, scr3, scr4, scr5, scr6, src7
#             ^
#          istart
#  =>      scr2, scr3, scr4
#          src5, src6, src7
#          src1
#        + ________________
#  =>      res1, res2, res3
function _periodize!(dest::AbstractArray{T}, src::AbstractArray{T}, istart::Int, step::Int=1) where {T}
    L = length(dest)
    srclength = length(src)
    for i in 1:L
        t = T(0)
        for m in mod(istart-1+step*(i-1),step*L)+1:step*L:srclength
            t += src[m]
        end
        dest[i] = t
    end
end

function _periodize_add!(dest::AbstractArray{T}, src::AbstractArray{T}, istart::Int, step::Int) where {T}
    L = length(dest)
    srclength = length(src)
    for i in 1:L
        t = T(0)
        loop_start = istart+step*(i-1)
        for m in mod(istart-1+step*(i-1),step*L)+1:step*L:srclength
            t += src[m]
        end
        dest[i] += t
    end
end

function _periodize_add!(dest::AbstractArray{T}, src::AbstractArray{T}, istart::Int) where {T}
    L = length(dest)
    srclength = length(src)
    if L <= srclength
        for i in 1:L
            t = T(0)
            for m in mod(istart-1+(i-1),L)+1:L:srclength
                t += src[m]
            end
            dest[i] += t
        end
    else
        dest_start = mod(1-istart, L)+1
        dest_end = dest_start+srclength-1
        if dest_end > L
            add!(dest, dest_start, src, 1, L-dest_start+1)
            add!(dest, 1, src, L-dest_start+2, dest_end-L)
        else
            add!(dest, dest_start, src, 1, srclength)
        end
    end
end

@inline function add!(dest::Array{T}, dest_o::Int, src::Array{T}, src_o::Int, N::Int) where T
    for i in 0:N-1
        dest[dest_o + i] += src[src_o+i]
    end
end
