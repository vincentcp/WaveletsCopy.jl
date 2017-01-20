# periodize.jl
in_support{T}(x, support::Tuple{T,T}) = support[1] <= x < support[2]

function in_periodic_support{T}(x, support::Tuple{T,T}...)
  for s in support
    (s[1] <= x <= s[2]) && (return true)
  end
  false
end

function periodic_support{T}(side::Side, kind::Kind, w::DiscreteWavelet{T}, j, k, a=0, b=1)
  s = DWT.support(side, kind, w, j, k)
  _periodize(s, a, b)
end

function _periodize{T}(s::Tuple{T,T}, a=0, b=1)
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
function _periodize!{T}(dest::AbstractArray{T}, src::AbstractArray{T}, istart)
  L = length(dest)
  j = mod(istart, L)
  srclength = length(src)
  (j == 0) && (j = L)
  for i in 1:length(dest)
    t = T(0)
    for m in i:L:srclength
      t += src[m]
    end
    dest[j] = t
    j += 1
    (j > L) && (j = 1)
  end
end
