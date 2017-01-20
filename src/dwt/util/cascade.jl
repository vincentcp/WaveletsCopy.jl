# cascade.jl

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

"""
The scaling function of a wavelet with filtercoeficients h evaluated in diadic points

`x[k] = 2^(-L)k for k=0..K`

where K is the length of h. Thus L=0, gives the evaluation of the
scaling function in the points [0,1,...,K], and L=1, the points [0,.5,1,...,K].
"""
function cascade_algorithm{T}(side::Side, kind::Kind, w::DiscreteWavelet{T}, L=0; options...)
  f = zeros(T, cascade_length(side, kind, w, L))
  cascade_algorithm!(f, side, kind, w, L; options...)
  f
end

function cascade_algorithm{T}(s::CompactSequence{T}, L; options...)
  f = zeros(T, cascade_length(s, L))
  cascade_algorithm!(f, s, L; options...)
  f
end

cascade_algorithm!{T}(f::AbstractArray{T,1}, side::Side, kind::Kind, w::DiscreteWavelet{T}, L=0; options...) =
    cascade_algorithm!(f::AbstractArray, filter(side, kind, w), L; options...)

cascade_algorithm!{T}(f::AbstractArray{T,1}, s::CompactSequence{T}, L; options...) = cascade_algorithm!(f::AbstractArray, s.a, L; options...)

function cascade_algorithm!{T}(f::AbstractArray{T,1}, h::AbstractArray{T,1}, L; tol = sqrt(eps(T)), options...)
  @assert L >= 0
  @assert sum(h)≈sqrt(T(2))
  N = length(h)
  # find ϕ(0), .. , ϕ(N-1) by solving a eigenvalue problem
  # The eigenvector with eigenvalue equal to 1/√2 is the vector containting ϕ(0), .. , ϕ(N-1)
  # see http://archive.cnx.org/contents/d279ef3c-661f-4e14-bba8-70a4eb5c0bcf@7/computing-the-scaling-function-the-cascade-algorithm for more mathematics

  # Create matrix from filter coefficients
  H = DWT._get_H(h)
  # Find eigenvector eigv
  E = eigfact(H)
  index = find(abs(E[:values]-1/sqrt(T(2))).<tol)
  @assert length(index) > 0
  i = index[1]
  V = E[:vectors][:,i]
  @assert norm(imag(V)) < tol
  eigv = real(V)
  (abs(sum(eigv)) < tol*100) && (warn("Cascade algorithm is not convergent"))

  eigv /= sum(eigv)
  eigv_length = N
  f[1:1<<L:end] = eigv
  # # Find intermediate values ϕ(1/2), .. ,ϕ(N-1 -1/2)
  K = 2
  for m in 1:L
    for n in 1+1<<(L-m):1<<(L-m+1):length(f)
      t = T(0)
      for i in max(1,ceil(Int, (2n-1-length(f))/(1<<L)+1)):min(N, 1+(n-1)>>(L-1))
        t += h[i].*f[2*n-1-(1<<L)*(i-1)]
      end
      f[n] = sqrt(T(2))*t
    end
  end
  nothing
end

cascade_length{T}(side::Side, kind::Kind, w::DiscreteWavelet{T}, L::Int) =
    cascade_length(support_length(side, kind, w),L)
cascade_length(f::CompactSequence, L::Int) = cascade_length(sublength(f)-1, L)
cascade_length(H::Int, L::Int) = (L >= 0) ? (1<<L)*H+1 : 1

dyadicpointsofcascade{T}(side::Side, kind::Kind, w::DiscreteWavelet{T}, j::Int, k::Int, d::Int) =
    T(2)^(-j)*(k+dyadicpointsofcascade(side, kind, w, d-j))
function dyadicpointsofcascade{T}(side::Side, kind::Kind, w::DiscreteWavelet{T}, L::Int)
  s = support(side, kind, w)
  H = support_length(side, kind, w)
  if L >= 0
    linspace(T(s[1]), T(s[2]), (1<<L)*H+1)
  else
    linspace(T(s[1]), T(s[2]), H+1)[1:1<<-L:end]
  end
end

function _get_H(h)
  N = length(h)
  T = eltype(h)
  M = zeros(T,N,N)
  for l in 1:2:N
    hh = h[l]
    for m in 1:ceil(Int,N/2)
      M[l>>1+m,1+(m-1)*2] = hh
    end
  end
  for l in 2:2:N
    hh = h[l]
    for m in 1:N>>1
      M[l>>1+m,2+(m-1)*2] = hh
    end
  end
  M
end
