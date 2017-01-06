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
cascade_algorithm(side::Side, kind::Kind, w::DiscreteWavelet, L=0; options...) =
    cascade_algorithm(filter(side, kind, w), L; options...)
#
cascade_algorithm(s::CompactSequence, L; options...) = cascade_algorithm(s.a, L; options...)

function cascade_algorithm{T}(h::AbstractArray{T}, L; tol = sqrt(eps(T)), options...)
  @assert sum(h)≈sqrt(T(2))
  N = length(h)
  # find ϕ(0), .. , ϕ(N-1) by solving a eigenvalue problem
  # The eigenvector with eigenvalue equal to 1/√2 is the vector containting ϕ(0), .. , ϕ(N-1)
  # see http://cnx.org/contents/0nnvPGYf@7/Computing-the-Scaling-Function for more mathematics

  # Create matrix from filter coefficients
  H = _get_H(h)
  # Find eigenvector eigv
  E = eigfact(H)
  index = find(abs(E[:values]-1/sqrt(T(2))).<tol)
  @assert length(index) > 0
  for i in index
    V = E[:vectors][:,i]
    @assert norm(imag(V)) < tol
    eigv = real(V)
    abs(sum(eigv)) > tol*100 ? break : nothing
  end
  eigv /= sum(eigv)

  # Find intermediate values ϕ(1/2), .. ,ϕ(N-1 -1/2)
  K = 2
  for m in 1:L
    eigv = vcat(eigv, zeros(T,length(eigv)-1))
    eigv_length = length(eigv)
    for n in eigv_length:-2:3
      eigv[n] = eigv[ceil(Int,n/2)]
    end
    eigv[2:2:end] = 0
    for n in 2:2:eigv_length
      I = max(1,ceil(Int,(2n-1-eigv_length)/K+1)):min(N,floor(Int, (2n-1)/K+1))
      eigv[n] = sqrt(T(2))*sum(h[i].*eigv[2*n-1-K*(i-1)] for i in I)
    end
    K = K << 1
  end
  eigv
end

dyadicpointsofcascade{T}(side::Side, kind::Kind, w::DiscreteWavelet{T}, j::Int, k::Int, d::Int) =
    T(2)^(-j)*(k+dyadicpointsofcascade(side, kind, w, d-j))
function dyadicpointsofcascade{T}(side::Side, kind::Kind, w::DiscreteWavelet{T}, L::Int)
  s = support(side, kind, w)
  H = support_length(side, kind, w)
  linspace(T(s[1]), T(s[2]), (1<<L)*H+1)
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
