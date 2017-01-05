# cascade.jl
export cascade_algorithm, dyadicpointsofcascade

function eval_periodic_wavelet{T, S<:Real}(side::Symbol, kind::Symbol, w::DiscreteWavelet{T}, j::Int, k::Int, x::S, xtol::S; options...)
  offset = floor(x)
  x -= offset
  d, kpoint = closest_dyadic_point(x, xtol; options...)
  (!in_support(kpoint/2^d, periodic_support(side, kind, w, j, k)...)) && return T(0)
  f = periodic_function_in_dyadic_points(side, kind, w, j, k, d)
  f[kpoint+1]
end

function in_support{T}(x, support::Tuple{T,T}...)
  for s in support
    (s[1] <= x <= s[2]) && (return true)
  end
  false
end

function periodic_support{T}(side::Symbol, kind::Symbol, w::DiscreteWavelet{T}, j, k)
  if side == :primal
    if kind == :scaling
      s =  primal_scalingsupport(w, j, k)
    elseif kind == :wavelet
      s =  primal_waveletsupport(w, j, k)
    end
  elseif side == :dual
    if kind == :scaling
      s = dual_scalingsupport(w, j, k)
    elseif kind == :wavelet
      s = dual_waveletsupport(w, j, k)
    end
  end
  _periodize(s)
end

function _periodize{T}(s::Tuple{T,T}, a=T(0), b=T(1))
  a = T(a)
  b = T(b)
  p = b-a
  ((s[2]-s[1]) >= p) && (return ((a,b),))
  offset = -p*fld(s[1]-a, p)
  s = (s[1]+offset, s[2]+offset)
  (s[2] <= b) && (return (s,))
  (s[2] > b) && (return (s[1], b), (a, s[2]-(b-a)))
end

function periodic_function_in_dyadic_points(side::Symbol, kind::Symbol, w::DiscreteWavelet, j::Int, k::Int, d::Int; options...)
  if side == :primal
    if kind == :scaling
      return periodic_primal_scalingfunction_in_dyadic_points(w, j, k, d; options...)
    elseif kind == :wavelet
      return periodic_primal_waveletfunction_in_dyadic_points(w, j, k, d; options...)
    end
  elseif side == :dual
    if kind == :scaling
      return periodic_dual_scalingfunction_in_dyadic_points(w, j, k, d; options...)
    elseif kind == :wavelet
      return periodic_dual_waveletfunction_in_dyadic_points(w, j, k, d; options...)
    end
  end
end

"""
  A dyadic point given by `k/2^d` at most xtol separated from x.
"""
function closest_dyadic_point(x, xtol; dmax = 20)
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

function _periodize{T}(n::Int, src::AbstractArray{T}, istart)
  dest = zeros(T,n)
  _periodize!(dest, src, istart)
  dest
end

function _periodize!{T}(dest::AbstractArray{T}, src::AbstractArray{T}, istart)
  L = length(dest)
  j = mod(istart, L)
  (j == 0) && (j = L)
  for i in 1:length(dest)
    dest[j] = sum( src[i:L:end] )
    j += 1
    (j > L) && (j = 1)
  end

end
for (periodicfunction, simplefunction, support) in
    ((:periodic_primal_waveletfunction_in_dyadic_points,  :primal_waveletfunction_in_dyadic_points,   :primal_waveletsupport),
     (:periodic_dual_waveletfunction_in_dyadic_points,    :dual_waveletfunction_in_dyadic_points,     :dual_waveletsupport),
     (:periodic_primal_scalingfunction_in_dyadic_points,  :primal_scalingfunction_in_dyadic_points,   :primal_scalingsupport),
     (:periodic_dual_scalingfunction_in_dyadic_points,    :dual_scalingfunction_in_dyadic_points,     :dual_scalingsupport))
  @eval begin
    function $periodicfunction{T}(w::DiscreteWavelet{T}, j=0, k=0, d=10; a = T(0), b = T(1), points=false, options...)
      s = $simplefunction(w,j,k,d)
      L = round(Int, (1<<d)*(b-a))
      @assert abs(L-(1<<d)*(b-a)) < eps(real(T))
      f = _periodize(L, s, Int((1<<d)*$support(w,j,k)[1])+1)
      if points
        f, linspace(a,b,L+1)[1:end-1]
      else
        f
      end
    end
  end
end
for (waveletfunction, scalingfunction, waveletfilter, scalingfilter, waveletsupportlength, waveletpoints, scalingpoints) in
    ((:primal_waveletfunction_in_dyadic_points, :primal_scalingfunction_in_dyadic_points, :primal_waveletfilter,  :primal_scalingfilter,  :primal_waveletsupport_length,  :primal_wavelet_dyadicpointsofcascade,  :primal_scaling_dyadicpointsofcascade),
     (:dual_waveletfunction_in_dyadic_points,   :dual_scalingfunction_in_dyadic_points,   :dual_waveletfilter,    :dual_scalingfilter,    :dual_waveletsupport_length,    :dual_wavelet_dyadicpointsofcascade,    :dual_scaling_dyadicpointsofcascade))
  @eval begin
    function $scalingfunction{T}(w::DiscreteWavelet{T}, j=0, k=0, d=10; points=false, options...)
      @assert (d-j) >= 0
      f = scalingfunction_in_dyadic_points($scalingfilter(w), j, k, d; options...)
      if points
        f, $scalingpoints(w,j,k,d; options...)
      else
        f
      end
    end

    function $waveletfunction{T}(w::DiscreteWavelet{T}, j=0, k=0, d=10; points=false, options...)
      s = $scalingfunction(w, j+1, k, d; options...)
      filter = DWT.$waveletfilter(w)
      L = DWT.$waveletsupportlength(w)
      f = zeros(T, (1<<(d-j))*L+1)
      for (i,l) in enumerate(firstindex(filter):lastindex(filter))
        offset = (1<<(d-1-j))*(i-1)
        f[offset+1:offset+length(s)] += filter[l]*s
      end
      if points
        f, $waveletpoints(w,j,k,d; options...)
      else
        f
      end
    end
  end
end

function scalingfunction_in_dyadic_points{T}(s::CompactSequence{T}, j=0, k=0, d=0; options...)
  T(2)^(T(j)/2)*cascade_algorithm(s, (d-j); options...)
end

# scaling_function_in_dyadic_points(DWT.db1)
"""
The scaling function of a wavelet with filtercoeficients h evaluated in diadic points

`x[k] = 2^(-L)k for k=0..K`

where K is the length of h. Thus L=0, gives the evaluation of the
scaling function in the points [0,1,...,K], and L=1, the points [0,.5,1,...,K].
"""
cascade_algorithm(w::DiscreteWavelet, L=0; dual=true, options...) = dual?
    cascade_algorithm(dual_scalingfilter(w), L; options...) :
    cascade_algorithm(primal_scalingfilter(w), L; options...)
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

# dyadicpointsofcascade{T}(w::DiscreteWavelet{T}, j::Int, k::Int, d::Int; options...) =
#   T(2)^(-j)*(k+dyadicpointsofcascade(w, d-j; options...))
#
# dyadicpointsofcascade{T}(w::DiscreteWavelet{T}, L::Int=10; dual=true, options...) = dual ?
#     _dyadicpointsofcascade_dual(w, L; options...) :
#     _dyadicpointsofcascade_primal(w, L; options...)


# create methods
# [primal,dual]_[scaling,wavelet]_dyadicpointsofcascade{T}(w::DiscreteWavelet{T}, j::Int, k::Int, d::Int; options...)
# [primal,dual]_[scaling,wavelet]_dyadicpointsofcascade{T}(w::DiscreteWavelet{T}, L::Int; options...)
for kind in (:dual, :primal)
  for f in (:scaling, :wavelet)
    start = string(kind,"_",f)
    f = Symbol(string(start,"_dyadicpointsofcascade"))
    support = Symbol(string(start,"support"))
    support_length = Symbol(string(start,"support_length"))
    @eval begin
      $f{T}(w::DiscreteWavelet{T}, j::Int, k::Int, d::Int; options...) =
          T(2)^(-j)*(k+$f(w, d-j; options...))
      function $f{T}(w::DiscreteWavelet{T}, L::Int; options...)
        s = $support(w)
        H = $support_length(w)
        linspace(T(s[1]), T(s[2]), (1<<L)*H+1)
      end
    end
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
