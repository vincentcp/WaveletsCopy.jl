# cascade.jl
export cascade_algorithm, dyadicpointsofcascade
"""
The scaling function of a wavelet with filtercoeficients h evaluated in diadic points

`x[k] = 2^(-L)k for k=0..K`

where K is the length of h. Thus L=0, gives the evaluation of the
scaling function in the points [0,1,...,K], and L=1, the points [0,.5,1,...,K].
"""
cascade_algorithm(w::DiscreteWavelet; dual=true, options...) = dual?
    cascade_algorithm(dual_scalingfilter(w); options...) :
    cascade_algorithm(primal_scalingfilter(w); options...)
cascade_algorithm(s::CompactSequence; options...) = cascade_algorithm(s.a; options...)
function cascade_algorithm(h; L = 0, tol = sqrt(eps(eltype(h))), options...)
  T = eltype(h)
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

dyadicpointsofcascade(s::DiscreteWavelet, L::Int=10) = dyadicpointsofcascade(primal_scalingfilter(s), L)
dyadicpointsofcascade(s::CompactSequence, L::Int=10) = dyadicpointsofcascade(s.a, L, s.offset)
dyadicpointsofcascade(h, L::Int=10, offset::Int=0) = dyadicpointsofcascade(eltype(h), length(h), L, offset)
dyadicpointsofcascade(T::Type, H::Int, L::Int, offset::Int=0) = L==0 ?
      T(offset)+linspace(T(0),T(H-1),H) : T(offset)+linspace(T(0),T(H-1),2^L*(H-1)+1)

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
