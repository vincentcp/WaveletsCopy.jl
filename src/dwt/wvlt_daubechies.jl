# daubechies.jl
struct DaubechiesWavelet{N,T} <: DiscreteWavelet{T}
end

HaarWavelet{T} = DaubechiesWavelet{1,T}

is_symmetric{T}(::Type{DaubechiesWavelet{1,T}}) = True

is_orthogonal{N,T}(::Type{DaubechiesWavelet{N,T}}) = True
is_biorthogonal{N,T}(::Type{DaubechiesWavelet{N,T}}) = True
is_semiorthogonal{N,T}(::Type{DaubechiesWavelet{N,T}}) = True

db1_h = [1, 1]

db4_h = [0.23037781330889648, 0.7148465705529157, 0.6308807679298589, -0.027983769416860003,
    -0.1870348117190931, 0.030841381835560722, 0.03288301166688518, -0.010597401785069035]
_db1_h(T::Type) = CompactSequence(T(1)/sqrt(T(2))*db1_h, 0)
_db2_h(T::Type) = CompactSequence(1/sqrt(T(2))*[(1+sqrt(T(3)))/4, (3+sqrt(T(3)))/4, (3-sqrt(T(3)))/4, (1-sqrt(T(3)))/4], 0)

T0 = Float64
filter(::Prl, ::Scl, ::Type{DaubechiesWavelet{1,T0}}) = _db1_h(T0)
filter(::Prl, ::Scl, ::Type{DaubechiesWavelet{2,T0}}) = _db2_h(T0)
filter(::Prl, ::Scl, ::Type{DaubechiesWavelet{4,T0}}) = CompactSequence(db4_h, 0)
filter{N}(::Prl, ::Scl, ::Type{DaubechiesWavelet{N,T0}}) = CompactSequence(daubechies(N), 0)
filter{T}(::Prl, ::Cof, ::Type{DaubechiesWavelet{1,T}}) = CompactSequence(db1_h)

IMPLEMENTED_DB_WAVELETS = []
for N in 1:10
  db = Symbol(string("db",N))
  T = DaubechiesWavelet{N,T0}
  @eval begin
    $db = $T()
    class(::$T) = string($T)
    push!(IMPLEMENTED_DB_WAVELETS,$db)
  end
end
name{N,T}(::Type{DaubechiesWavelet{N,T}}) = string("db",N,"_",T)
name{N}(::Type{DaubechiesWavelet{N,Float64}}) = string("db",N)

support{N,T}(::Prl, ::Scl, ::Type{DaubechiesWavelet{N,T}}) = (0,2N-1)
vanishingmoments{N,T}(::Prl, ::Type{DaubechiesWavelet{N,T}}) = N

using .Cardinal_b_splines
# Commented for testing purposes
# evaluate{T,S<:Real}(::Prl, ::Scl, w::DWT.HaarWavelet{T}, j, k, x::Number; options...) =
#       evaluate(Prl(), Scl(), CDFWavelet{1,1,T}(), j, k, x; options...)

##########################################################################
# From original Wavelets package
##########################################################################
function daubechies(N::Int)
    @assert N > 0
    # Create polynomial
    C = Array{Int}(N)
    @inbounds for n = 0:N-1
        C[N-n] = binomial(N-1+n, n)
    end

    # Find roots in y domain (truncated binomial series; (1 - y)^{-N})
    Y = roots(C)

    # Find roots in z domain:
    # z + z^{-1} = 2 - 4*y
    # where y is a root from above
    Z = zeros(Complex128, 2*N-2)
    @inbounds for i = 1:N-1
        Yi = Y[i]
        d = 2*sqrt( Yi*Yi - Yi )
        y2 = 1 - 2*Yi
        Z[i] = y2 + d
        Z[i+N-1] = y2 -d
    end

    # Retain roots inside unit circle
    nr = 0  # count roots
    @inbounds for i = eachindex(Z)
        if abs(Z[i]) <= 1 + eps()
            nr += 1
        end
    end

    # Find coefficients of the polynomial
    # (1 + z)^N * \prod_i (z - z_i)
    R = Array{Complex128}(N+nr)
    @inbounds for i = 1:N
        R[i] = -1
    end
    k = N
    @inbounds for i = eachindex(Z)
        if abs(Z[i]) <= 1 + eps()
            k += 1
            R[k] = Z[i]
        end
    end
    HH = vieta( R )

    # Normalize coefficients
    scale!(HH, 1/norm(HH))
    return real(HH)
end

# Compute roots of polynomial
# Input is a coefficient vector with highest powers first
function roots(C::AbstractVector)
    A = compan(C)
    return eigvals(A)
end

# Create companion matrix for a polynomial
# Input is a coefficient vector with highest powers first
function compan(C::AbstractVector)
    n = length(C)
    A = zeros(n-1, n-1)

    if n > 1
        @inbounds A[1,:] = -C[2:end] ./ C[1]
        @inbounds A[2:n:end] = 1
    end
    return A
end

# Vieta-like formula for computing polynomial coefficients from roots
# See
# http://www.mathworks.se/help/matlab/ref/poly.html
function vieta(R::AbstractVector)
    n = length( R )
    C = zeros(Complex128, n+1)
    C[1] = 1
    Ci::Complex128 = 0
    Cig::Complex128 = 0

    @inbounds for k = 1:n
        Ci = C[1]
        for i = 1:k
            Cig = C[i+1]
            C[i+1] = Cig - R[k] * Ci
            Ci = Cig
        end
    end
    return C
end
