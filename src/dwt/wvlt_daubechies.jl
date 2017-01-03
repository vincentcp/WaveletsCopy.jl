# daubechies.jl
using ..WT: daubechies
immutable DaubechiesWavelet{N,T} <: DiscreteWavelet{T}
end

typealias HaarWavelet{T} DaubechiesWavelet{1,T}

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
primal_scalingfilter(::Type{DaubechiesWavelet{1,T0}}) = _db1_h(T0)
primal_scalingfilter(::Type{DaubechiesWavelet{2,T0}}) = _db2_h(T0)
primal_scalingfilter(::Type{DaubechiesWavelet{4,T0}}) = CompactSequence(db4_h, 0)
primal_scalingfilter{N}(::Type{DaubechiesWavelet{N,T0}}) = CompactSequence(daubechies(N), 0)

primal_coefficientfilter{T}(::Type{DaubechiesWavelet{1,T}}) = CompactSequence(db1_h)

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

primal_scaling_support{N,T}(::Type{DaubechiesWavelet{N,T}}) = (0,2N-1)
primal_vanishingmoments{N,T}(::Type{DaubechiesWavelet{N,T}}) = N

using .Cardinal_b_splines
evaluate_primal_scalingfunction{T}(w::DWT.HaarWavelet{T}, x::Number) =
      evaluate_Bspline(0, x, eltype(x, w))
