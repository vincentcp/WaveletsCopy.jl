# daubechies.jl

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

primal_scalingfilter{T}(w::DaubechiesWavelet{1,T}) = CompactSequence(sqrt(T(2))*db1_h, 0)

primal_scalingfilter{T}(w::DaubechiesWavelet{2,T}) =
    CompactSequence(1/sqrt(T(2))*[(1+sqrt(T(3)))/4, (3+sqrt(T(3)))/4, (3-sqrt(T(3)))/4, (1-sqrt(T(3)))/4], 0)

primal_scalingfilter(w::DaubechiesWavelet{4,Float64}) =
    CompactSequence(db4_h, 0)

db1 = DaubechiesWavelet{1,Float64}()
db2 = DaubechiesWavelet{2,Float64}()
db3 = DaubechiesWavelet{3,Float64}()
db4 = DaubechiesWavelet{4,Float64}()
db5 = DaubechiesWavelet{5,Float64}()
db6 = DaubechiesWavelet{6,Float64}()

