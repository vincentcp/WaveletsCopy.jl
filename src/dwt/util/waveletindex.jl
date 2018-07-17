const WaveletIndex = Tuple{Kind,Int,Int}
kind(idx::WaveletIndex) = idx[1]
level(idx::WaveletIndex) = idx[2]
offset(idx::WaveletIndex) = idx[3]
value(wind::WaveletIndex) = coefficient_index(wind)

Base.show(io::IO, wind::WaveletIndex) =
	print(io, "$(name(kind(wind))) index with level $(level(wind)), and offset $(offset(wind)) ($(value(wind)))")
# $(kind(wind)) $(offset(wind))   $(level(wind)) $(value(wind))

"""
`DWTIndexList` defines the map from native indices to linear indices
for a finite wavelet basis, when the indices are ordered in the way they
are expected in the DWT routine.
"""
struct DWTIndexList <: AbstractVector{WaveletIndex}
	l	::	Int
end

# Assume linear indexing,
Base.IndexStyle(list::DWTIndexList) = Base.IndexLinear()

Base.length(list::DWTIndexList) = 1<<list.l
Base.size(list::DWTIndexList) = (1<<list.l,)

Base.getindex(m::DWTIndexList, idx::Int) ::WaveletIndex = wavelet_index(m.l, idx)

Base.getindex(list::DWTIndexList, wind::WaveletIndex)::Int = value(wind)

" All wavelet indices on a certain level"
wavelet_indices(l::Int) = DWTIndexList(l)

"""
  The index ([:scaling/:wavelet], j, k) in the (scaling+wavelet) sequence for coefficient i in a sequence of length n after l dwt synthesis_lowpassfilter

  For example, the indices of a sequence with 4 elements after
  0 dwt steps
    (Scl(), 2, 0),    (Scl(), 2, 1),    (Scl(), 2, 2),     (Scl(), 2, 3)
  1 dwt step
    (Scl(), 1, 0),    (Scl(), 1, 1),    (Wvl(), 1, 0),     (Wvl(), 1, 1)
  2 dwt steps
    (Scl(), 0, 0),    (Wvl(), 0, 0),    (Wvl(), 1, 0),     (Wvl(), 1, 1)
"""
function wavelet_index(n::Int, i::Int, l::Int)
    if i > n/(1<<l)
        j = level(n,i)
        k = mod(i-1,1<<j)
        Wvl(), j, k
    else
        Scl(), Int(log2(n))-l, i-1
    end
end

wavelet_index(level::Int, index::Int) = wavelet_index(1<<level, index, level)

"""
  The coefficient corresponding to the wavelet index (kind, j, k). The integer
  is the index in the array of the DWT.
"""
coefficient_index(index::WaveletIndex) = coefficient_index(kind(index), level(index), offset(index))
coefficient_index(kind::Scl, j::Int, k::Int)::Int = k+1
coefficient_index(kind::Wvl, j::Int, k::Int)::Int = 1<<j+k+1

"Return the wavelet level of a coefficient in an array of length `n` with index `i`"
function level(n::Int, i::Int)
    (i == 1 || i == 2) && (return 0)
    for l in 1:round(Int,log2(n))
        if i <= (1<<(l+1))
            return l
        end
    end
end

"""
`DWaveletsCopycalingIndexList` defines the map from native indices to linear indices
for a finite wavelet basis, when the indices are ordered in the way they
are expected in the DWT routine.
"""
struct DWaveletsCopycalingIndexList <: AbstractVector{WaveletIndex}
	l	::	Int
end

# Assume linear indexing,
Base.IndexStyle(list::DWaveletsCopycalingIndexList) = Base.IndexLinear()

Base.length(list::DWaveletsCopycalingIndexList) = 1<<list.l
Base.size(list::DWaveletsCopycalingIndexList) = (1<<list.l,)

Base.getindex(m::DWaveletsCopycalingIndexList, idx::Int) ::WaveletIndex = scaling_index(m.l, idx)

Base.getindex(list::DWaveletsCopycalingIndexList, wind::WaveletIndex)::Int = scaling_value(wind)

" All wavelet indices on a certain level"
scaling_indices(l::Int) = DWaveletsCopycalingIndexList(l)

scaling_index(level::Int, index::Int) = (Scl(), level, index-1)

scaling_value(wind::WaveletIndex)::Int = offset(wind)+1
