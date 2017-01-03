# sequence.jl

module Sequences

using FixedSizeArrays

import Base: eltype, getindex, setindex!, eachindex, collect

import Base: &, |, *, transpose, ctranspose, conj, sum, +, -, /

import Base: convert, widen

import Base: reverse

import ..Util: upsample, downsample

# Main abstract types
export Sequence, ExtensionSequence, DerivedSequence

# Utility function
export promote_eltype
# Traits
export True, False, hascompactsupport

export moment, ztransform, fouriertransform, support

export evenpart, oddpart, alternating_flip, alternating

export firstindex, lastindex, each_nonzero_index

# Traits
export True, False, hascompactsupport

# Extension sequences
export PeriodicExtension, ZeroPadding, ConstantPadding, ShiftedExtension,
    SymmetricExtension, UndefinedExtension

export each_subindex, subvector, sublength

# Compactly supported sequences
export CompactSequence, FixedSequence

# Derived sequences
export UpsampledSequence, DownsampledSequence, ReversedSequence, ShiftedSequence

# Embedding sequences
export EmbeddingSequence, PeriodicEmbedding, SymmetricEmbedding, FunctionEmbedding

export shift, reverse, upsample, downsample


typealias True Val{true}
typealias False Val{false}

(&){T1,T2}(::Type{Val{T1}}, ::Type{Val{T2}}) = Val{T1 & T2}
(|){T1,T2}(::Type{Val{T1}}, ::Type{Val{T2}}) = Val{T1 | T2}

(&){T1,T2}(::Val{T1}, ::Val{T2}) = Val{T1 & T2}()
(|){T1,T2}(::Val{T1}, ::Val{T2}) = Val{T1 | T2}()



"""
Any subtype of Sequence is a bi-infinite list of values, with one value for each integer in Z.

Each Sequence has an element type, given by `eltype`.
"""
abstract Sequence


"For a compactly supported sequence, `firstindex` returns the index of the first non-zero element."
firstindex(s::Sequence) = -inf

"For a compactly supported sequence, `lastindex` returns the index of the last non-zero element."
lastindex(s::Sequence) = inf

# By convention, eachindex sums only over nonzero elements of the sequence.
eachindex(s::Sequence) = firstindex(s):lastindex(s)

# Range of indices that covers all nonzero elements of both sequences.
eachindex(s1::Sequence, s2::Sequence) =
    min(firstindex(s1),firstindex(s2)):max(lastindex(s1),lastindex(s2))

"An iterator over all non-zero values in the sequence."
each_nonzero_index(s::Sequence) = eachindex(s)

#getindex(s::Sequence, r::Range) = eltype(s)[s[k] for k in range]

hascompactsupport{S <: Sequence}(::Type{S}) = False
hascompactsupport(s::Sequence) = hascompactsupport(typeof(s))()


transpose(s::Sequence) = reverse(s)

"The j-th discrete moment of a sequence is defined as `\sum_k h_k k^j`."
function moment(s::Sequence, j)
    z = zero(eltype(s))
    for k in each_nonzero_index(s)
        z += s[k] * k^j
    end
    z
end


"""
The Z transform of a sequence is a continuous function of `z`, defined by
`S(z) = \sum_k s_k z^{-k}`.
"""
function ztransform(s::Sequence, z)
    T = promote_type(eltype(s), eltype(z), eltype(1/z))
    S = zero(T)
    for k in each_nonzero_index(s)
        S += s[k] * z^(-k)
    end
    S
end

"""
The Fourier transform of a sequence is defined by `S(ω) = \sum_k s_k e^{-i ω k}`. It is like
the Z transform with `z = e^{i ω}`. The Fourier transform is a `2π`-periodic continuous function
of `ω`.
"""
fouriertransform(s::Sequence, ω) = ztransform(s, exp(im*ω))


"ZTransform is a wrapper type for the `ztransform` function."
immutable ZTransform{S}
    seq     ::  S
end

(f::ZTransform)(z) = ztransform(z.seq, z)

"FourierTransform is a wrapper type for the `fouriertransform` function."
immutable FourierTransform{S}
    seq     ::  S
end

(f::FourierTransform)(ω) = fouriertransform(z.seq, ω)



immutable Convolution{S1,S2} <: Sequence
    s1  ::  S1
    s2  ::  S2
end

(*)(s1::Sequence, s2::Sequence) = Convolution(s1, s2)

ztransform(s::Convolution, z) = ztransform(z.s1, z) * ztransform(z.s2, z)

function getindex(s::Convolution, k)
    #TODO: implement
end


include("compactsequences.jl")
include("extensionsequences.jl")
include("derivedsequences.jl")
include("embeddingsequences.jl")


promote_eltype{T}(s::CompactSequence{T}, ::Type{T}) = s

end # module
