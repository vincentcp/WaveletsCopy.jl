# derivedsequences.jl


######################################
# Type hierarchy and general interface
######################################

"""
A DerivedSequence is any sequence that is derived from another sequence.
Examples include a DownsampledSequence and an UpsampledSequence.
DerivedSequences are lazy. No computation is performed until the sequence is
evaluated. This can be done for compactly supported sequences using `collect`.
"""
abstract type DerivedSequence{S} <: Sequence end

eltype{S}(::Type{DerivedSequence{S}}) = eltype(S)
eltype{DS <: DerivedSequence}(::Type{DS}) = eltype(super(DS))

# If the subtypes only transform indices, then they can implement mapindex and imapindex,
# where mapindex maps an index 'k' of the derived sequence into an index 'l' of the original sequence,
# and imapindex is the inverse.
# We always use the symbol 'l' for an index of the original sequence and 'k' for the derived sequence.

# Default indexing rules are:
getindex{M}(s::DerivedSequence{M}, k) = s.seq[mapindex(s, k)]

setindex!{M}(s::DerivedSequence{M}, val, k) = s.seq[mapindex(s, k)] = val

"Return the original sequence of the derived sequence."
sequence(s::DerivedSequence) = s.seq

# Traits
hascompactsupport{S}(::Type{DerivedSequence{S}}) = hascompactsupport(S)


collect(s::DerivedSequence) = _collect(s, hascompactsupport(s))

_collect(s::DerivedSequence, compactsupport::True) = CompactSequence(eltype(s)[s[i] for i in fistindex(s):lastindex(s)], firstindex(s))

_collect(s::DerivedSequence, compactsupport::False) = throw(BoundsError())



######################################
# Definitions of derived sequences
######################################


"""
A DownsampledSequence is determined by a sequence `s`, a downsampling factor `M`
and a `shift`. It is defined by `ds_k = s_{shift+M*k}`.

A DownsampledSequence acts as a mutating view. One can set elements of the
sequence, and the corresponding entry of 's' will be modified.
"""
struct DownsampledSequence{M,S} <: DerivedSequence{S}
    seq     ::  S
    shift   ::  Int

    DownsampledSequence{M,S}(seq::Sequence, shift) where {M,S} = new(seq, shift)
end

# Default downsampling factor is 2.
DownsampledSequence(s) = DownsampledSequence(s, Val{2})

# This one is convenient but not type-stable.
DownsampledSequence(s, m::Int, shift) = DownsampledSequence(s, Val{m}, shift)

# Default shift is 0.
DownsampledSequence{M,S}(s::S, ::Type{Val{M}}, shift = 0) = DownsampledSequence{M,S}(s, shift)

downsample(s::Sequence, args...) = DownsampledSequence(s, args...)


mapindex{M}(s::DownsampledSequence{M}, k) = s.shift+M*k

imapindex{M}(s::DownsampledSequence{M}, l) = mod(l-s.shift, M) == 0 ? div(l-s.shift, M) : throw(BoundsError())



"""
An UpsampledSequence is determined by a sequence `s`, an upsampling factor `M`
and a `shift`. It is defined by `us_k = s_l`, for `k = shift+M*l`, and `us_k = 0` otherwise.

An UpsampledSequence acts as a mutating view. One can set elements of the
sequence, and the corresponding entry of 's' will be modified.
"""
struct UpsampledSequence{M,S} <: DerivedSequence{S}
    seq     ::  S
    shift   ::  Int

    UpsampledSequence{M,S}(seq::Sequence, shift) where {M,S} = new(seq, shift)
end

# Default upsampling factor is 2.
UpsampledSequence(s) = UpsampledSequence(s, Val{2})

# Default shift is 0.
UpsampledSequence{M,S}(s::S, ::Type{Val{M}}, shift = 0) = UpsampledSequence{M,S}(s, shift)

# This one is convenient but not type-stable.
UpsampledSequence(s, m::Int, shift) = UpsampledSequence(s, Val{m}, shift)

upsample(s::Sequence, args...) = UpsampledSequence(s, args...)


mapindex{M}(s::UpsampledSequence{M}, k) = mod(k-s.shift, M) == 0 ? k : throw(BoundsError())

getindex{M}(s::UpsampledSequence{M}, k) = mod(k-s.shift, M) == 0 ? s.seq[k] : zero(eltype(s))

each_nonzero_index{M}(s::UpsampledSequence{M}) =
    s.shift+imapindex(s, firstindex(s.seq)):M:s.shift+imapindex(s, lastindex(s.seq))


"""
A ReversedSequence is the time-reversal of a given sequence `s`, i.e. `rs_k = s_{-k}`.

A ReversedSequence acts as a mutating view. One can set elements of the
sequence, and the corresponding entry of 's' will be modified.
"""
struct ReversedSequence{S} <: DerivedSequence{S}
    seq     ::  S
end

mapindex(s::ReversedSequence, k) = -k

imapindex(s::ReversedSequence, l) = -l

firstindex(s::ReversedSequence) = -lastindex(sequence(s))

lastindex(s::ReversedSequence) = -firstindex(sequence(s))

"Reverse the given sequence."
reverse(s::Sequence) = ReversedSequence(s)

reverse(s::ReversedSequence) = sequence(s)




"""
A ShiftedSequence is determined by a `shift` and a given sequence `s`. It is
defined by `ss_k = s_{k+shift}`.

A ShiftedSequence acts as a mutating view. One can set elements of the
sequence, and the corresponding entry of 's' will be modified.
"""
struct ShiftedSequence{S} <: DerivedSequence{S}
    seq     ::  S
    shift   ::  Int
end

mapindex(s::ShiftedSequence, k) = k + s.shift

impaindex(s::ShiftedSequence, l) = l - s.shift

"Shift a sequence forward by `shift` positions."
shift(s::Sequence, shift::Int) = ShiftedSequence(s, shift)

shift(s::ShiftedSequence, shift::Int) = ShiftedSequence(sequence(s), shift+s.shift)

firstindex(s::ShiftedSequence) = firstindex(sequence(s)) + s.shift

lastindex(s::ShiftedSequence) = lastindex(sequence(s)) + s.shift
