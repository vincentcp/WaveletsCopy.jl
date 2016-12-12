# compactsequence.jl


#########################################
# CompactSequence
#########################################

"""
A CompactSequence is a compactly supported sequence of a certain length starting at a
given offset index.
"""
immutable CompactSequence{T} <: Sequence
    a       ::  Vector{T}
    offset  ::  Int
    n       ::  Int

    CompactSequence(a, offset) = new(a, offset, length(a))
end

CompactSequence{T}(a::Vector{T}, offset = 0) = CompactSequence{T}(a, offset)

CompactSequence{T}(a::AbstractVector{T}, offset = 0) = CompactSequence{T}(collect(a), offset)

eltype{T}(::Type{CompactSequence{T}}) = T

shift(s::CompactSequence, k::Int) = CompactSequence(s.a, s.offset+k)

sublength(s::CompactSequence) = s.n

mapindex(s::CompactSequence, k) = k - s.offset + 1

imapindex(s::CompactSequence, l) = l + s.offset - 1


# We override getindex to return zero outside our embedded vector.
# TODO: test whether there is overhead in calling length(s.a), rather than having the s.n field.
getindex(s::CompactSequence, k::Int) =
    (k < s.offset) || (k >= s.offset+s.n) ? zero(eltype(s)) : getindex(s.a, k-s.offset+1)
function getindex(s::CompactSequence, c::Range)
  e = Array(eltype(s), length(c))
end

function Base.getindex(s::CompactSequence, c::UnitRange)
  e = zeros(eltype(s), length(c))
  i1 = max(s.offset, c[1])
  i2 = min(s.offset + s.n-1, c[end])
  e1 = max(1, 1-c[1]+s.offset)
  e2 = min(length(e),i2-i1+1-c[1]+s.offset)
  e[e1:e2] = s.a[i1-s.offset+1:i2-s.offset+1]
  e
end

firstindex(s::CompactSequence) = imapindex(s, 1)

lastindex(s::CompactSequence) = imapindex(s, sublength(s))

hascompactsupport{T}(::Type{CompactSequence{T}}) = True

"""
From a given filter h_i, compute a new filter satisfying the alternating flip relation,
centered around the given pivot:

`g_k = (-1)^k h_{pivot-k}`

The default pivot is 1.
"""
function alternating_flip(s::CompactSequence, pivot = 1)
    hflip = similar(s.a)
    # The first element of hflip has index pivot-lastindex(f).
    # Whether or not we need to flip its sign depend on the parity of this index:
    isodd(pivot-lastindex(s)) ? t = -1 : t = 1
    hflip[1:2:end] =  t * s.a[end:-2:1]
    hflip[2:2:end] = -t * s.a[end-1:-2:1]
    CompactSequence(hflip, pivot - lastindex(s))
end

"Compute a new 'alternating' filter satisfying `g_k = (-1)^k h_k`."
function alternating(s::CompactSequence)
    halt = similar(s.a)
    t = (-1)^firstindex(s)
    for k in eachindex(s.a)
        halt[k] = t * s.a[k]
        t = -t
    end
    CompactSequence(halt, firstindex(s))
end


"The first even number greater than or equal to n."
nexteven(n) = isodd(n) ? n+1 : n

"The last even number, smaller than or equal to n."
previouseven(n) = isodd(n) ? n-1 : n

"The first odd number greater than or equal to n."
nextodd(n) = isodd(n) ? n : n+1

"Return the even part of a sequence `s`, defined by `s_e[k] = s[2k]`."
evenpart(s::CompactSequence) =
    CompactSequence(eltype(s)[s[j] for j in nexteven(firstindex(s)):2:lastindex(s)], div(nexteven(firstindex(s)),2))

"Return the odd part of a sequence `s`, defined by `s_o[k] = s[2k+1]`."
oddpart(s::CompactSequence) =
    CompactSequence(eltype(s)[s[j] for j in nextodd(firstindex(s)):2:lastindex(s)], div(previouseven(firstindex(s)),2))


ctranspose(s::CompactSequence) = CompactSequence(conj(flipdim(s.a,1)), -lastindex(s))

reverse(s::CompactSequence) = CompactSequence(flipdim(s.a, 1), -lastindex(s))

conj(s::CompactSequence) = CompactSequence(conj(s.a), firstindex(s))





#########################################
# FixedSequence
#########################################


"""
A FixedSequence has a fixed length `L` starting at a certain offset `OFS`.
It contains its values in a fixed size array.
"""
immutable FixedSequence{L,OFS,T} <: Sequence
    a       ::  Vec{L,T}
end

FixedSequence{L,T,OFS}(a::Vec{L,T}, ::Type{Val{OFS}}) = FixedSequence{L,OFS,T}(a)

# type unstable constructor
FixedSequence(a::AbstractVector, offset = 0) = FixedSequence(Vec(a...), Val{offset})

# type unstable constructor
FixedSequence(s::CompactSequence) = FixedSequence(s.a, s.offset)

convert{L,OFS,T}(::Type{CompactSequence{T}}, s::FixedSequence{L,OFS,T}) = CompactSequence(Array(s.a), OFS)

eltype{L,OFS,T}(::Type{FixedSequence{L,OFS,T}}) = T

sublength{L,OFS,T}(s::FixedSequence{L,OFS,T}) = L

getindex{L,OFS,T}(s::FixedSequence{L,OFS,T}, i) = s.a[i-OFS+1]

firstindex{L,OFS,T}(s::FixedSequence{L,OFS,T}) = OFS

lastindex{L,OFS,T}(s::FixedSequence{L,OFS,T}) = OFS+L-1

shift{L,OFS,T}(s::FixedSequence{L,OFS,T}, k::Int) = FixedSequence(s.a, Val{OFS+k})


hascompactsupport{L,OFS,T}(::Type{FixedSequence{L,OFS,T}}) = True


for op in (:ctranspose, :evenpart, :oddpart, :alternating_flip, :reverse, :conj, :alternating)
    @eval $op{L,OFS,T}(s::FixedSequence{L,OFS,T}) = FixedSequence($op(CompactSequence{T}(s)))
end
