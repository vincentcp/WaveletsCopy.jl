# compactsequence.jl


#########################################
# CompactSequence
#########################################

"""
A CompactSequence is a compactly supported sequence of a certain length starting at a
given offset index.
"""
struct CompactSequence{T} <: Sequence
    a       ::  Vector{T}
    offset  ::  Int
    n       ::  Int

    CompactSequence{T}(a, offset) where T = new(a, offset, length(a))
end

CompactSequence(a::Vector{T}, offset = 0) where {T} = CompactSequence{T}(a, offset)

CompactSequence(a::AbstractVector{T}, offset = 0) where {T} = CompactSequence{T}(collect(a), offset)

eltype(::Type{CompactSequence{T}}) where {T} = T

promote_eltype(s::CompactSequence{T1}, ::Type{T2}) where {T1,T2} = CompactSequence(convert(Vector{promote_type(T1,T2)}, s.a), s.offset)

shift(s::CompactSequence, k::Int) = CompactSequence(s.a, s.offset+k)

sublength(s::CompactSequence) = s.n

offset(s::CompactSequence) = s.offset

mapindex(s::CompactSequence, k) = k - s.offset + 1

imapindex(s::CompactSequence, l) = l + s.offset - 1


# We override getindex to return zero outside our embedded vector.
# TODO: test whether there is overhead in calling length(s.a), rather than having the s.n field.
getindex(s::CompactSequence, k::Int) =
    (k < s.offset) || (k >= s.offset+s.n) ? zero(eltype(s)) : getindex(s.a, k-s.offset+1)
if VERSION < v"0.7-"
    AbstractRange = Range
end
function getindex(s::CompactSequence, c::AbstractRange)
  e = Array(eltype(s), length(c))
end

function Base.getindex(s::CompactSequence, c::UnitRange)
    e = zeros(eltype(s), length(c))
    i1 = max(s.offset, c[1])
    i2 = min(s.offset + s.n-1, c[end])
    e[(1:i2-i1+1)+max(0,s.offset-c[1])] = s.a[i1-s.offset+1:i2-s.offset+1]
    e
end

firstindex(s::CompactSequence) = imapindex(s, 1)

lastindex(s::CompactSequence) = imapindex(s, sublength(s))

hascompactsupport(::Type{CompactSequence{T}}) where {T} = True

Base.isapprox(s1::CompactSequence, s2::CompactSequence) = s1.a≈s2.a

function shifted_conv(c1::CompactSequence{ELT}, c2::CompactSequence{ELT}, shift::Int) where {ELT}
    l1 = sublength(c1)
    l2 = sublength(c2)
    o1 = c1.offset
    o2 = c2.offset
    offset = o1+shift*o2
    L = (l2-1)+shift*(l1-1)+1
    a = zeros(ELT, L)
    for ai in 0:L-1
        t = ELT(0)
        # for k in max(0,floor(Int,(firstindex(c1)-l2+1)//shift)):max(0,floor(Int,(lastindex(c1)+l2)//shift))
        for k in firstindex(c1):lastindex(c1)
            t += c1.a[k-c1.offset+1]*c2[ai-shift*k]
        end
        a[ai+1] = t
    end
    CompactSequence(a, offset)
end

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

if VERSION < v"0.7-"
    ctranspose(s::CompactSequence) = CompactSequence(conj(reverse(s.a,1)), -lastindex(s))
else
    adjoint(s::CompactSequence) = CompactSequence(conj(reverse(s.a,1)), -lastindex(s))
end

reverse(s::CompactSequence) = CompactSequence(reverse(s.a, 1), -lastindex(s))

conj(s::CompactSequence) = CompactSequence(conj(s.a), firstindex(s))

moment(s::CompactSequence,i) = sum(s[k]*k^i for k in firstindex(s):lastindex(s))

support(s::CompactSequence) = (s.offset, s.offset + s.n - 1)

support(s::CompactSequence, j::Int, k::Int) = (1/(1<<j)*(support(s)[1]+k), 1/(1<<j)*(support(s)[2]+k))

Base.sum(s::CompactSequence) = sum(s.a)

widen(s::CompactSequence) = CompactSequence(widen(s.a), s.offset)
widen(a::Array{T,N}) where {T,N} = convert(Array{widen(T),N}, a)

for op in (:+, :-, :/, :*)
  @eval ($op)(x::Number, s::CompactSequence) = CompactSequence(($op)(s.a,x),s.offset)
  @eval ($op)(s::CompactSequence, x::Number) = ($op)(x, s)
end

#########################################
# FixedSequence
#########################################


"""
A FixedSequence has a fixed length `L` starting at a certain offset `OFS`.
It contains its values in a fixed size array.
"""
struct FixedSequence{L,OFS,T} <: Sequence
    a       ::  SVector{L,T}
end

FixedSequence(a::SVector{L,T}, ::Type{Val{OFS}}) where {L,T,OFS} = FixedSequence{L,OFS,T}(a)

# type unstable constructor
FixedSequence(a::AbstractVector, offset = 0) = FixedSequence(SVector(a...), Val{offset})

# type unstable constructor
FixedSequence(s::CompactSequence) = FixedSequence(s.a, s.offset)

convert(::Type{CompactSequence{T}}, s::FixedSequence{L,OFS,T}) where {L,OFS,T} = CompactSequence(Array(s.a), OFS)

eltype(::Type{FixedSequence{L,OFS,T}}) where {L,OFS,T} = T

sublength(s::FixedSequence{L,OFS,T}) where {L,OFS,T} = L

getindex(s::FixedSequence{L,OFS,T}, i) where {L,OFS,T} = s.a[i-OFS+1]

firstindex(s::FixedSequence{L,OFS,T}) where {L,OFS,T} = OFS

lastindex(s::FixedSequence{L,OFS,T}) where {L,OFS,T} = OFS+L-1

shift(s::FixedSequence{L,OFS,T}, k::Int) where {L,OFS,T} = FixedSequence(s.a, Val{OFS+k})


hascompactsupport(::Type{FixedSequence{L,OFS,T}}) where {L,OFS,T}= True

for op in (:adjoint, :evenpart, :oddpart, :alternating_flip, :reverse, :conj, :alternating)
    @eval $op(s::FixedSequence{L,OFS,T}) where {L,OFS,T} = FixedSequence($op(CompactSequence{T}(s)))
end
