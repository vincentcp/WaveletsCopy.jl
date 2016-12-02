# extensionsequences.jl


########################################
# Type hierarchy and general interface
########################################

"""
Any subtype of ExtensionSequence embeds an indexable vectorlike object with finite length
and extends it into a bi-infinite sequence. Indexing for the subvector starts
at 1, like it is for an actual Vector.

The ExtensionSequence acts as a mutating view. One can set elements of the
extension, and the corresponding entry of the subvector will be modified.
"""
abstract ExtensionSequence{A} <: Sequence


# We assume that the embedded vector is in the field 'a' and its length is in the field 'n'.
# Variable 'k' refers to indices of the sequence, variable 'i' to indices of the embedded vector.
#
# All subtypes should implement:
# mapindex(s::SomeSubType, k) -> return the corresponding index of the embedded vector
# imapindex(s::SomeSubType, i) -> the inverse map

eltype{A}(::Type{ExtensionSequence{A}}) = eltype(A)
eltype{E <: ExtensionSequence}(::Type{E}) = eltype(super(E))

"The subvector of the extension sequence."
subvector(s::ExtensionSequence) = s.a

"The length of the subvector of the extension sequence."
sublength(s::ExtensionSequence) = s.n
# We can't really use length(s) for this, since the length of the sequence itself is infinite.


# Default mapindex and imapindex: shift by one
# We use l as index for the subvector, and k for the extension sequence.
mapindex(s::ExtensionSequence, k) = k+1
imapindex(s::ExtensionSequence, l) = l-1

getindex(s::ExtensionSequence, k) = getindex(s.a, mapindex(s, k))

checkbounds(s::ExtensionSequence, k) = 0 <= k < s.n || BoundsError()

setindex!(s::ExtensionSequence, val, k) = setindex!(s.a, val, mapindex(s, k))


"The first index of the subvector."
first_subindex(s::ExtensionSequence) = imapindex(s, 1)

"The last index of the subvector."
last_subindex(s::ExtensionSequence) = imapindex(s, sublength(s)-1)

"Iterator over each of the embedded elements of the sequence."
each_subindex(s::ExtensionSequence) = first_subindex(s):last_subindex(s)

# Invoke a constructor of an ExtensionSequence with default values.
# Default extension: ZeroPadding.
extend{EXT <: ExtensionSequence}(a, ::Type{EXT} = ZeroPadding{A}) = EXT(a)



######################################
# Definitions of specific extensions
######################################

"""
A PeriodicExtension extends a vector 'a' periodically to a bi-infinite sequence.

The indices [0...length(a)-1] map to a. Indices outside this range are mapped
by periodization modulo length(a).

The periodic extension acts as a mutating view. One can set elements of the
extension, and the corresponding entry of 'a' will be modified.
"""
immutable PeriodicExtension{A} <: ExtensionSequence{A}
    a :: A
    n :: Int
    
    PeriodicExtension(a) = new(a, length(a))
end

PeriodicExtension{A}(a::A) = PeriodicExtension{A}(a)

# It is faster to check whether k is in the proper range first, to avoid the expensive `mod`
mapindex(s::PeriodicExtension, k) = 0 <= k < s.n ? k+1 : mod(k, s.n) + 1




"""
ZeroPadding extends a vector `a` with zeros to a bi-infinite sequence.

The indices `[0...length(a)-1]` map to `a`. Indices outside this range correspond
to zero values.

A ZeroPadding acts as a mutating view. One can set elements of the
sequence, and the corresponding entry of `a` will be modified.
"""
immutable ZeroPadding{A} <: ExtensionSequence{A}
    a :: A
    n :: Int
    
    ZeroPadding(a) = new(a, length(a))
end

ZeroPadding{A}(a::A) = ZeroPadding{A}(a)

# We override getindex to return zero outside our embedded vector.
getindex(s::ZeroPadding, k::Int) = (k < 0) || (k >= s.n) ? zero(eltype(s)) : getindex(s.a, k+1)

firstindex(s::ZeroPadding) = 0

lastindex(s::ZeroPadding) = s.n-1

hascompactsupport{A}(::Type{ZeroPadding{A}}) = True

"""
A ConstantPadding extends a vector 'a' with a given `constant` to a bi-infinite sequence.

The indices [0...length(a)-1] map to a. Indices outside this range correspond to the constant
value.

A ConstantPadding acts as a mutating view. One can set elements of the
sequence, and the corresponding entry of 'a' will be modified.
"""
immutable ConstantPadding{T,A} <: ExtensionSequence{A}
    a           ::  A
    constant    ::  T
    n           ::  Int
    
    ConstantPadding(a, constant) = new(a, constant, length(a))
end

ConstantPadding{T,A}(a::A, constant::T) = ConstantPadding{T,A}(a, constant)

# We override getindex to return the constant outside our embedded vector.
getindex(s::ConstantPadding, k::Int) = (k < 0) || (k >= s.n) ? s.constant : getindex(s.a, k+1)

firstindex(s::ConstantPadding) = 0

lastindex(s::ConstantPadding) = s.n-1





"""
UndefinedExtension is an extension sequence that is undefined outside the range of the
embedded vector. The indices [0...length(a)-1] map to a. The effect of the UndefinedExtension
is that its indices start at 0.

UndefinedExtension acts as a mutating view. One can set elements of the
sequence, and the corresponding entry of 'a' will be modified.
"""
immutable UndefinedExtension{A} <: ExtensionSequence{A}
    a   ::  A
    n   ::  Int

    UndefinedExtension(a) = new(a, length(a))
end

UndefinedExtension{A}(a::A) = UndefinedExtension{A}(a)

# We simply inherit the defaults for getindex and setindex!.




"""
A SymmetricExtension extends a vector 'a' symmetrically to a bi-infinite sequence.

The indices [0...length(a)-1] map to a. Indices outside this range are mapped
by symmetrizing.

The symmetry around each of the endpoints can be whole-point (the endpoint is not repeated)
or half-point (the endpoint is repeated). The symmetry can also be even (symmetric)
or odd (anti-symmetric).

The symmetric extension acts as a mutating view. One can set elements of the
extension, and the corresponding entry of 'a' will be modified.

Definition:

immutable SymmetricExtension{A,PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT} <: ExtensionSequence{A}

Parameters:
- A:    the type of the embedded vector
- PT_LEFT:  either :wp (whole point) or :hp (half point) near left endpoint
- PT_RIGHT: either :wp or :hp for right endpoint
- SYM_LEFT: either :odd or :even symmetry near left endpoint
- SYM_RIGHT: also :odd or :even

"""
immutable SymmetricExtension{A,PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT} <: ExtensionSequence{A}
    a :: A
    n :: Int
    
    SymmetricExtension(a) = new(a, length(a))
end

SymmetricExtension{A}(a::A) = symmetric_extension_wholepoint_even(a)

# Provide four of the sixteen combinations for convenience. Construct the others by explicitly
# specifying the type parameters symbols.
symmetric_extension_wholepoint_even{A}(a::A) = SymmetricExtension{A,:wp,:wp,:even,:even}(a)

symmetric_extension_halfpoint_even{A}(a::A) = SymmetricExtension{A,:hp,:hp,:even,:even}(a)

symmetric_extension_wholepoint_odd{A}(a::A) = SymmetricExtension{A,:wp,:wp,:odd,:odd}(a)

symmetric_extension_halfpoint_even{A}(a::A) = SymmetricExtension{A,:hp,:hp,:odd,:odd}(a)


left_parity{A,PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT}(s::SymmetricExtension{A,PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT})    = SYM_LEFT
right_parity{A,PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT}(s::SymmetricExtension{A,PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT})   = SYM_RIGHT
left_symmetry{A,PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT}(s::SymmetricExtension{A,PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT})  = PT_LEFT
right_symmetry{A,PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT}(s::SymmetricExtension{A,PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT}) = PT_RIGHT


# Compute the index by mapping any index outside the range of the embedded vector to
# an index that is closer to the interval and repeat. The recursion ends when the index
# lands inside the interval.
# For far away points this is not the fastest way, but that won't happen very often. The
# alternative is to reduce the extension by exploiting periodicity of the symmetric extension,
# but the period depends on all symmetries too and this will requires a modulo operation. For
# nearby points, the implementation below is probably faster (though I haven't actually checked this).
function mapindex(s::SymmetricExtension, k)
    if k >= s.n
        # We are on the right of the interval: use symmetry wrt right endpoint
        mapindex_right(s, k)
    elseif k < 0
        # We are on the left of the interval: use symmetry wrt left endpoint
        mapindex_left(s, k)
    else
        k+1
    end
end

# Right whole point symmetry: map across the midpoint at index n-1, then convert to index of a.
# Index in k-space is:  n-1 - (k - (n-1)) = 2n - k - 2
mapindex_right{A,PT_LEFT}(s::SymmetricExtension{A,PT_LEFT,:wp}, k) = mapindex(s, 2s.n-k-2)

# Right half point symmetry: map across the midpoint at index n-1/2, then convert to index of a.
# Index in k-space is:  n-1/2 - (k - (n-1/2)) = 2n - k - 1
mapindex_right{A,PT_LEFT}(s::SymmetricExtension{A,PT_LEFT,:hp}, k) = mapindex(s, 2s.n-k-1)

# Left whole point symmetry: map across the midpoint at index 0, then convert to index of a.
# Index in k-space is:  0 + (0-k) = -k
mapindex_left{A,PT_RIGHT}(s::SymmetricExtension{A,:wp,PT_RIGHT}, k) = mapindex(s, -k)

# Left half point symmetry: map across the midpoint at index -1/2, then convert to index of a.
# Index in k-space is:  -1/2 + (-1/2-k) = -k-1
mapindex_left{A,PT_RIGHT}(s::SymmetricExtension{A,:hp,PT_RIGHT}, k) = mapindex(s, -k-1)


# For getindex we have to use the same logic as mapindex, but now we also have to trace
# the odd/even-ness of the symmetries.
function getindex(s::SymmetricExtension, k)
    if k >= s.n
        # We are on the right of the interval: use symmetry wrt right endpoint
        getindex_right(s, k)
    elseif k < 0
        # We are on the left of the interval: use symmetry wrt left endpoint
        getindex_left(s, k)
    else
        getindex(s.a, k+1)
    end
end

# For even symmetry on both endpoints, it is simple: the sign never flips. Short circuit.
# This probably does not gain much, as mapindex still does the recursion anyway...
getindex{A,PT_LEFT,PT_RIGHT}(s::SymmetricExtension{A,PT_LEFT,PT_RIGHT,:even,:even}, k) = getindex(s.a, mapindex(s, k))

# Right whole point symmetry
getindex_right{A,PT_LEFT,SYM_LEFT}(s::SymmetricExtension{A,PT_LEFT,:wp,SYM_LEFT,:even}, k) = getindex(s, 2s.n-k-2)
getindex_right{A,PT_LEFT,SYM_LEFT}(s::SymmetricExtension{A,PT_LEFT,:wp,SYM_LEFT,:odd}, k) = -getindex(s, 2s.n-k-2)

# Right half point symmetry
getindex_right{A,PT_LEFT,SYM_LEFT}(s::SymmetricExtension{A,PT_LEFT,:hp,SYM_LEFT,:even}, k) = getindex(s, 2s.n-k-1)
getindex_right{A,PT_LEFT,SYM_LEFT}(s::SymmetricExtension{A,PT_LEFT,:hp,SYM_LEFT,:odd}, k) = -getindex(s, 2s.n-k-1)

# Left whole point symmetry
getindex_left{A,PT_RIGHT}(s::SymmetricExtension{A,:wp,PT_RIGHT,:even}, k) = getindex(s, -k)
getindex_left{A,PT_RIGHT}(s::SymmetricExtension{A,:wp,PT_RIGHT,:odd}, k) = -getindex(s, -k)

# Left half point symmetry
getindex_left{A,PT_RIGHT}(s::SymmetricExtension{A,:hp,PT_RIGHT,:even}, k) = getindex(s, -k-1)
getindex_left{A,PT_RIGHT}(s::SymmetricExtension{A,:hp,PT_RIGHT,:odd}, k) = -getindex(s, -k-1)


