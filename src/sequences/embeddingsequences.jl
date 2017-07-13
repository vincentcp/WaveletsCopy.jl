# embeddingsequence.jl

"""
An EmbeddingSequence is like an ExtensionSequence, except that it does not store a vector.
Instead, one passes a vector along with each `getindex` call.
"""
abstract type EmbeddingSequence <: Sequence end


struct PeriodicEmbedding <: EmbeddingSequence
end

getindex(s::PeriodicEmbedding, x, k) = 0 <= k < length(x) ? x[k+1] : x[mod(k, length(x)) + 1]


struct SymmetricEmbedding{PT_LEFT,PT_RIGHT,SYM_LEFT,SYM_RIGHT} <: EmbeddingSequence
end

struct FunctionEmbedding <: EmbeddingSequence
  f :: Function
end

getindex(s::FunctionEmbedding, x, k) = s.f(k)


SymmetricEmbedding_wholepoint_even() = SymmetricEmbedding{:wp,:wp,:even,:even}()

SymmetricEmbedding_halfpoint_even() = SymmetricEmbedding{:hp,:hp,:even,:even}()


function getindex(s::SymmetricEmbedding, x, k)
    if k >= length(x)
        # We are on the right of the interval: use symmetry wrt right endpoint
        getindex_right(s, x, k)
    elseif k < 0
        # We are on the left of the interval: use symmetry wrt left endpoint
        getindex_left(s, x, k)
    else
        x[k+1]
    end
end

# Right whole point symmetry
getindex_right{PT_LEFT,SYM_LEFT}(s::SymmetricEmbedding{PT_LEFT,:wp,SYM_LEFT,:even}, x, k) = getindex(s, x, 2*length(x)-k-2)
getindex_right{PT_LEFT,SYM_LEFT}(s::SymmetricEmbedding{PT_LEFT,:wp,SYM_LEFT,:odd}, x, k) = -getindex(s, x, 2*length(x)-k-2)

# Right half point symmetry
getindex_right{PT_LEFT,SYM_LEFT}(s::SymmetricEmbedding{PT_LEFT,:hp,SYM_LEFT,:even}, x, k) = getindex(s, x, 2*length(x)-k-1)
getindex_right{PT_LEFT,SYM_LEFT}(s::SymmetricEmbedding{PT_LEFT,:hp,SYM_LEFT,:odd}, x, k) = -getindex(s, x, 2*length(x)-k-1)

# Left whole point symmetry
getindex_left{PT_RIGHT}(s::SymmetricEmbedding{:wp,PT_RIGHT,:even}, x, k) = getindex(s, x, -k)
getindex_left{PT_RIGHT}(s::SymmetricEmbedding{:wp,PT_RIGHT,:odd}, x, k) = -getindex(s, x, -k)

# Left half point symmetry
getindex_left{PT_RIGHT}(s::SymmetricEmbedding{:hp,PT_RIGHT,:even}, x, k) = getindex(s, x, -k-1)
getindex_left{PT_RIGHT}(s::SymmetricEmbedding{:hp,PT_RIGHT,:odd}, x, k) = -getindex(s, x, -k-1)


struct CompactEmbedding <: EmbeddingSequence
    offset  ::  Int
end

getindex(s::CompactEmbedding, x, i) = x[s.offset+i+1]
