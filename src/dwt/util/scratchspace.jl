###############################################################################
## scratch util function ######################################################
###############################################################################

# Length functions
@inline evaluate_periodic_in_dyadic_points_length(side::Side, kind::Kind, w, j, k, d) =
    evaluate_periodic_in_dyadic_points_length(d)
@inline evaluate_periodic_in_dyadic_points_length(d) = 1<<d
@inline evaluate_periodic_in_dyadic_points_scratch_length(side::Side, kind::Kind, w, j, k, d) =
    evaluate_in_dyadic_points_length(side, kind, w, j, k, d)
@inline evaluate_periodic_in_dyadic_points_scratch2_length(side::Side, kind::Kind, w, j, k, d) =
    evaluate_in_dyadic_points_scratch_length(side, kind, w, j, k, d)
@inline evaluate_periodic_in_dyadic_points_scratch3_length(side::Side, kind::Kind, w, j, k, d) =
    evaluate_in_dyadic_points_scratch2_length(side, kind, w, j, k, d)
@inline evaluate_in_dyadic_points_length(side::Side, kind::Kind, w, j, k, d) = ((d-j) >= 0) ?
    (1<<(d-j))*support_length(side, kind, w)+1 :
    cld(support_length(side, kind, w)+1, 1<<(j-d))
@inline evaluate_in_dyadic_points_scratch_length(side::Side, kind::Kind, w, j, k, d) = ((d-j) >= 0) ?
    0 :
    support_length(side, kind, w)+1
@inline evaluate_in_dyadic_points_scratch_length(side::Side, kind::Wvl, w, j, k, d) = ((d-j) >= 1) ?
    evaluate_in_dyadic_points_length(side, Scl(), w, j+1, k, d) :
    evaluate_in_dyadic_points_scratch_length(side, kind, w, j, k, j+1)
@inline evaluate_in_dyadic_points_scratch2_length(side::Side, kind::Wvl, w, j, k, d) = ((d-j) >= 1) ?
    evaluate_in_dyadic_points_scratch_length(side, Scl(), w, j+1, k, d) :
    evaluate_in_dyadic_points_length(side, kind, w, j, k, j+1)

# combining length funtions
@inline function get_evaluate_periodic_scratch_space_length(side::Side, kind::Scl, w::DiscreteWavelet, j::Int, k::Int, d::Int)
    f = evaluate_periodic_in_dyadic_points_length(side, kind, w, j, k, d)
    scratch = evaluate_periodic_in_dyadic_points_scratch_length(side, kind, w, j, k, d)
    scratch2 = evaluate_periodic_in_dyadic_points_scratch2_length(side, kind, w, j, k, d)
    f, scratch, scratch2
end

@inline function get_evaluate_periodic_scratch_space_length(side::Side, kind::Wvl, w::DiscreteWavelet, j::Int, k::Int, d::Int)
    f = evaluate_periodic_in_dyadic_points_length(side, kind, w, j, k, d)
    scratch = evaluate_periodic_in_dyadic_points_scratch_length(side, kind, w, j, k, d)
    scratch2 = evaluate_periodic_in_dyadic_points_scratch2_length(side, kind, w, j, k, d)
    scratch3 = evaluate_periodic_in_dyadic_points_scratch3_length(side, kind, w, j, k, d)
    f, scratch, scratch2, scratch3
end

@inline function get_evaluate_scratch_space_length(side::Side, kind::Kind, w::DiscreteWavelet{T}, j::Int, k::Int, d::Int) where {T}
    f = evaluate_in_dyadic_points_length(side, kind, w, j, k, d)
    scratch = evaluate_in_dyadic_points_scratch_length(side, kind, w, j, k, d)
    f, scratch
end

@inline function get_evaluate_scratch_space_length(side::Side, kind::Wvl, w::DiscreteWavelet{T}, j::Int, k::Int, d::Int) where {T}
    f = evaluate_in_dyadic_points_length(side, kind, w, j, k, d)
    scratch = evaluate_in_dyadic_points_scratch_length(side, kind, w, j, k, d)
    scratch2 = evaluate_in_dyadic_points_scratch2_length(side, kind, w, j, k, d)
    f, scratch, scratch2
end

"All lengths for a wavelet evaluation with l levels"
function get_evaluate_scratch_space_lengths(side::Side, w::DiscreteWavelet{T}, l::Int, d::Int) where {T}
    lengths = Set{Int}()
    for index in wavelet_indices(l)
        push!(lengths,get_evaluate_scratch_space_length(side, kind(index), w, level(index), offset(index), d)...)
    end
    sort(collect(lengths))
end

"All lengths for a periodic wavelet evaluation with l levels"
function get_evaluate_periodic_scratch_space_lengths(side::Side, w::DiscreteWavelet{T}, l::Int, d::Int) where {T}
    lengths = Set{Int}()
    for index in wavelet_indices(l)
        push!(lengths,get_evaluate_periodic_scratch_space_length(side, kind(index), w, level(index), offset(index), d)...)
    end
    sort(collect(lengths))
end


# allocating functions
@inline get_evaluate_periodic_scratch_space(side::Side, kind::Kind, w::DiscreteWavelet{T}, j::Int, k::Int, d::Int) where {T} =
    [zeros(T, l) for l in get_evaluate_periodic_scratch_space_length(side, kind, w, j, k, d)]

@inline get_evaluate_scratch_space(side::Side, kind::Kind, w::DiscreteWavelet{T}, j::Int, k::Int, d::Int) where {T} =
    [zeros(T, l) for l in get_evaluate_scratch_space_length(side, kind, w, j, k, d)]

"All scratch space ever necessary for a wavelet evaluation with l levels"
get_evaluate_scratch_space(side::Side, w::DiscreteWavelet{T}, l::Int, d::Int) where {T} =
    [zeros(T, l) for l in get_evaluate_scratch_space_lengths(side, w, l, d)]

"All scratch space ever necessary for a periodic wavelet evaluation with l levels"
get_evaluate_periodic_scratch_space(side::Side, w::DiscreteWavelet{T}, l::Int, d::Int) where {T} =
    [zeros(T, l) for l in get_evaluate_periodic_scratch_space_lengths(side, w, l, d)]

"A container for scratch space vectors"
struct ScratchSpace{T}
    scratch::Vector{Vector{T}}
    lengths::Vector{Int}
end

ScratchSpace(s::Vector{Vector{T}}) where {T} = ScratchSpace{T}(s, map(length, s))

"Retrieve vectors of given lengths of the scratch space container."
get_scratch_space(SS::ScratchSpace, lengths) =
    [SS.scratch[index]  for index in _get_indices(SS.lengths, lengths)]

function _get_indices(lengths::Vector{Int}, input)
    L = length(lengths)
    LL = length(input)
    r = Array{Int}(undef, LL)
    @inbounds for i in 1:LL
        s = input[i]
        for j in 1:L
            if lengths[j] == s
                r[i] = j
                break
            end
        end
    end
    r
end

"Scratch space for a wavelet evaluation with l levels"
EvalScratchSpace(side::Side, w::DiscreteWavelet, l::Int, d::Int) =
    ScratchSpace(get_evaluate_scratch_space(side, w, l, d))

"Scratch space for periodic wavelet evaluation with l levels"
EvalPeriodicScratchSpace(side::Side, w::DiscreteWavelet, l::Int, d::Int) =
    ScratchSpace(get_evaluate_periodic_scratch_space(side, w, l, d))

"Scratch space for one particular wavelet evaluation"
EvalScratchSpace(side::Side, kind::Kind, w::DiscreteWavelet, j::Int, k::Int, d::Int) =
    ScratchSpace(get_evaluate_scratch_space(side, kind, w, j, k, d))

"Scratch space for one particular periodic wavelet evaluation"
EvalPeriodicScratchSpace(side::Side, kind::Kind, w::DiscreteWavelet, j::Int, k::Int, d::Int)  =
    ScratchSpace(get_evaluate_periodic_scratch_space(side, kind, w, j, k, d))

"Retrieve scratch space needed for a wavelet evaluation from the structure ScratchSpace"
get_evaluate_scratch_space(SS::ScratchSpace, side::Side, kind::Kind, w::DiscreteWavelet, j::Int, k::Int, d::Int) =
    get_scratch_space(SS, get_evaluate_scratch_space_length(side, kind, w, j, k , d))

"Retrieve scratch space needed for a periodic wavelet evaluation from the structure ScratchSpace"
get_evaluate_periodic_scratch_space(SS::ScratchSpace, side::Side, kind::Kind, w::DiscreteWavelet, j::Int, k::Int, d::Int) =
    get_scratch_space(SS, get_evaluate_periodic_scratch_space_length(side, kind, w, j, k , d))
