# iterators.jl

module Iterators

"An iterator that is longer than the vector it iterates over, and iterates it periodically between given bounds."
immutable PeriodizingIterator
# TODO: try out whether this might be faster than embedding a vector in a periodic sequence
end


"An iterator that is longer than the vector it iterates over, and iterates it symmetrically between given bounds."
immutable SymmetrizingIterator
# TODO
end



end
