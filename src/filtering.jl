# filtering.jl

module Sampling

export  filter_and_downsample!, filter_and_downsample_bruteforce!,
        upsample_and_filter!, upsample_and_filter_bruteforce!,
        fstate_length_downsampling, fstate_length_upsampling,
        fstate_downsampling, fstate_upsampling,
        filter_and_downsample_history, upsample_and_filter_history


# The next couple of routines are convenience functions for easier use of the sampling routines.
# For an explanation of the formulas, see the comments in filter_and_downsample! and in
# upsample_and_filter!.

"The length of the filter state vector for a filter of length L and downsampling by a factor of M."
fstate_length_downsampling(L, M) = div(L+M-2, M) # or ceil(Int,(L-1)/M)

"A vector suitable to use as fstate for filter_and_downsample."
fstate_downsampling(T, L, M) =
    zeros(T, fstate_length_downsampling(L, M))



"The number of entries in x that are history when filtering and downsampling. These are used
in the computation of the first output element. After these, every M samples in x yield another
output value."
filter_and_downsample_history(L, M) = L-1


"""
`function filter_and_downsample!(y, x, filter, M[, iter_y, iter_x, fstate])`

Convolve the given `filter` with the data in vector `x`, and store the result in `y`. The output
is subsampled by a factor of `M`. The filter is assumed to be causal and given by `L` coefficients
`filter[1..L]`. Note that the output of a non-causal filter can be obtained by shifting the input
values.

The first value of `y` uses the first `L` values of `x`. The second value of `y` uses `M` more
samples, and so on. Thus, vector `x` should have length `L + (N-1)M` in total. The first `L-1`
values of `x` may be considered to be history before the first sample point of `y`.

In that setting, the precise computation of this routine is:

`y[k] = sum([ h[l] * x[L+(k-1)M-l+1] for l = 1:L])`

Alternatively, the values of `y` and `x` to be used can be specified with iterators `iter_y` and
`iter_x`. In that case, the iterators should have length `N` and `L + (N-1)M`, while `y` and `x`
can have any size and dimension. The iterators will be traversed exactly once.
(Note: this interface may change once views become more efficient in Julia.)

A filter state vector `fstate` can be supplied, in which case this routine does not allocate memory.
The length of the vector can be computed with the routine `fstate_length_downsampling`.
"""
function filter_and_downsample!(y, x, filter, M, iter_y = 1:length(y), iter_x = 1:length(x),
    fstate = fstate_downsampling(eltype(x), length(filter), M))

    # We will compute N coefficients
    N = length(iter_y)
    L = length(filter)

    if L == 1
        # Save ourselves the trouble of the algorithm below when the filter has length 1.
        x_state = start(iter_x)
        y_state = start(iter_y)

        # We let m range from 0 to M-1 over and over again, and we write x to y whenever m == 0.
        m = 0
        for n = 1:length(iter_x)
            (i_x,x_state) = next(iter_x, x_state)

            if m == 0
                xk = x[i_x]
                (i_y,y_state) = next(iter_y, y_state)
                y[i_y] = filter[1] * xk
            end

            m <= M-2 ? m += 1 : m = 0
        end
        return
    end

    # Let's be strict in the size of x for the time being, to catch indexing bugs quicker.
    @assert length(iter_x) == L+(N-1)*M

    # Each value of y is a linear combination of values of x. In order to sample each value
    # of x only once, we keep track of partial sums in the state vector fstate.
    # We have Q bins in fstate, containing a partial sum of (at most) M x-values each.
    Q = fstate_length_downsampling(L, M)
    # Allow longer state vectors, so that one sufficiently large vector can be used for all
    @assert length(fstate) >= Q

    # The last bin contains only 1<=R<=M values
    R = (L-1) - (Q-1)*M
    # Hence: L-1 = (Q-1)*M + R

    # We have to initialize the filter state before we can start. We will need L
    # samples of x for the computation of the first value of y.
    # We read in the first set of L-1 values of x and distribute them over fstate.
    #
    # In our description below, we use the notation y_k with indices starting from zero.
    # So, in terms of the Julia vectors we have: y_k = y[k+1].
    #
    # The first element of y will be y_0:
    # y_0 = \sum_{l=0}^{L-1} h_l x_{p-1}
    # with p = L-1.
    # The next element of y is M samples down the road:
    # y_1 = \sum_{l=0}^{L-1} h_l x_{p+M-1}
    # and so on.
    #
    # We keep the current partial sum for y_0 in fstate_0, the partial sum for y_1 in fstate_1,
    # and so on. The first L-1 samples of x are x_k for k=0..L-2.
    # Thus, we want to reach the following state:
    # fstate_0 = \sum_{l=1}^{L-1} h_l x_{p-l}
    # fstate_1 = \sum_{l=M+1}^{L-1} h_l x_{p+M-l}
    # ...
    # fstate_{Q-1} = \sum_{l=(Q-1)M+1}^{L-1} h_l x_{p+(Q-1)M-l}
    #
    # But we want to indices in these sums also reversed, since we compute the samples x_i
    # one by one, i=0,1,...,p-1. Where should x_k go, and with which filter coefficient?
    # To that end, after a substition we obtain:
    # fstate_0 = \sum_{k=p-L+1}^{p-1} h_{p-k} x_k
    # fstate_1 = \sum_{k=p-L+1+M}^{p-1} h_{p+M-k} x_k
    # ...
    # fstate_{Q-1} = \sum_{k=p-L+1+(Q-1)M}^{p-1} h_{p+(Q-1)M-k} x_k
    # For p = L-1, the summation for fstate_q simply starts at q*M.
    p = L-1

    # Zero out fstate first
    for q = 0:Q-1
        fstate[q+1] = 0
    end

    # Loop over the first L-1 samples of x
    # We use that L-1 = (Q-1)*M+R in order to avoid having to compute modulo's
    x_state = start(iter_x)
    k = 0
    for i = 0:Q-1
        for j = 0:M-1
            # Break out of the loop if we have reached (Q-1)*M+R
            if (i == Q-1) && (j==R)
                break
            end

            # Starting at index 0, this is sample x_{k} = x_{iQ+j}
            (i_x,x_state) = next(iter_x, x_state)
            xk = x[i_x]

            # Since the n-th bin starts with sample x_{qM}, add it to the first i bins.
            for q = 0:i
                fstate[q+1] += filter[p+q*M-k+1] * xk
            end
            k += 1
        end
    end
    # Helpful output for debugging:
    # println(fstate)

    # Now we can start the loop over the elements of y
    y_state = start(iter_y)
    for n in 0:N-1
        # In this loop we will take M more samples from x.
        # The first sample results in a new value for y. The course of action for
        # this sample is similar to the M=1 case (without subsampling). The other samples
        # only have to be distributed over fstate in the right way.
        # For the last value of n, we only need one more sample of x instead of M.

        for j = 0:M-1
            # Sample the next value of x
            # Starting at index 0, this is sample x_{k} = x_{p+nM+j} = x_{L-1+nM+j}
            (ix,x_state) = next(iter_x, x_state)
            xk = x[ix]

            # As announced, the first sample requires special treatment for each n
            if j == 0
                # First, compute the next value for y
                (i_y,y_state) = next(iter_y, y_state)
                y[i_y] = fstate[1] + filter[1] * xk

                # If this was the last value of y, we are done and we can stop right here
                if n == N-1
                    break
                end

                # Now we can shift the values in fstate, the first one is no longer needed
                for q = 0:Q-2
                    fstate[q+1] = fstate[q+2]
                end
                fstate[Q] = 0
            end

            # At the start (j=0) of this loop, the vector fstate is like the equations above
            # with p = L-1 + nM. Now we want to make fstate ready for p = L-1 + (n+1)M.
            # At this point we have k = L-1+nM+j.
            # We have the same summation as before, relative to p:
            # fstate_q = \sum_{k=p-L+1+qM}^{p-1} h_{p+qM-k} x_k
            # and we see that the coefficient with x_k is h_{p+qM-k} = h_{(q+1)M-j}.
            for q = 0:Q-2
                fstate[q+1] += filter[(q+1)*M-j+1] * xk
            end
            if j >= M-R
                fstate[Q] += filter[Q*M-j+1] * xk
            end
            # We don't actually need k:
            # k += 1

            # Helpful output for debugging:
            # println(fstate)
        end
    end
end


"""
Brute force computation of applying a filter followed by downsampling. The result should be the same as
the (more efficient) routine filter_and_downsample!. This function is implemented for testing purposes.
"""
function filter_and_downsample_bruteforce!(y, x, filter, M, iter_y = 1:length(y), iter_x = 1:length(x))
    L = length(filter)

    # Collect the values of x from the iterator
    x_values = eltype(x)[x[i] for i in iter_x]

    # Apply the filter to compute values of y
    y_values = zeros(eltype(y), length(iter_y))
    for k in 1:length(iter_y)
        y_values[k] = 0
        for l = 1:L
            y_values[k] += filter[l] * x_values[L+(k-1)*M-l+1]
        end
    end

    # And move the results to the given vector y
    k = 0
    for i in iter_y
        k += 1
        y[i] = y_values[k]
    end
end



"The number of entries in x that are history when upsampling and filtering. These are used
in the computation of the first output element. After these, every sample of x yields M additional
output values."
function upsample_and_filter_history(L, M, M0 = 0)
    # An example layout of z after time k=0 is as follows, for M = 3 and M0 = 1:
    # k  :  0    1    2    3    4
    # z_k:  0    x0   0    0    x1
    # We want to compute the number of samples of x in the L-1 samples of z preceding k = 0.
    # First, we compute L-1 = Q1 * M + R1. There are at least Q1 groups of M samples, with
    # one sample of x in each. We have R1 additional samples at the start, and M0 additional
    # samples (all zero) starting at time 0. Thus, if R1+M0 is greater than or equal to M, there
    # will be an additional sample of x in the history at the start of z.
    Q1 = div(L-1,M)
    R1 = L - 1 - Q1*M
    R1+M0 >= M ? Q1+1 : Q1
end

"The number of leading zeros in the upsampled vector z that is used in upsample_and_filter!."
function upsample_and_filter_zshift(L, M, M0 = 0)
    # We compute the number of zeros at the start of the upsampled vector z. For the reasoning,
    # see the comments in upsample_and_filter_history.
    # There are Q1 groups of M samples of z, and R1+M0 additional samples. The latter value modulo
    # M is the number of leading zeros.
    Q1 = div(L-1,M)
    R1 = L-1 - Q1*M
    mod(R1+M0,M)
end

"The lenght of the filter state vector for a filter of length L and upsampling by a factor of M."
fstate_length_upsampling(L, M) = L-1

"A vector suitable to use as fstate for upsample_and_filter."
fstate_upsampling(T, L, M) =
    zeros(T, fstate_length_upsampling(L, M))



"""
`function upsample_and_filter!(y, x, filter, M[, M0, iter_y, iter_x, fstate])`

Upsample the data in vector `x` by a factor of `M`, and convolve with the given `filter`. The result
is *added* to the vector `y`. The filter is assumed to be causal and given by `L` coefficients
`filter[1..L]`.

Upsampling is defined as placing `M-1` zeros between two consecutive samples of `x`. This leads to an
upsampled vector z. The first output value of `y` is computed from the first `L` samples of `z`, and
each additional value in `z` yields another value in `y`. Thus, the first `L-1` samples of `z` can be
considered to be history.

The upsampling is aligned such that `z[L]` corresponds to a value of `x`, say `x_0`. It is followed by
`M-1` zeros in `z[L+1...L+M-1]`. Optionally, the location of `x_0` can be shifted by `M0` places. The number
of samples of `x` in history, i.e. in `z[1...L-1]` can be computed with `upsample_and_filter_history(L,M,M0)`.
It is assumed that this many values are given before `x1` in the input vector `x`.

To achieve this alignment, the vector `z` starts with a number of leading zeros, rather than with the first
value in `x`. The number of zeros can be computed with `upsample_and_filter_zshift`.

For completeness, we precisely describe the computations of this function. Conceptually the vector `x` is
upsampled to a vector `z` as follows:

`z[1+r+iM] = x[1+i]`

for i = 0,1,... , while all other values are 0. Here, `r` is the number of leading zeros (`zshift`) of `z`.

Let `L` be the length of the given `filter`, and let `L = QM + R` with `0 <= R < M`. Each value `y[k]`
uses `L` values of `z`, which depending on the alignment corresponds to either `Q` or
`Q+1` values of `x`. The precise computation of this routine is:

`y[k] = sum([ h[l] * z[L+k-l] for l = 1:L])`

Alternatively, the values of `y` and `x` to be used can be specified with iterators `iter_y`
and `iter_x`. These iterators will be traversed exactly once.

A filter state vector `fstate` can be supplied, in which case this routine does not allocate memory.
The minimal length of the vector can be computed with the routine `fstate_length_upsample`.
"""
function upsample_and_filter!(y, x, filter, M, M0 = 0, iter_y = 1:length(y), iter_x = 1:length(x),
    fstate = fstate_upsampling(eltype(x), length(filter), M))

    @assert 0 <= M0 < M

    # We will compute N coefficients
    N = length(iter_y)
    L = length(filter)

    # The number of zeros at the start of the z-vector
    zshift = upsample_and_filter_zshift(L, M, M0)

    if L == 1
        # Save ourselves the trouble of the algorithm below when the filter has length 1.
        x_state = start(iter_x)
        y_state = start(iter_y)

        # The upsampling complicates things a little. We let m range from 0 to M-1 over and over again,
        # taking a sample from x whenever m == 0. The starting value of m is the zshift, i.e. the number of
        # leading zeros in the z-vector.
        m = zshift
        for n = 1:length(iter_y)
            (i_y,y_state) = next(iter_y, y_state)

            if m == 0
                (i_x,x_state) = next(iter_x, x_state)
                xk = x[i_x]
                y[i_y] += filter[1] * xk
            else
                # Of course the line below does not do anything, but we keep it to remember
                # that y[i_y] should be set to zero if we are not adding to y, but storing the
                # result in y instead. I.e., for a non-additive version the line below should
                # be: y[i_y] = 0.
                y[i_y] += 0
            end

            m <= M-2 ? m += 1 : m = 0
        end
        return
    end

    # Read the documention of filter_and_downsample! first. Here, we focus on the differences.
    #
    # As before, each value of y is a linear combination of values of x, and we keep track of partial
    # sums in the state vector fstate. This time around, the length of the state vector is L-1, because
    # we are not downsampling y.
    fslen = L-1
    @assert length(fstate) >= fslen

    # L = QM + R
    Q = div(L, M)
    R = L - Q*M

    # and L-1 = Q1*M + R1
    Q1 = div(L-1,M)
    R1 = L - 1 - Q1*M


    # We have to initialize the filter state before we can start. In this case, each bin corresponds to
    # a future value of y, containing either Q or Q+1 values of x. At the start, we want to distribute
    # L-1 samples of the upsampled vector z over the bins. How many x samples are included is calculated
    # in a separate routine:
    Nx = upsample_and_filter_history(L, M, M0)

    # TODO: check the length of iter_x. What should be its minimal length?
#    @assert length(iter_x) >= Nx1 + div(N,M)

    # The bins contain partial sums of z, say up to z[p] = z_{p-1} where p = L-1.
    p = L-1

    # We now have
    # fstate_0 = \sum_{l=1}^{L-1} h_l z_{p-l}
    # fstate_1 = \sum_{l=2}^{L-1} h_l z_{p+1-l}
    # ...
    # fstate_s = \sum_{l=s+1}^{L-1} h_l z_{p+s-l}
    # ...
    # fstate_{L-2} = \sum_{l=L-1}^{L-1} h_l z_{p+(L-2)-l}
    #
    # Substituting indices (p+s-l=t) we sum over consecutive values of z instead:
    # fstate_0 = \sum_{t=p-L+1}^{p-1} h_{p-t} z_t
    # fstate_1 = \sum_{l=p-L+2}^{p-1} h_{p+1-t} z_t
    # ...
    # fstate_s = \sum_{t=p-L+s+1}^{p-1} h_{p+s-t} z_t
    # ...
    # fstate_{L-2} = \sum_{t=p-1}^{p-1} h_{p+L-2-t} z_t
    #
    # That was not too bad, but it is more complicated in terms of x. We have z_{zshift+M*k} = x_k, so
    # so we need to select from the above sums only those terms where t is a multiple of M plus zshift.

    # Zero out fstate first
    for s = 1:fslen
        fstate[s] = 0
    end

    # Now draw the first Nx samples of x
    x_state = start(iter_x)
    for k = 0:Nx-1
        # We draw sample x_k
        (i_x,x_state) = next(iter_x, x_state)
        xk = x[i_x]

        # Sample x_k corresponds to z_{r+M*k}, where r is the zshift. It goes in bin s if r+M*k
        # is in the range of the sum for fstate_s above.
        # Note that at this stage p = L-1, and hence the starting index p-L+s+1 is simply s. This means
        # that bin s = r+M*k is the last bin that includes our current sample x_k.
        for s = 0:zshift+M*k
            fstate[s+1] += filter[p+s-M*k-zshift+1] * xk
        end
    end

    # Helpful output for debugging:
    # println(fstate)

    # Now we can start the main loop over the elements of y
    y_state = start(iter_y)
    # We have to sample x only once every M values of y. We let a variable m go from 0 to M-1,
    # and sample x whenever m == 0. By construction, we have M0 zeros to go before the first
    # sample, so we initialize m to M-M0 (modulo M).
    m = mod(M-M0,M)
    for n in 0:N-1
        if m == 0
            # We draw sample z[p] = z[L+n] (= z_{p-1} = z_{L-1+n}, remember we start from 0 for z_j)
            # This equals x_k.
            (i_x,x_state) = next(iter_x, x_state)
            xk = x[i_x]

            # We don't actually need k anymore, but here it is:
            # k += 1

            # We can compute the y value and shift the bins
            (i_y,y_state) = next(iter_y, y_state)
            y[i_y] += fstate[1] + filter[1] * xk
            for s = 0:fslen-2
                fstate[s+1] = fstate[s+2]
            end
            fstate[fslen] = 0

            # In which bin does xk go? At this stage, the formulas higher up for fstate_s are valid for p = L-1+n.
            # Now, they should hold for p=L+n (ready for the next sample). That means we have to shift the bins,
            # and then add the sample with the right filter coefficient to each bin.
            # In the sum for fstate_s, t ranges from p-L+s+1 = n+s+1 to p-1 = L+n-1.
            # The current index of z equals the upper bound. In other words, the value goes in _all_ the bins.
            # The filter coefficient is h_{p+s-t} for t=p-1
            for s = 0:fslen-1
                fstate[s+1] += filter[s+2] * xk
            end
        else
            # No need to update the bins since we have no new x-value, but assign y and shift the bins
            (i_y,y_state) = next(iter_y, y_state)
            y[i_y] += fstate[1]
            for s = 0:fslen-2
                fstate[s+1] = fstate[s+2]
            end
            fstate[fslen] = 0
        end

        # Update m and set to zero when the new value would equal M.
        m < M-1 ? m += 1 : m = 0

        # Helpful output for debugging:
        # println(fstate)
    end
end


"""
Brute force computation of upsampling followed by applying a filter. The result should be the same as
the (more efficient) routine upsample_and_filter!. This function is implemented for testing purposes.
"""
function upsample_and_filter_bruteforce!(y, x, filter, M, M0 = 0,
    iter_y = 1:length(y), iter_x = 1:length(x))

    L = length(filter)

    # Collect the values of x from the iterator
    x_values = eltype(x)[x[i] for i in iter_x]

    # Brute force compute the upsampled vector
    zshift = upsample_and_filter_zshift(L, M, M0)
    z = zeros(eltype(x), 1 + zshift + M * (length(x_values)-1))
    z[1+zshift:M:end] = x_values

    # Apply the filter to obtain values of y
    y_values = zeros(eltype(y), length(iter_y))
    for k in 1:length(iter_y)
        y_values[k] = 0
        for l = 1:L
            y_values[k] += filter[l] * z[L+k-l]
        end
    end

    # Copy them to the right place in y
    k = 0
    for i in iter_y
        k += 1
        y[i] = y_values[k]
    end
end




end #module
