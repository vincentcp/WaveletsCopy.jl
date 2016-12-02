# wvlt_cdf.jl

"The Cohen-Daubechies-Feauveau family of compactly supported biorthogonal wavelets."
immutable CDFWavelet{P,Q,T} <: DiscreteWavelet{T}
end

is_symmetric{P,Q,T}(::Type{CDFWavelet{P,Q,T}}) = True

is_biorthogonal{P,Q,T}(::Type{CDFWavelet{P,Q,T}}) = True

T0 = Float64

cdf11 = CDFWavelet{1,1,T0}()
cdf13 = CDFWavelet{1,3,T0}()
cdf15 = CDFWavelet{1,5,T0}()

cdf22 = CDFWavelet{2,2,T0}()
cdf24 = CDFWavelet{2,4,T0}()
cdf26 = CDFWavelet{2,6,T0}()

cdf31 = CDFWavelet{3,1,T0}()
cdf33 = CDFWavelet{3,3,T0}()
cdf35 = CDFWavelet{3,5,T0}()

cdf42 = CDFWavelet{4,2,T0}()
cdf44 = CDFWavelet{4,4,T0}()
cdf46 = CDFWavelet{4,6,T0}()

cdf51 = CDFWavelet{5,1,T0}()
cdf53 = CDFWavelet{5,3,T0}()
cdf55 = CDFWavelet{5,5,T0}()

cdf62 = CDFWavelet{6,2,T0}()
cdf64 = CDFWavelet{6,4,T0}()
cdf66 = CDFWavelet{6,6,T0}()


# Explicit listing of some coefficient follows.
# The format is a tuple: (n,h). The filter coefficients are sqrt(2)/n * h.

cdf11_htilde = (2, [1, 1])
cdf13_htilde = (16, [-1, 1, 8, 8, 1, -1])
cdf15_htilde = (256, [3, -3, -22, 22, 128, 128, 22, -22, -3, 3])


cdf22_htilde = (8, [-1, 2, 6, 2, -1])
cdf24_htilde = (128, [3, -6, -16, 38, 90, 38, -16, -6, 3])
cdf26_htilde = (1024, [-5, 10, 34, -78, -123, 324, 700, 324, -123, -78, 34, 10, -5])


cdf31_htilde = (4, [-1, 3, 3, -1])
cdf33_htilde = (64, [3, -9, -7, 45, 45, -7, -9, 3])
cdf35_htilde = (512, [-5, 15, 19, -97, -26, 350, 350, -26, -97, 19, 15, -5])


cdf42_htilde = (32, [3, -12, 5, 40, 5, -12, 3])
cdf44_htilde = (512, [-10, 40, -2, -192, 140, 560, 140, -192, -2, 40, -10])
cdf46_htilde = (8192, [35, -140, -55, 920, -557, -2932, 2625, 8400, 2625, -2932, -557, 920, -55, -140, 35])


cdf51_htilde = (16, [3, -15, 20, 20, -15, 3])
cdf53_htilde = (128, [-5, 25, -26, 70, 140, 140, -70, -26, 25, -5])
cdf55_htilde = (4096, [35, -175, 120, 800, -1357, -1575, 4200, 4200, -1575, -1357, 800, 120, -175, 35])


cdf62_htilde = (64, [-5, 30, -56, -14, 154, -14, -56, 30, -5])
cdf64_htilde = (2048, [35, -210, 330, 470, -1827, 252, 3948, 252, -1827, 470, 330, -210, 35])
cdf66_htilde = (16384, [-63, 378, -476, -1554, 4404, 1114, -13860, 4158, 28182, 4158, -13860, 1114, 4404, -1554, -476, 378, -63])


# Determine the offset of the filter such that the filter is even symmetric w.r.t 
# 0 (for filters of odd length) or 1/2 (for filters of even length).
symmetric_offset(n) = iseven(n) ? -n>>1+1 : -(n-1)>>1

# The coefficients of the primal scaling function are simply a binomial sequence for all n.
cdf_coef(n, T) = sqrt(T(2))/T(2)^n*[binomial(n,k) for k in 0:n]


primal_scalingfilter{P,Q,T}(w::CDFWavelet{P,Q,T}) = CompactSequence(cdf_coef(P,T), symmetric_offset(P+1))

for (p,q,htilde) in ( (1, 1, :cdf11_htilde), (1, 3, :cdf13_htilde), (1, 5, :cdf15_htilde),
                      (2, 2, :cdf22_htilde), (2, 4, :cdf24_htilde), (2, 6, :cdf26_htilde),
                      (3, 1, :cdf31_htilde), (3, 3, :cdf33_htilde), (3, 5, :cdf35_htilde),
                      (4, 2, :cdf42_htilde), (4, 4, :cdf44_htilde), (4, 6, :cdf46_htilde),
                      (5, 1, :cdf51_htilde), (5, 3, :cdf53_htilde), (5, 5, :cdf55_htilde),
                      (6, 2, :cdf62_htilde), (6, 4, :cdf64_htilde), (6, 6, :cdf66_htilde)  )
    @eval dual_scalingfilter{T}(w::CDFWavelet{$p,$q,T}) = CompactSequence(sqrt(T(2))/$htilde[1]*$htilde[2], symmetric_offset(length($htilde[2])))
end




