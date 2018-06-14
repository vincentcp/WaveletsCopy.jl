# wvlt_cdf.jl
using CardinalBSplines
"The Cohen-Daubechies-Feauveau family of compactly supported biorthogonal wavelets."
struct CDFWavelet{P,Q,T} <: DiscreteWavelet{T}
end

is_symmetric{P,Q,T}(::Type{CDFWavelet{P,Q,T}}) = True

is_biorthogonal{P,Q,T}(::Type{CDFWavelet{P,Q,T}}) = True

T0 = Float64

IMPLEMENTED_CDF_WAVELETS = []
for ET in (Float16,      Float32,      Float64,      BigFloat)
    implemented = Symbol(string("IMPLEMENTED_CDF_WAVELETS_",ET))
    @eval $implemented = []
    for N1 in 1:6
        iseven(N1) ? (N2I = 2:2:6) : (N2I = 1:2:6)
        for N2 in N2I
            if N1 == 2 && N2 == 2
                continue
            end
            fname =
            cdf = Symbol(string("cdf",N1,N2,"_",ET))
            T = CDFWavelet{N1,N2,ET}
            @eval begin
                $cdf = $T()
                class(::$T) = string($T)
                push!($implemented,$cdf)
                export $cdf
            end
            if ET == T0
                fname = string("cdf",N1,N2)
                cdf64 = Symbol(fname)
                @eval begin
                    $cdf64 = $cdf
                    export $cdf64
                end
            end
        end
    end
    @eval IMPLEMENTED_CDF_WAVELETS = [IMPLEMENTED_CDF_WAVELETS...,$implemented...]
end
name{P,Q,T}(::Type{CDFWavelet{P,Q,T}}) = string("cdf",P,Q,"_",T)
name{P,Q}(::Type{CDFWavelet{P,Q,Float64}}) = string("cdf",P,Q)

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
cdf53_htilde = (128, [-5, 25, -26, -70, 140, 140, -70, -26, 25, -5])
cdf55_htilde = (4096, [35, -175, 120, 800, -1357, -1575, 4200, 4200, -1575, -1357, 800, 120, -175, 35])


cdf62_htilde = (64, [-5, 30, -56, -14, 154, -14, -56, 30, -5])
cdf64_htilde = (2048, [35, -210, 330, 470, -1827, 252, 3948, 252, -1827, 470, 330, -210, 35])
cdf66_htilde = (16384, [-63, 378, -476, -1554, 4404, 1114, -13860, 4158, 28182, 4158, -13860, 1114, 4404, -1554, -476, 378, -63])


# Determine the offset of the filter such that the filter is even symmetric w.r.t
# 0 (for filters of odd length) or 1/2 (for filters of even length).
symmetric_offset(n) = -((n-1)>>1)

# The coefficients of the primal scaling function are simply a binomial sequence for all n.
cdf_coef(n) = 1//(1<<(n-1))*[binomial(n,k) for k in 0:n]

filter{P,Q,T}(side::Prl, kind::Scl, ::Type{CDFWavelet{P,Q,T}}) = CompactSequence(T(1)/sqrt(T(2))*convert(Array{T,1},cdf_coef(P)), symmetric_offset(P+1))
filter{P,Q,T}(side::Prl, kind::Cof, ::Type{CDFWavelet{P,Q,T}}) = CompactSequence(cdf_coef(P), symmetric_offset(P+1))

for (p,q,htilde) in ( (1, 1, :cdf11_htilde), (1, 3, :cdf13_htilde), (1, 5, :cdf15_htilde),
                      (2, 2, :cdf22_htilde), (2, 4, :cdf24_htilde), (2, 6, :cdf26_htilde),
                      (3, 1, :cdf31_htilde), (3, 3, :cdf33_htilde), (3, 5, :cdf35_htilde),
                      (4, 2, :cdf42_htilde), (4, 4, :cdf44_htilde), (4, 6, :cdf46_htilde),
                      (5, 1, :cdf51_htilde), (5, 3, :cdf53_htilde), (5, 5, :cdf55_htilde),
                      (6, 2, :cdf62_htilde), (6, 4, :cdf64_htilde), (6, 6, :cdf66_htilde)  )
    @eval filter{T}(side::Dul, kind::Cof, ::Type{CDFWavelet{$p,$q,T}}) = CompactSequence(2//$htilde[1]*$htilde[2], symmetric_offset(length($htilde[2])))
    @eval filter{T}(side::Dul, kind::Scl, ::Type{CDFWavelet{$p,$q,T}}) = CompactSequence(sqrt(T(2))/$htilde[1]*convert(Array{T,1}, $htilde[2]), symmetric_offset(length($htilde[2])))
    @eval support{T}(side::Dul, kind::Scl, ::Type{CDFWavelet{$p,$q,T}}) = (symmetric_offset(length($htilde[2])),symmetric_offset(length($htilde[2]))+length($htilde[2])-1)
    @eval support_length{T}(side::Dul, kind::Scl, ::Type{CDFWavelet{$p,$q,T}}) = length($htilde[2])-1
end

vanishingmoments{N1,N2,T}(::Prl, ::Type{CDFWavelet{N1,N2,T}}) = N1
vanishingmoments{N1,N2,T}(::Dul, ::Type{CDFWavelet{N1,N2,T}}) = N2
support{N1,N2,T}(::Prl, ::Scl, ::Type{CDFWavelet{N1,N2,T}}) = (symmetric_offset(N1+1), symmetric_offset(N1+1) + N1)
support_length(::Prl, ::Scl, ::Type{CDFWavelet{N1,N2,T}}) where {N1,N2,T} = N1

evaluate{N1,N2,T,S<:Real}(side::Prl, kind::Scl, w::CDFWavelet{N1,N2,T}, j::Int, k::Int, x::S; options...) =
      T(2)^(j/2)*evaluate_Bspline(N1-1, T(2)^j*x-T(k)-T(DWT.symmetric_offset(N1+1)), promote_type(T, eltype(x)))

evaluate{N1,N2,T,S<:Real}(side::Prl, kind::Wvl, w::CDFWavelet{N1,N2,T}, j::Int, k::Int, x::S; options...) =
      mother_relation(Prl(), w, j, k, x; options...)
# Periodic transformed ϕ,  ϕjk,  is the tranformed version of periodized spline with period 2^j
# This consturction of methods is necesarry to avoid conflicts with the methods in discretewavelets.jl
evaluate_periodic{N1,N2,T,S<:Real}(side::Prl, kind::Scl, w::CDFWavelet{N1,N2,T}, j::Int, k::Int, x::S; options...) =
      T(2)^(j/2)*evaluate_periodic_Bspline(N1-1, T(2)^j*x-T(k)-T(DWT.symmetric_offset(N1+1)), T(1<<j), T)

evaluate_periodic{N1,N2,T,S<:Real}(side::Prl, kind::Wvl, w::CDFWavelet{N1,N2,T}, j::Int, k::Int, x::S; options...) =
      mother_relation_periodic(Prl(), w, j, k, x; options...)

function mother_relation{T,S<:Real}(side::Side, w::DiscreteWavelet{T}, j::Int, k::Int, x::S; options...)
    flt = filter(side, Wvl(), w)
    res = T(0)
    for l in  firstindex(flt):lastindex(flt)
        res += flt[l]*evaluate(side, Scl(), w, j+1, 2k, x-l/(1<<(j+1)); options...)
    end
    res
end

function mother_relation_periodic{T,S<:Real}(side::Side, w::DiscreteWavelet{T}, j::Int, k::Int, x::S; options...)
    flt = filter(side, Wvl(), w)
    res = T(0)
    for l in firstindex(flt):lastindex(flt)
        res += flt[l]*evaluate_periodic(side, Scl(), w, j+1, 2k, x-l/(1<<(j+1)); options...)
    end
    res
end
