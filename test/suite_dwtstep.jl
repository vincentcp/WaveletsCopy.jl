# suit_dwtstep.jl
using WaveletsCopy
using Base.Test
W = WaveletsCopy
P = 80
W.DWT.perbound
W.DWT.symbound
W.DWT.zerobound
rng = MersenneTwister(3000)
jumpfunction(x) = (-0.5 < x < 0.5 ? 1.0: 0.0) + (-0.25 < x < .75 ? 1.0 : 0.0)
characteristicfunction(x) = (0<x<1) ? 1.0 : 0.0
randomfunction(x) = rand(rng)

@testset "$(rpad("Inversibility of dwtstep",P))"  begin
  for N in 10:10:1000
    t = linspace(-1,1,N)
    for f in (sin, characteristicfunction, randomfunction)
      x = map(f,t)
      for w in (W.DWT.IMPLEMENTED_DB_WAVELETS..., W.DWT.IMPLEMENTED_CDF_WAVELETS...)
        fb = W.Filterbank(w)
        for bound in (DWT.perbound, )
          y = dwtstep(x, fb, bound)
          xx = idwtstep(y..., fb, bound)
          @test (norm(xx-x)) < 1e-10
        end
      end
    end
  end
end

@testset "$(rpad("Inversibility of dwt",P))"  begin
  T = Float64
  for n in 1:10
    N = 2^n
    t = linspace(-1,1,N)
    for f in (sin, characteristicfunction, randomfunction)
      x = map(f,t)
      for w in (DWT.IMPLEMENTED_DB_WAVELETS..., DWT.IMPLEMENTED_CDF_WAVELETS...)
        for bound in (DWT.perbound,)
          for L in 1:n
            y = DWT.dwt(x, w, bound, L)
            xx = DWT.idwt(y, w, bound, L)
            @test (norm(xx-x)) < 1e-10
          end
        end
      end
    end
  end
end

@testset "$(rpad("Inversibility of full_dwt using constant function (periodic)",P))"  begin
  for l in 0:10
    x0 = ones(1<<l)
    for w in DWT.IMPLEMENTED_WAVELETS
      x1 = full_dwt(x0, w, DWT.perbound)
      y  = full_idwt(x1, w, DWT.perbound)
      @test (norm(x0-y)< 1e-10)
    end
  end
end

@testset "$(rpad("Inversibility of full_dwt (periodic)",P))"  begin
  for l in 0:7
    for p in 0:9
      x0 = x0 = [Float64(i^p)/(1<<l) for i in 1:1<<l]; x0/=sum(x0)
      for i in p:9
        w = DWT.DaubechiesWavelet{i+1,Float64}()
        offset = (max([support_length(side, kind, w) for side in (Primal, Dual) for kind in (scaling, DWT.wavelet)]...))
        x1 = full_dwt(x0, w, DWT.perbound)
        y  = full_idwt(x1, w, DWT.perbound)
        d = abs.(y-x0)
        @test (sum(d[offset+1:end-offset])<1e-10)
      end
    end
  end
end
