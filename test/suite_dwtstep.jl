# suit_dwtstep.jl
using Wavelets
using Base.Test
typealias W Wavelets
P = 80
DWT.perbound
DWT.symbound
DWT.zerobound
rng = MersenneTwister(3000)
jumpfunction(x) = (-0.5 < x < 0.5 ? 1.0: 0.0) + (-0.25 < x < .75 ? 1.0 : 0.0)
characteristicfunction(x) = (0<x<1) ? 1.0 : 0.0
randomfunction(x) = rand(rng)

@testset "$(rpad("Inversibility of dwtstep",P))"  begin
  for N in 10:10:1000
    t = linspace(-1,1,N)
    for f in (sin, characteristicfunction, randomfunction)
      x = map(f,t)
      for w in (DWT.IMPLEMENTED_DB_WAVELETS..., DWT.IMPLEMENTED_CDF_WAVELETS...)
        fb = Filterbank(w)
        for bound in (DWT.perbound, DWT.symbound)
          y = dwtstep(x, fb, DWT.perbound)
          xx = idwtstep(y..., fb, DWT.perbound)
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
        for bound in (DWT.perbound, DWT.symbound)
          for L in 1:n
            y = DWT.dwt(x, w, DWT.perbound, L)
            xx = DWT.idwt(y, w, DWT.perbound, L)
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
        offset = (max([support_length(side, kind, w) for side in (primal, dual) for kind in (scaling, DWT.wavelet)]...))
        x1 = full_dwt(x0, w, DWT.perbound)
        y  = full_idwt(x1, w, DWT.perbound)
        d = abs(y-x0)
        @test (sum(d[offset+1:end-offset])<1e-10)
      end
    end
  end
end
