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
      for filter in (DWT.IMPLEMENTED_DB_WAVELETS..., DWT.IMPLEMENTED_CDF_WAVELETS...)
        fb = Filterbank(filter)
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
      for filter in (DWT.IMPLEMENTED_DB_WAVELETS..., DWT.IMPLEMENTED_CDF_WAVELETS...)
        fb = Filterbank(filter)
        for bound in (DWT.perbound, DWT.symbound)
          for L in 1:n
            y = DWT.dwt(x, fb, DWT.perbound, L)
            xx = DWT.idwt(y, fb, DWT.perbound, L)
            @test (norm(xx-x)) < 1e-10
          end
        end
      end
    end
  end
end

# w = DWT.cdf66
# j = 1; k = 0; d = 10;
# f(x) = cos(2*pi*x)
# for d in j:10
#   dx = 1/(1<<d)
#   coefs = DWT.scaling_coefficients(f, w, d, DWT.perbound)
#   ff, xx = DWT.scaling_coefficients_to_dyadic_grid(coefs, w, grid=true)
#   println(norm(ff-map(f,xx))*dx)
# end
#
# w = DWT.cdf15
# j = 1; k = 0; d = 10;
# g(x) = x < .5 ? x : 1-x
# for d in j:10
#   dx = 1/(1<<d)
#   coefs = DWT.scaling_coefficients(g, w, d, DWT.perbound)
#   ff, xx = DWT.scaling_coefficients_to_dyadic_grid(coefs, w, grid=true)
#   println(norm(ff-map(g,xx))*dx)
# end
