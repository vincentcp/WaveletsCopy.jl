# suit_dwtstep.jl
using Wavelets
using Base.Test
typealias W Wavelets

DWT.perbound
DWT.symbound
DWT.zerobound

jumpfunction(x) = (-0.5 < x < 0.5 ? 1 : 0) + (-0.25 < x < .75 ? 1 : 0)
characteristicfunction(x) = (0<x<1) ? 1 : 0
randomfunction(x) = rand()

@testset "Inversibility of dwtstep" begin
  for N in 10:10:1000
    t = linspace(-1,1,N)
    for f in (sin, characteristicfunction, randomfunction)
      x = map(f,t)
      for filter in (DWT.IMPLEMENTED_DB_WAVELETS..., DWT.IMPLEMENTED_CDF_WAVELETS...)
        fb = Filterbank(filter)
        for bound in (DWT.perbound, DWT.symbound)
          y = dwtstep(x, fb, DWT.perbound)
          xx = idwtstep(y...,fb,DWT.perbound)
          @test (norm(xx-x)) < 1e-10
        end
      end
    end
  end
end
# exit()
# x = rand(9)
# println(x')
# y = dwtstep(x, Filterbank(DWT.db1), DWT.perbound)
# #println(y)
# xx = idwtstep(y...,Filterbank(DWT.db1),DWT.perbound)
# println(xx')
# println((norm(xx-x)))
#
