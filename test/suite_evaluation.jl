# test_evaluation.jl
using Base.Test
using Wavelets

WTS = Wavelets
function constantbsplinetest()
  @testset "constant B spline" begin
    for x in linspace(-1.3,1.3)
      for w in (:(WT.db1), :(WT.haar), :(WT.cdf11), :(DWT.db1), :(DWT.cdf11))
        @eval begin
          s = Wavelets.evaluate_primal_scalingfunction($w,$x)
          @test 0.0 <= $x < 1. ? s==1. : s==0.
        end
      end
    end
  end
end

function linearbsplintest()
  @testset "linear B spline" begin
    for x in linspace(-1.3,2.3)
      s = Wavelets.evaluate_primal_scalingfunction(WT.cdf22,x)
      if 0.0 <= x < 1.
        @test x≈s
      elseif 1. <= x < 2.
        @test -x+2≈s
      else
        @test 0.0≈s
      end
    end
  end
end

function elementarypropsofsplinetest()
  @testset "Elementary properties" begin
    tol = 1e-8
    S = 20
    for N in 1:10
      f = x->WTS.evaluate_Bspline(N-1, x, Float64)
      # Integral should be 1
      I,e = quadgk(f, 0, N, reltol = tol)
      @test I≈1 && abs(e) < tol
      # Infinite summation of shifted versions is 1
      xx = linspace(N-1, N, S)[1:end-1]
      g = zeros(xx)
      for k in 0:N-1
        g += map(x->f(x-k), xx)
      end
      @test (norm(g-1) < tol)
      # Two scale relation
      x = linspace(-1, N+1, S)
      g = zeros(x)
      for k in 0:N
        g += binomial(N,k)*map(x->f(2x-k), x)
      end
      g *= 2.0^(-N+1)
      G = map(f, x)
      @test (norm(g-G)) < tol
    end
  end
end

function cascadetest()
  @testset "cascade_algorithm" begin
    T = Float64
    tol = sqrt(eps(T))
    for L in 1:5
      for N in 1:8
        f = x->WTS.evaluate_Bspline(N-1, x, Float64)
        h = DWT.cdf_coef(N,T)
        x = Wavelets.dyadicpointsofcascade(h,L)
        g = Wavelets.cascade_algorithm(h,L)
        F = map(f,x)
        @test (norm(g-F))<tol
      end
    end
  end
end

constantbsplinetest()
linearbsplintest()
elementarypropsofsplinetest()
cascadetest()

# using Plots
# gr()
# plot()
# for i in 1:5
#   f = Symbol(string("db",i))
#   @eval h = DWT.primal_scalingfilter(DWT.$f).a
#   x = Wavelets.dyadicpointsofcascade(h)
#   g = Wavelets.cascade_algorithm(h)
#   plot!(x,g)
# end
# plot!()
#
