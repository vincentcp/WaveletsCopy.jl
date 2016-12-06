# test_evaluation.jl
using Base.Test
using Wavelets

WTS = Wavelets

@testset "constant B spline" begin
  for x in linspace(-1.3,1.3)
    s = Wavelets.evaluate_primal_scalingfunction(WT.db1,x)
    @test 0.0 <= x < 1. ? s==1. : s==0.
    t = Wavelets.evaluate_primal_scalingfunction(WT.haar,x)
    @test 0.0 <= x < 1. ? t==1. : t==0.
    t = Wavelets.evaluate_primal_scalingfunction(WT.cdf11,x)
    @test 0.0 <= x < 1. ? t==1. : t==0.
  end
end

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


@testset "Elementary properties" begin
  tol = 1e-8
  S = 20
  for N in 1:10
    f = x->WTS.evaluate_primal_scalingfunction(N, x, Float64)
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
