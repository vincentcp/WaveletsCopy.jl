# test_evaluation.jl
using Base.Test
using Wavelets

WTS = Wavelets

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
        g = Wavelets.cascade_algorithm(h,L=L)
        F = map(f,x)
        @test (norm(g-F))<tol
      end
    end
  end
end

function primalfunctiontest()
  @testset "primal function of cdf" begin
    T = Float64
    tol = sqrt(eps(T))
    for L in 1:5
      for N in 1:6
        w = DWT.CDFWavelet{N,N,T}()
        f = x->WTS.evaluate_primal_scalingfunction(w, x)
        x = Wavelets.dyadicpointsofcascade(w,L)
        g = Wavelets.cascade_algorithm(w,dual=false,L=L)
        F = map(f,x)
        @test (norm(g-F))<tol
      end
    end
  end
end

function scalingtest()
  @testset "scaling_coefficients of constant function" begin
    T = Float64
    tol = sqrt(eps(T))
    for w in (DWT.IMPLEMENTED_DB_WAVELETS..., DWT.IMPLEMENTED_CDF_WAVELETS...)
      for L in 4:10
        c = Wavelets.scaling_coefficients(x->one(T), w, PeriodicEmbedding(), L=L)
        for i in 1:length(c)
          @test (abs(c[i]-2.0^(-L/2)) < tol)
        end
      end
    end
  end
  @testset "scaling_coefficients of triangle function with cdf(2,x)" begin
    T = Float64
    tol = sqrt(eps(T))
    for w in (DWT.cdf22, DWT.cdf24, DWT.cdf26)
      for L in 4:10
        c = Wavelets.scaling_coefficients(x->T(x), DWT.cdf22, nothing, L=L)
        c_control = 2.0^(-3L/2)*[0:(1<<L-1)...]
        for i in 0:length(c)-1
          @test abs(c[i+1]-2.0^(-3L/2)*i) < tol
        end
      end
    end
  end
end

# elementarypropsofsplinetest()
# cascadetest()
# primalfunctiontest()
scalingtest()


# Plot Daubechies wavelets
# using Plots
# gr()
# plot()
# for i in 2:2:6
#   f = Symbol(string("db",i))
#   @eval h = DWT.primal_scalingfilter(DWT.$f).a
#   println(h)
#   x = Wavelets.dyadicpointsofcascade(h,10)
#   g = Wavelets.cascade_algorithm(h,L=10)
#   plot!(x,g)
# end
# plot!()

# Plot spline wavelets
# using Plots; gr(); plot()
# x = linspace(-5,5,1000)
# for i in 1:6
#   c = Symbol(string("cdf",i,i))
#   @eval f = x->Wavelets.evaluate_primal_scalingfunction(DWT.$c, x)
#   ff = map(f, x)
#   plot!(x, ff)
# end
# plot!()
