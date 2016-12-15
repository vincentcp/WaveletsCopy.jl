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
        c = Wavelets.scaling_coefficients(x->one(T), DWT.primal_scalingfilter(w), PeriodicEmbedding(), L=L)
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
        for i in 0:length(c)-1
          @test abs(c[i+1]-2.0^(-3L/2)*i) < tol
        end
        c = Wavelets.scaling_coefficients(x->T(x), DWT.primal_scalingfilter(DWT.cdf22), nothing, L=L)
        for i in 0:length(c)-1
          @test abs(c[i+1]-2.0^(-3L/2)*i) < tol
        end
      end
    end
  end
  @testset "scaling coefficients of triangle function with cdf(2,x) - shifted" begin
    T = Float64
    tol = sqrt(eps(T))
    for w in (DWT.cdf22, DWT.cdf24, DWT.cdf26)
      for L in 4:10
        c = Wavelets.scaling_coefficients(x->T(1/2*(x+1)), DWT.cdf22, nothing, -1, 1, L=L)
        for i in 0:length(c)-1
          @test abs(c[i+1]-2*2.0^(-3L/2)*i) < tol
        end
      end
    end
  end
  @testset "order of spline wavelet should equal degree+1" begin
    L = 1
    for T in (Float16, Float32, Float64)
      println(string(T))
      tol = 1e3*eps(T)
      for k in 1:5
        implemented = Symbol(string("IMPLEMENTED_CDF_WAVELETS_",T))
        @eval begin
          for w in DWT.$implemented
            if DWT.primal_vanishingmoments(w) > $k
              c = Wavelets.scaling_coefficients(x->$T(x)^$k, DWT.primal_scalingfilter(w), nothing, L=$L)
              s = DWT.primal_support(w)
              nodes = [s[1]:s[2]...]/2
              if iseven(DWT.primal_vanishingmoments(w)) && isodd($k)
                c1_c=0; c1_c_e = 0
              else
                c1_c, c1_c_e = quadgk(x->sqrt($T(2))*Wavelets.evaluate_primal_scalingfunction(w, 2*x)*x^$k, nodes...)
              end
              nodes = nodes + .5
              c2_c, c2_c_e = quadgk(x->sqrt($T(2))*Wavelets.evaluate_primal_scalingfunction(w, 2*x-1)*x^$k, nodes...)
              @test (norm(c[1]-c1_c) < $tol)
              @test (norm(c[2]-c2_c) < $tol)
            end
          end
        end
      end
    end
  end
end

function supporttest()
  @testset "Support" begin
    for w in (DWT.IMPLEMENTED_CDF_WAVELETS..., DWT.IMPLEMENTED_DB_WAVELETS...)
      f = DWT.primal_scalingfilter(w)
      s = DWT.primal_support(w)
      @test ((firstindex(f), lastindex(f)) == s)
      f = DWT.dual_scalingfilter(w)
      s = DWT.dual_support(w)
      @test ((firstindex(f), lastindex(f)) == s)
    end
  end
end

function vanishing_moments_test()
  @testset "Vanishing moments" begin
    for (w,moments) in ((DWT.cdf11,(1,1)), (DWT.cdf24,(2,4)), (DWT.cdf51,(5,1)), (DWT.db4, (4,4)))
      @test DWT.primal_vanishingmoments(w) == moments[1]
      @test DWT.dual_vanishingmoments(w) == moments[2]
    end
  end
end

elementarypropsofsplinetest()
cascadetest()
primalfunctiontest()
scalingtest()
supporttest()
vanishing_moments_test()

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
