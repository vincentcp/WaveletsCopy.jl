# test_evaluation.jl
using Base.Test
using Wavelets

WTS = Wavelets
function pad(s,i=80)
   "$(rpad(s,i))"
end
P = 80

function elementarypropsofsplinetest()
  @testset "$(rpad("Elementary properties",P))" begin
    tol = 1e-8
    S = 20
    for N in 1:10
      f = x->DWT.Cardinal_b_splines.evaluate_Bspline(N-1, x, Float64)
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
  @testset "$(rpad("cascade_algorithm",P))"  begin
    T = Float64
    tol = sqrt(eps(T))
    for L in 0:5
      for N in 1:6
        f = x->WTS.DWT.Cardinal_b_splines.evaluate_Bspline(N-1, x+(N>>1), Float64)
        h = DWT.CDFWavelet{N,N,T}()
        x = Wavelets.DWT.primal_scaling_dyadicpointsofcascade(h,L)
        g = Wavelets.DWT.cascade_algorithm(h,L,dual=false)
        F = map(f,x)
        @test (norm(g-F))<tol
      end
    end
  end
end

function primalfunctiontest()
  @testset "$(rpad("primal function of cdf",P))"  begin
    T = Float64
    tol = sqrt(eps(T))
    for L in 1:5
      for N in 1:6
        w = DWT.CDFWavelet{N,N,T}()
        f = x->DWT.evaluate_primal_scalingfunction(w, x)
        x = DWT.primal_scaling_dyadicpointsofcascade(w,L)
        g = DWT.cascade_algorithm(w, L, dual=false)
        F = map(f,x)
        @test (norm(g-F))<tol
      end
      w = DWT.HaarWavelet{T}()
      f = x->WTS.DWT.evaluate_primal_scalingfunction(w, x)
      x = DWT.dual_scaling_dyadicpointsofcascade(w,L)
      g = DWT.cascade_algorithm(w,L,dual=false)
      F = map(f,x)
      @test (norm(g-F))<tol
    end
  end
end

function scalingtest()
  @testset "$(rpad("scaling_coefficients of constant function",P))" begin
    T = Float64
    tol = sqrt(eps(T))
    for w in (DWT.IMPLEMENTED_DB_WAVELETS..., DWT.IMPLEMENTED_CDF_WAVELETS...)
      for L in 4:10
        c = Wavelets.DWT.scaling_coefficients(x->one(T), w, L, PeriodicEmbedding())
        for i in 1:length(c)
          @test (abs(c[i]-2.0^(-L/2)) < tol)
        end
        c = Wavelets.DWT.scaling_coefficients(x->one(T), DWT.primal_scalingfilter(w), L, PeriodicEmbedding())
        for i in 1:length(c)
          @test (abs(c[i]-2.0^(-L/2)) < tol)
        end
      end
    end
  end
  @testset "$(rpad("scaling_coefficients of triangle function with cdf(2,x)",P))"  begin
    T = Float64
    tol = sqrt(eps(T))
    for w in (DWT.cdf24, DWT.cdf26)
      for L in 4:10
        c = Wavelets.DWT.scaling_coefficients(x->T(x), w, L, nothing)
        for i in 0:length(c)-1
          @test abs(c[i+1]-2.0^(-3L/2)*i) < tol
        end
        c = Wavelets.DWT.scaling_coefficients(x->T(x), DWT.primal_scalingfilter(w), L, nothing)
        for i in 0:length(c)-1
          @test abs(c[i+1]-2.0^(-3L/2)*i) < tol
        end
      end
    end
  end
  @testset "$(rpad("scaling coefficients of triangle function with cdf(2,x) - shifted",P))"  begin
    T = Float64
    tol = sqrt(eps(T))
    for w in (DWT.cdf24, DWT.cdf26)
      for L in 4:10
        c = Wavelets.DWT.scaling_coefficients(x->T(1/2*(x+1)), w, L, nothing, -1, 1)
        for i in 0:length(c)-1
          @test abs(c[i+1]-2*2.0^(-3L/2)*i) < tol
        end
      end
    end
  end
  @testset "$(rpad("order of spline wavelet should equal degree+1",P))"  begin
    L = 1
    for T in (Float32, Float64)
      println(string(T))
      tol = 1e3*eps(T)
      for k in 1:5
        implemented = Symbol(string("IMPLEMENTED_CDF_WAVELETS_",T))
        @eval begin
          for w in DWT.$implemented
            if DWT.primal_vanishingmoments(w) > $k
              c = Wavelets.DWT.scaling_coefficients(x->$T(x)^$k, DWT.primal_scalingfilter(w), $L, nothing)
              s = DWT.primal_scalingsupport(w)
              nodes = [s[1]:s[2]...]/2
              if iseven(DWT.primal_vanishingmoments(w)) && isodd($k)
                c1_c=0; c1_c_e = 0
              else
                c1_c, c1_c_e = quadgk(x->sqrt($T(2))*Wavelets.DWT.evaluate_primal_scalingfunction(w, 2*x)*x^$k, nodes...)
              end
              nodes = nodes + .5
              c2_c, c2_c_e = quadgk(x->sqrt($T(2))*Wavelets.DWT.evaluate_primal_scalingfunction(w, 2*x-1)*x^$k, nodes...)
              @test (norm(c[1]-c1_c) < $tol)
              @test (norm(c[2]-c2_c) < $tol)
            end
          end
        end
      end
    end
  end
end

function vanishing_moments_test_dual()
  @testset "$(rpad("Vanishing moments of dual scaling functions",P))"  begin
    tol = sqrt(eps(Float64))
    for w in DWT.IMPLEMENTED_WAVELETS
      for d in 1:1
        p = DWT.dual_vanishingmoments(w)-1
        scaling_coefficients = DWT.scaling_coefficients(x->x^p, w, d, nothing)
        for k in 0:(1<<d)-1
          D = 10
          scaling_function, x = DWT.dual_scalingfunction_in_dyadic_points(w,d,k,D, points = true)
          estimate = sum(scaling_function.*map(x->x^p,x))/(1<<D)
          @test (norm(estimate-scaling_coefficients[k+1])<tol)
        end
      end
    end
  end
  @testset "$(rpad("Vanishing moments of primal scaling functions",P))"  begin
    tol = sqrt(eps(Float64))
    for w in DWT.IMPLEMENTED_WAVELETS
      for d in 1:1
        p = DWT.primal_vanishingmoments(w)-1
        scaling_coefficients = DWT.scaling_coefficients(x->x^p, w, d, nothing, dual=false)
        for k in 0:(1<<d)-1
          D = 10
          scaling_function, x = DWT.primal_scalingfunction_in_dyadic_points(w,d,k,D, points = true)
          estimate = sum(scaling_function.*map(x->x^p,x))/(1<<D)
          @test (norm(estimate-scaling_coefficients[k+1])<tol)
        end
      end
    end
  end
end

function supporttest()
  @testset "$(rpad("Support",P))"  begin
    for w in (DWT.IMPLEMENTED_CDF_WAVELETS..., DWT.IMPLEMENTED_DB_WAVELETS...)
      f = DWT.primal_scalingfilter(w)
      s = DWT.primal_scalingsupport(w)
      @test ((firstindex(f), lastindex(f)) == s)
      f = DWT.dual_scalingfilter(w)
      s = DWT.dual_scalingsupport(w)
      @test ((firstindex(f), lastindex(f)) == s)
    end
  end
end

function vanishing_moments_test()
  @testset "$(rpad("Vanishing moments",P))"  begin
    for (w,moments) in ((DWT.cdf11,(1,1)), (DWT.cdf24,(2,4)), (DWT.cdf51,(5,1)), (DWT.db4, (4,4)))
      @test DWT.primal_vanishingmoments(w) == moments[1]
      @test DWT.dual_vanishingmoments(w) == moments[2]
    end
  end
end

function filter_tests()
  @testset "$(rpad("Sum of filters",P))"  begin
    for w in DWT.IMPLEMENTED_WAVELETS
      @test sum(DWT.primal_scalingfilter(w)) ≈ sqrt(2)
      @test sum(DWT.dual_scalingfilter(w)) ≈ sqrt(2)
      @test sum(DWT.primal_coefficientfilter(w)) ≈ 2
      @test sum(DWT.dual_coefficientfilter(w)) ≈ 2
    end
  end
end

function implementation_test()
  @testset "$(rpad("Some simple tests",P))" begin
    @test DWT.is_biorthogonal(DWT.cdf11)==True()
    @test DWT.is_biorthogonal(DWT.db4)==True()
    @test DWT.is_orthogonal(DWT.db1)==True()
    @test DWT.is_semiorthogonal(DWT.db5)==True()
    @test DWT.is_symmetric(DWT.cdf11)==True()
    @test DWT.is_symmetric(DWT.cdf33)==True()
    @test DWT.is_symmetric(DWT.db1)==True()
    @test DWT.name(DWT.db1) == "db1"
    @test DWT.name(DWT.cdf11) == "cdf11"
    @test DWT.name(DWT.cdf11_Float16) == "cdf11_Float16"
    @test DWT.class(DWT.db1) == "Wavelets.DWT.DaubechiesWavelet{1,Float64}"
    @test DWT.class(DWT.cdf11) == "Wavelets.DWT.CDFWavelet{1,1,Float64}"
    @test DWT.primal_scaling_dyadicpointsofcascade(DWT.cdf13,1,1,3)==linspace(.5,1,5)
    @test DWT.primal_wavelet_dyadicpointsofcascade(DWT.cdf13,1,1,3)==linspace(0,1.5,13)
    @test DWT.dual_scaling_dyadicpointsofcascade(DWT.cdf13,1,1,3)==linspace(-.5,2,21)
    @test DWT.dual_wavelet_dyadicpointsofcascade(DWT.cdf13,1,1,3)==linspace(0,1.5,13)
    @test norm(DWT.dual_waveletfunction_in_dyadic_points(DWT.cdf24,1,3,4)-[0.0,8.3234e-6,0.000177566,-0.0014307,0.00378807,0.00350662,
                                                                  -0.0305216,0.0262683,0.0808122,-0.0458831,0.0762876,-0.206873,
                                                                    -0.61956,0.20897,0.306045,1.09198,2.39743,0.456957,-0.351988,
                                                                      -1.5335,-3.72494,-1.5335,-0.351988,0.456957,2.39743,1.09198,
                                                                        0.306045,0.20897,-0.61956,-0.206873,0.0762876,-0.0458831,
                                                                          0.0808122,0.0262683,-0.0305216,0.00350662,0.00378807,
                                                                            -0.0014307,0.000177566,8.3234e-6,0.0]) <= 1e-5
    @test norm(DWT.primal_waveletfunction_in_dyadic_points(DWT.cdf24,1,3,4)-[0.0,-0.0165728,-0.0331456,-0.0497184,-0.0662913,-0.0828641,
                                                                    -0.0994369,-0.11601,-0.132583,-0.0110485,0.110485,0.232019,0.353553,
                                                                      0.475087,0.596621,0.718155,0.839689,0.132583,-0.574524,-1.28163,
                                                                        -1.98874,-1.28163,-0.574524,0.132583,0.839689,0.718155,0.596621,
                                                                          0.475087,0.353553,0.232019,0.110485,-0.0110485,-0.132583,-0.11601,
                                                                            -0.0994369,-0.0828641,-0.0662913,-0.0497184,-0.0331456,-0.0165728,0.0]) <= 1e-5
    @test norm(DWT.primal_scalingfunction_in_dyadic_points(DWT.db1)[1:end-1]-
        DWT.periodic_primal_scalingfunction_in_dyadic_points(DWT.db1))<1e-14
    @test norm(DWT.primal_waveletfunction_in_dyadic_points(DWT.db1)[1:end-1]-
        DWT.periodic_primal_waveletfunction_in_dyadic_points(DWT.db1))<1e-14
    @test norm(DWT.dual_scalingfunction_in_dyadic_points(DWT.db1)[1:end-1]-
        DWT.periodic_dual_scalingfunction_in_dyadic_points(DWT.db1))<1e-14
    @test norm(DWT.dual_waveletfunction_in_dyadic_points(DWT.db1)[1:end-1]-
        DWT.periodic_dual_waveletfunction_in_dyadic_points(DWT.db1))<1e-14
  end
end

implementation_test()
elementarypropsofsplinetest()
cascadetest()
primalfunctiontest()
vanishing_moments_test_dual()
scalingtest()
supporttest()
vanishing_moments_test()
filter_tests()


# # Plot Daubechies wavelets
# using Plots
# gr()
# plot()
# for i in 2:2:6
#   f = Symbol(string("db",i))
#   @eval begin
#     x = Wavelets.DWT.primal_scaling_dyadicpointsofcascade(DWT.$f,10)
#     g = Wavelets.DWT.cascade_algorithm(DWT.$f,10)
#   end
#   plot!(x,g)
# end
# plot!()

# # Plot spline scaling functions
# using Plots; gr(); plot()
# x = linspace(-5,5,1000)
# for i in 1:6
#   i==2 ? c = Symbol(string("cdf",i,4)): c = Symbol(string("cdf",i,i))
#   @eval f = x->Wavelets.DWT.evaluate_primal_scalingfunction(DWT.$c, x)
#   ff = map(f, x)
#   plot!(x, ff)
# end
# plot!()


# # Plot spline wavelets
# using Plots; gr(legend=false); plot(layout=(6,6)); i=1;
# for p in 1:6
#   qs = 2:2:6
#   isodd(p) && (qs = 1:2:6 )
#   for q in qs
#     f, x = Wavelets.DWT.primal_waveletfunction_in_dyadic_points(DWT.CDFWavelet{p,q,Float64}(); points=true)
#     plot!(x, -f, subplot=i)
#     f, x = Wavelets.DWT.primal_scalingfunction_in_dyadic_points(DWT.CDFWavelet{p,q,Float64}(); points=true)
#     plot!(x, f, subplot=i)
#     i += 1
#   end
#   for q in qs
#     f, x = Wavelets.DWT.dual_waveletfunction_in_dyadic_points(DWT.CDFWavelet{p,q,Float64}(); points=true)
#     plot!(x, -f, subplot=i)
#     f, x = Wavelets.DWT.dual_scalingfunction_in_dyadic_points(DWT.CDFWavelet{p,q,Float64}(); points=true)
#     plot!(x, f, subplot=i)
#     i += 1
#   end
# end
# plot!()

# # Easy way to plot wavelets
# using Wavelets
# using Plots
# gr()
# w = DWT.CDFWavelet{2,2,Float64}()
# plot(w; side=:primal, j=0, k=-1, periodic=false)
# plot!(w; side=:primal, j=0, k=-1, periodic=true)
