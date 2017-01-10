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

function periodicbsplinetest()
  @testset "$(rpad("periodic B splines",P))"  begin
  T = Float64
  for N in 0:4
    period = real(N+1)
    for x in linspace(0,period,10)[1:end-1]
      @test (DWT.Cardinal_b_splines.evaluate_periodic_Bspline(N, x, period, T) ≈ DWT.Cardinal_b_splines.evaluate_Bspline(N, x, T))
      for k in -2:2
        @test (DWT.Cardinal_b_splines.evaluate_periodic_Bspline(N, x, period, T) ≈ DWT.Cardinal_b_splines.evaluate_periodic_Bspline(N, x+k*period, period, T))
      end
    end
    period = real(N+1)/3
    for x in linspace(0,period,10)[1:end-1]
      for k in -2:2
        @test (DWT.Cardinal_b_splines.evaluate_periodic_Bspline(N, x, period, T) ≈ DWT.Cardinal_b_splines.evaluate_periodic_Bspline(N, x+k*period, period, T))
      end
    end
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
        x = Wavelets.DWT.dyadicpointsofcascade(primal, scaling, h,L)
        g = Wavelets.DWT.cascade_algorithm(primal, scaling, h,L)
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
        f = x->DWT.evaluate(primal, scaling, w, 0, 0, x)
        x = DWT.dyadicpointsofcascade(primal, scaling, w,L)
        g = DWT.cascade_algorithm(primal, scaling, w, L)
        F = map(f,x)
        @test (norm(g-F))<tol
      end
      w = DWT.HaarWavelet{T}()
      f = x->WTS.DWT.evaluate(primal, scaling, w, 0, 0, x)
      x = DWT.dyadicpointsofcascade(dual, scaling, w,L)
      g = DWT.cascade_algorithm(primal, scaling, w,L)
      F = map(f,x)
      @test (norm(g-F))<tol
    end
  end
  @testset "$(rpad("closest dyadic point",P))"  begin
    for c in -1:-1:-14
      xtol = 10.0^c
      for xtest in linspace(-1,2)
        d,k = Wavelets.DWT.closest_dyadic_point(xtest, xtol; dmax=200)
        @test abs(k/(1<<d)-xtest)<xtol
      end
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
        c = Wavelets.DWT.scaling_coefficients(x->one(T), DWT.filter(primal, scaling, w), L, PeriodicEmbedding())
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
        c = Wavelets.DWT.scaling_coefficients(x->T(x), DWT.filter(primal, scaling, w), L, nothing)
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
            if DWT.vanishingmoments(primal, w) > $k
              c = Wavelets.DWT.scaling_coefficients(x->$T(x)^$k, DWT.filter(primal, scaling, w), $L, nothing)
              s = DWT.support(primal, scaling, w)
              nodes = [s[1]:s[2]...]/2
              if iseven(DWT.vanishingmoments(primal, w)) && isodd($k)
                c1_c=0; c1_c_e = 0
              else
                c1_c, c1_c_e = quadgk(x->sqrt($T(2))*Wavelets.DWT.evaluate(primal, scaling, w, 0, 0, 2*x)*x^$k, nodes...)
              end
              nodes = nodes + .5
              c2_c, c2_c_e = quadgk(x->sqrt($T(2))*Wavelets.DWT.evaluate(primal, scaling, w, 0, 0, 2*x-1)*x^$k, nodes...)
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
        p = DWT.vanishingmoments(dual, w)-1
        scaling_coefficients = DWT.scaling_coefficients(x->x^p, w, d, nothing)
        for k in 0:(1<<d)-1
          D = 10
          scaling_function, x = DWT.evaluate_in_dyadic_points(dual, scaling, w, d, k, D, points = true)
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
        p = DWT.vanishingmoments(primal, w)-1
        scaling_coefficients = DWT.scaling_coefficients(x->x^p, w, d, nothing, side=primal)
        for k in 0:(1<<d)-1
          D = 10
          scaling_function, x = DWT.evaluate_in_dyadic_points(primal, scaling, w, d, k, D, points = true)
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
      f = DWT.filter(primal, scaling, w)
      s = DWT.support(primal, scaling, w)
      @test ((firstindex(f), lastindex(f)) == s)
      f = DWT.filter(dual, scaling, w)
      s = DWT.support(dual, scaling, w)
      @test ((firstindex(f), lastindex(f)) == s)
    end
  end
end

function vanishing_moments_test()
  @testset "$(rpad("Vanishing moments",P))"  begin
    for (w,moments) in ((DWT.cdf11,(1,1)), (DWT.CDFWavelet{1,1,Float16}, (1,1)), (DWT.cdf24,(2,4)), (DWT.cdf51,(5,1)), (DWT.db4, (4,4)), (DWT.DaubechiesWavelet{1,Float16}, (1,1)))
      @test DWT.vanishingmoments(primal, w) == moments[1]
      @test DWT.vanishingmoments(dual, w) == moments[2]
    end
  end
end

function filter_tests()
  @testset "$(rpad("Sum of filters",P))"  begin
    for w in DWT.IMPLEMENTED_WAVELETS
      @test sum(DWT.filter(primal, scaling, w)) ≈ sqrt(2)
      @test sum(DWT.filter(dual, scaling, w)) ≈ sqrt(2)
      @test sum(DWT.filter(primal, coefficient, w)) ≈ 2
      @test sum(DWT.filter(primal, coefficient, w)) ≈ 2
    end
  end
end

function implementation_test()
  @testset "$(rpad("Some simple tests",P))" begin
    @test DWT.name(DWT.scaling) == "scaling"
    @test DWT.name(DWT.wavelet) == "wavelet"
    @test DWT.name(DWT.primal) == "primal"
    @test DWT.name(DWT.dual) == "dual"
    @test DWT.is_symmetric(DWT.TestWavelet{Float16}) == False
    @test DWT.is_biorthogonal(DWT.TestWavelet{Float16}) == False
    @test DWT.is_orthogonal(DWT.TestWavelet{Float16}) == False
    @test DWT.is_semiorthogonal(DWT.TestWavelet{Float16}) == False
    @test eltype(DWT.TestWavelet{Float16}) == Float16
    @test_throws String vanishingmoments(primal, DWT.TestWavelet)
    @test DWT.is_biorthogonal(DWT.cdf11)==True()
    @test DWT.is_biorthogonal(DWT.db4)==True()
    @test DWT.is_orthogonal(DWT.db1)==True()
    @test DWT.is_semiorthogonal(DWT.db5)==True()
    @test DWT.is_symmetric(DWT.cdf11)==True()
    @test DWT.is_symmetric(DWT.cdf33)==True()
    @test DWT.is_symmetric(DWT.db1)==True()
    @test DWT.is_biorthogonal(DWT.CDFWavelet{1,1,Float16})==True
    @test DWT.is_biorthogonal(DWT.DaubechiesWavelet{1,Float64})==True
    @test DWT.is_orthogonal(DWT.DaubechiesWavelet{1,Float16})==True
    @test DWT.is_semiorthogonal(DWT.DaubechiesWavelet{1,Float64})==True
    @test DWT.is_symmetric(DWT.CDFWavelet{1,1,Float16})==True
    @test DWT.is_symmetric(DWT.CDFWavelet{3,3,Float16})==True
    @test DWT.is_symmetric(DWT.DaubechiesWavelet{1,Float16})==True
    @test DWT.name(DWT.db1) == "db1"
    @test DWT.name(DWT.DaubechiesWavelet{1,Float16}()) == "db1_Float16"
    @test DWT.name(DWT.cdf11) == "cdf11"
    @test DWT.name(DWT.cdf11_Float16) == "cdf11_Float16"
    @test DWT.class(DWT.db1) == "Wavelets.DWT.DaubechiesWavelet{1,Float64}"
    @test DWT.class(DWT.cdf11) == "Wavelets.DWT.CDFWavelet{1,1,Float64}"
    @test DWT.dyadicpointsofcascade(primal, scaling, DWT.cdf11,1,0,0)==[0.]
    @test DWT.dyadicpointsofcascade(primal, scaling, DWT.cdf13,1,1,3)==linspace(.5,1,5)
    @test DWT.dyadicpointsofcascade(primal, DWT.wavelet, DWT.cdf13,1,1,3)==linspace(0,1.5,13)
    @test DWT.dyadicpointsofcascade(dual, scaling, DWT.cdf13,1,1,3)==linspace(-.5,2,21)
    @test DWT.dyadicpointsofcascade(dual, DWT.wavelet, DWT.cdf13,1,1,3)==linspace(0,1.5,13)
    @test norm(DWT.evaluate_in_dyadic_points(dual, DWT.wavelet, DWT.cdf24,1,3,4)-[0.0,8.3234e-6,0.000177566,-0.0014307,0.00378807,0.00350662,
                                                                  -0.0305216,0.0262683,0.0808122,-0.0458831,0.0762876,-0.206873,
                                                                    -0.61956,0.20897,0.306045,1.09198,2.39743,0.456957,-0.351988,
                                                                      -1.5335,-3.72494,-1.5335,-0.351988,0.456957,2.39743,1.09198,
                                                                        0.306045,0.20897,-0.61956,-0.206873,0.0762876,-0.0458831,
                                                                          0.0808122,0.0262683,-0.0305216,0.00350662,0.00378807,
                                                                            -0.0014307,0.000177566,8.3234e-6,0.0]) <= 1e-5
    @test norm(DWT.evaluate_in_dyadic_points(primal, DWT.wavelet, DWT.cdf24,1,3,4)-[0.0,-0.0165728,-0.0331456,-0.0497184,-0.0662913,-0.0828641,
                                                                    -0.0994369,-0.11601,-0.132583,-0.0110485,0.110485,0.232019,0.353553,
                                                                      0.475087,0.596621,0.718155,0.839689,0.132583,-0.574524,-1.28163,
                                                                        -1.98874,-1.28163,-0.574524,0.132583,0.839689,0.718155,0.596621,
                                                                          0.475087,0.353553,0.232019,0.110485,-0.0110485,-0.132583,-0.11601,
                                                                            -0.0994369,-0.0828641,-0.0662913,-0.0497184,-0.0331456,-0.0165728,0.0]) <= 1e-5
    @test norm(DWT.evaluate_in_dyadic_points(primal, scaling, DWT.db1)[1:end-1]-
        DWT.evaluate_periodic_in_dyadic_points(primal, scaling, DWT.db1))<1e-14
    @test norm(DWT.evaluate_in_dyadic_points(primal, DWT.wavelet, DWT.db1)[1:end-1]-
        DWT.evaluate_periodic_in_dyadic_points(primal, DWT.wavelet, DWT.db1))<1e-14
    @test norm(DWT.evaluate_in_dyadic_points(dual, scaling, DWT.db1)[1:end-1]-
        DWT.evaluate_periodic_in_dyadic_points(dual, scaling, DWT.db1))<1e-14
    @test norm(DWT.evaluate_in_dyadic_points(dual, DWT.wavelet, DWT.db1)[1:end-1]-
        DWT.evaluate_periodic_in_dyadic_points(dual, DWT.wavelet, DWT.db1))<1e-14

    @test norm(DWT.evaluate_in_dyadic_points(primal, scaling, DWT.cdf11)[1:end-1]-
        DWT.evaluate_periodic_in_dyadic_points(primal, scaling, DWT.cdf11))<1e-14
    @test norm(DWT.evaluate_in_dyadic_points(primal, DWT.wavelet, DWT.cdf11)[1:end-1]-
        DWT.evaluate_periodic_in_dyadic_points(primal, DWT.wavelet, DWT.cdf11))<1e-14
    @test norm(DWT.evaluate_in_dyadic_points(dual, scaling, DWT.cdf11)[1:end-1]-
        DWT.evaluate_periodic_in_dyadic_points(dual, scaling, DWT.cdf11))<1e-14
    @test norm(DWT.evaluate_in_dyadic_points(dual, DWT.wavelet, DWT.cdf11)[1:end-1]-
        DWT.evaluate_periodic_in_dyadic_points(dual, DWT.wavelet, DWT.cdf11))<1e-14

    @test norm(DWT.evaluate_in_dyadic_points(primal, scaling, DWT.db1)-
        DWT.evaluate_in_dyadic_points(primal, scaling, DWT.cdf11))<1e-14
    @test norm(DWT.evaluate_in_dyadic_points(primal, DWT.wavelet, DWT.db1)-
        DWT.evaluate_in_dyadic_points(primal, DWT.wavelet, DWT.cdf11))<1e-14
    @test norm(DWT.evaluate_in_dyadic_points(dual, scaling, DWT.db1)-
        DWT.evaluate_in_dyadic_points(dual, scaling, DWT.cdf11))<1e-14
    @test norm(DWT.evaluate_in_dyadic_points(dual, DWT.wavelet, DWT.db1)-
        DWT.evaluate_in_dyadic_points(dual, DWT.wavelet, DWT.cdf11))<1e-14

    @test norm(DWT.evaluate_periodic_in_dyadic_points(primal, scaling, DWT.db1)-
        DWT.evaluate_periodic_in_dyadic_points(primal, scaling, DWT.cdf11))<1e-14
    @test norm(DWT.evaluate_periodic_in_dyadic_points(primal, DWT.wavelet, DWT.db1)-
        DWT.evaluate_periodic_in_dyadic_points(primal, DWT.wavelet, DWT.cdf11))<1e-14
    @test norm(DWT.evaluate_periodic_in_dyadic_points(dual, scaling, DWT.db1)-
        DWT.evaluate_periodic_in_dyadic_points(dual, scaling, DWT.cdf11))<1e-14
    @test norm(DWT.evaluate_periodic_in_dyadic_points(dual, DWT.wavelet, DWT.db1)-
        DWT.evaluate_periodic_in_dyadic_points(dual, DWT.wavelet, DWT.cdf11))<1e-14
  end
end

function eval_wavelet_test()
  @testset "$(rpad("evaluate wavelet in point",P))" begin
    @test DWT._periodize((.25,.5),-1,1)==((.25,.5),)
    @test DWT._periodize((-1.25,1.5),-1,1)==((-1,1),)
    @test DWT._periodize((-1.25,.5),-1,1)==((.75,1.),(-1.,.5))
    @test DWT._periodize((-1.5,0.),-1,1)==((0.5,1.),(-1,0))
    for side in (primal, dual)
      for kind in (scaling, DWT.wavelet)
        @test DWT.in_periodic_support(1,DWT.periodic_support(side,kind, DWT.cdf11, 0,0)...)
        @test DWT.in_periodic_support(0,DWT.periodic_support(side,kind, DWT.cdf11, 3,0)...)
        @test !DWT.in_periodic_support(0,DWT.periodic_support(side,kind, DWT.cdf11, 1,1)...)
      end
    end
  end
end

implementation_test()
periodicbsplinetest()
elementarypropsofsplinetest()
cascadetest()
primalfunctiontest()
scalingtest()
supporttest()
vanishing_moments_test()
filter_tests()
eval_wavelet_test()
vanishing_moments_test_dual()

# # Plot Daubechies wavelets
# using Plots
# gr()
# plot(layout=(2,2))
# for i in 1:2:20
#   plot!(primal, DWT.scaling, DWT.DaubechiesWavelet{i,Float64}(),subplot=1)
# end
# plot!(xlims=[0,20],subplot=1)
# for i in 1:2:20
#   plot!(primal, DWT.wavelet, DWT.DaubechiesWavelet{i,Float64}(),subplot=2)
# end
# plot!(xlims=[-5,15],subplot=2)
# for i in 1:2:20
#   plot!(primal, DWT.scaling, DWT.DaubechiesWavelet{i,Float64}(),periodic=true,subplot=3)
# end
# plot!(xlims=[0,2],subplot=3)
# for i in 1:2:20
#   plot!(primal, DWT.wavelet, DWT.DaubechiesWavelet{i,Float64}(),periodic=true,subplot=4)
# end
# plot!(xlims=[0,2],subplot=4)


# # Plot spline scaling functions
# using Plots; gr(); plot()
# x = linspace(-5,5,1000)
# for i in 1:6
#   f = x->Wavelets.DWT.evaluate(primal, scaling, DWT.CDFWavelet{i,i,Float64}(), 0,0,x)
#   ff = map(f, x)
#   plot!(x, ff)
# end
# plot!()


# # Plot spline wavelets
# using Plots; gr(legend=false); plot(layout=(6,6)); i=1;j=0;k=0;d=6
# for p in 1:6
#   qs = 2:2:6
#   isodd(p) && (qs = 1:2:6 )
#   for q in qs
#     println(p,q)
#     f, x = Wavelets.DWT.evaluate_in_dyadic_points(primal, DWT.wavelet, DWT.CDFWavelet{p,q,Float64}(), j, k, d; points=true)
#     plot!(x, -f, subplot=i)
#     f, x = Wavelets.DWT.evaluate_in_dyadic_points(primal, scaling, DWT.CDFWavelet{p,q,Float64}(), j, k, d; points=true)
#     plot!(x, f, subplot=i)
#     i += 1
#   end
#   for q in qs
#     println(p,q,"dual")
#     # f, x = Wavelets.DWT.evaluate_in_dyadic_points(dual, DWT.wavelet, DWT.CDFWavelet{p,q,Float64}(), j, k, d; points=true)
#     # plot!(x, -f, subplot=i)
#     f, x = Wavelets.DWT.evaluate_in_dyadic_points(dual, scaling, DWT.CDFWavelet{p,q,Float64}(), j, k, d; points=true)
#     plot!(x, f, subplot=i)
#     i += 1
#   end
# end
# plot!()

# using Plots
# using Wavelets
# gr()
# w = DWT.CDFWavelet{1,3,Float64}()
# plot(w; j=0, k=-1, periodic=false)
# plot(w; side=dual, j=0, k=-1, periodic=true)

# using Plots
# gr()
# f = x-> DWT.evaluate_periodic_Bspline(3, x, 0.75, Float64)
# x = linspace(-1,2,200)
# ff = map(f, x)
# plot(x, ff)
