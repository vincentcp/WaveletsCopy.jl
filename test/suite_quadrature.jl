if VERSION < v"0.7-"
    using Base.Test
else
    using Test, LinearAlgebra
    linspace(a,b,c) = range(a, stop=b, length=c)
end

using WaveletsCopy.DWT: quad_trap, quad_sf, quad_sf_weights, quad_sf_N, quad_trap_N, db4, db3, scaling, Primal
function test_wavelet_quadrature()
    @testset begin
        @test quad_trap(x->x, db4, 0, 0, 8)^2≈quad_trap(x->x^2, db4, 0, 0, 8)
        @test quad_sf(x->x, db3, 5, 0, 0,8)^2≈quad_sf(x->x^2, db3, 5, 0, 0, 8)

        M1 = 5; M2 = 10; J = 0
        f = sin ; wav = db3;
        reference = .741104421925905
        @test abs(quad_sf(x->f(2x-1), wav, M1, J+1, 1, 6) - reference) < 1e-3
        @test abs(quad_sf(x->f(2x-1), wav, M2, J+1, 1, 3) - reference) < 1e-15
        @test abs(quad_trap(f, wav, J, 0, 13) - reference) < 1e-12
        @test abs(dot(quad_sf_weights(Primal, scaling, wav, M1, 7), f.(linspace(0,5,M1*2^7+1)))-reference) < 1e-14
        @test abs(dot(quad_sf_weights(Primal, scaling, wav, M2, 3), f.(linspace(0,5,M2*2^3+1)))-reference) < 1e-15

        g = x->sin(2pi*x);J = 3
        reference = quad_sf(g, wav, M2, J, 0, 5; periodic=true)
        @test abs(quad_sf(g, wav, M1, J, 0, 6; periodic=true)-reference) < 1e-14
        @test abs(quad_sf(x->g(x-7//8), wav, M2, J, 7, 4; periodic=true)-reference)<1e-15
        @test abs(quad_trap(g, wav, J, 0, 12; periodic=true) -reference) < 1e-9

        reference = [quad_sf(g, wav, M2, J, k, 5; periodic=true) for k in 0:7]
        @test norm(quad_sf_N(g, wav, M2, J, 5)-reference) < 1e-15
        @test norm(quad_sf_N(g, wav, M1, J, 7)-reference) < 1e-14
        @test norm(quad_trap_N(g, wav, J, 12)-reference) < 1e-9

    end
end

test_wavelet_quadrature()
