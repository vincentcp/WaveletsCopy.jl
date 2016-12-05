# test_evaluation.jl
using Base.Test
using Wavelets

@testset "characteristic" begin
  for x in linspace(-1.3,1.3)
    s = Wavelets.evaluate_primal_scalingfunction(WT.db1,x)
    @test 0.0 <= x < 1. ? s==1. : s==0.
    t = Wavelets.evaluate_primal_scalingfunction(WT.haar,x)
    @test 0.0 <= x < 1. ? t==1. : t==0.
  end
end
