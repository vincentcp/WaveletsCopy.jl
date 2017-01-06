## Plotting recipe for a wavelet:

@recipe function f(w::DWT.DiscreteWavelet)
  w, is_orthogonal(w)
end
@recipe function f(side::DWT.Side, kind::DWT.Kind, w::DWT.DiscreteWavelet; j=0, k=0, points=256, periodic=false)
  label -->string(DWT.name(side), " ", DWT.name(kind), " of ", DWT.name(w))
  if !periodic
    f,x = DWT.eval_in_dyadic_points(side,kind,w,j,k,ceil(Int,log2(points)); points=true)
  else
    f,x = DWT.eval_periodic_in_dyadic_points(side,kind,w,j,k,ceil(Int,log2(points)); points=true)
    f = [f...,f[1]]; x = [x...,1]
  end
  x, f
end

@recipe f(w::DWT.DiscreteWavelet, is_orthogonal::DWT.True) = primal, w

@recipe function f(side::DWT.Side, w::DWT.DiscreteWavelet)
  @series begin
    side, scaling, w
  end
  @series begin
    side, DWT.wavelet, w
  end
end

immutable Both <: DWT.Side end

@recipe function f(w::DWT.DiscreteWavelet, is_orthogonal::DWT.False; side=Both()::DWT.Side)
  side, w
end

@recipe function f(::Both, w::DWT.DiscreteWavelet)
  @series begin
    primal, w
  end
  @series begin
    dual, w
  end
end

@recipe f(x::LinSpace, f::AbstractVector) = collect(x), f
