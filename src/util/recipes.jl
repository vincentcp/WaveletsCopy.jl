## Plotting recipe for a wavelet:

immutable Primal end
immutable Dual end

@recipe function f(w::DWT.DiscreteWavelet)
  title --> DWT.name(w)
  w, is_orthogonal(w)
end

@recipe function f(w::DWT.DiscreteWavelet, is_orthogonal::DWT.True; j=0::Int, k=0::Int, points=256::Int, periodic=false)
  @series begin
    label -->"scaling"
    if !periodic
      f,x = DWT.primal_scalingfunction_in_dyadic_points(w,j,k,ceil(Int,log2(points)); points=true)
    else
      f,x = DWT.periodic_primal_scalingfunction_in_dyadic_points(w,j,k,ceil(Int,log2(points)); points=true)
      f = [f...,f[1]]; x = [x...,1]
    end
    x, f
  end
  @series begin
    label -->"wavelet"
    if !periodic
      f,x = DWT.primal_waveletfunction_in_dyadic_points(w,j,k,ceil(Int,log2(points)); points=true)
    else
      f,x = DWT.periodic_primal_waveletfunction_in_dyadic_points(w,j,k,ceil(Int,log2(points)); points=true)
      f = [f...,f[1]]; x = [x...,1]
    end
    x, f
  end
end

@recipe function f(w::DWT.DiscreteWavelet, is_orthogonal::DWT.False; side = :primal)
  if side==:primal
    w, Primal()
  elseif side==:dual
    w, Dual()
  elseif side==:both
    w, Both()
  end
end

@recipe function f(w::DWT.DiscreteWavelet, side::Primal)
  w, DWT.True()
end

@recipe function f(w::DWT.DiscreteWavelet, side::Dual; j=0::Int, k=0::Int, points=256::Int, periodic=false)
  @series begin
    label -->"scaling"
    if !periodic
      f,x = DWT.dual_scalingfunction_in_dyadic_points(w,j,k,ceil(Int,log2(points)); points=true)
    else
      f,x = DWT.periodic_dual_scalingfunction_in_dyadic_points(w,j,k,ceil(Int,log2(points)); points=true)
      f = [f...,f[1]]; x = [x...,1]
    end
    x, f
  end
  @series begin
    label -->"wavelet"
    if !periodic
      f,x = DWT.dual_waveletfunction_in_dyadic_points(w,j,k,ceil(Int,log2(points)); points=true)
    else
      f,x = DWT.periodic_dual_waveletfunction_in_dyadic_points(w,j,k,ceil(Int,log2(points)); points=true)
      f = [f...,f[1]]; x = [x...,1]
    end
    x, f
  end
end

@recipe f(x::LinSpace, f::AbstractVector) = collect(x), f
