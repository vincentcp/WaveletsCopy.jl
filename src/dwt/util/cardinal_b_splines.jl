# cardinal_b_splines.jl
module Cardinal_b_splines

export evaluate_Bspline

# Implementation of cardinal B splines of degree N
typealias Degree{N} Val{N}

evaluate_Bspline(N::Int, x, T::Type) = evaluate_Bspline(Degree{N}, x, T)
function evaluate_Bspline{N}(::Type{Degree{N}}, x, T::Type)
  T(x)/T(N)*evaluate_Bspline(Degree{N-1}, x, T) +
      (T(N+1)-T(x))/T(N)*evaluate_Bspline(Degree{N-1}, x-1, T)
end

evaluate_Bspline(::Type{Degree{0}}, x, T::Type) = (0 <= x < 1) ? T(1) : T(0)

function evaluate_Bspline(::Type{Degree{1}}, x, T::Type)
  if (0 <= x < 1)
    s = T[0, 1]
  elseif (1 <= x < 2)
    s = T[2, -1]
  else
    s = T[0, 0]
  end
  @eval @evalpoly $x $(s...)
end
function evaluate_Bspline(::Type{Degree{2}}, x, T::Type)
  if (0 <= x < 1)
    s = T[0, 0, 1/2]
  elseif (1 <= x < 2)
    s = T[-3/2, 3, -1]
  elseif (2 <= x < 3)
    s = T[9/2, -3, 1/2]
  else
    s = T[0, 0]
  end
  @eval @evalpoly $x $(s...)
end

function evaluate_Bspline(::Type{Degree{3}}, x, T::Type)
  if (0 <= x < 1)
    s = T[0, 0, 0, 1/6]
  elseif (1 <= x < 2)
    s = T[2/3, -2, 2, -1/2]
  elseif (2 <= x < 3)
    s = T[-22/3, 10, -4, 1/2]
  elseif (3 <= x < 4)
    s = T[32/3, -8, 2, -1/6]
  else
    s = T[0, 0]
  end
  @eval @evalpoly $x $(s...)
end

function evaluate_Bspline(::Type{Degree{4}}, x, T::Type)
  if (0 <= x < 1)
    s = T[0, 0, 0, 0, 1/24]
  elseif (1 <= x < 2)
    s = T[-5/24, 5/6, -5/4, 5/6, -1/6]
  elseif (2 <= x < 3)
    s = T[155/24, -25/2, 35/4, -5/2, 1/4]
  elseif (3 <= x < 4)
    s = T[-655/24, 65/2, -55/4, 5/2, -1/6]
  elseif (4 <= x < 5)
    s = T[625/24, -125/6, 25/4, -5/6, 1/24]
  else
    s = T[0, 0]
  end
  @eval @evalpoly $x $(s...)
end

end # module Cardinal_b_splines
