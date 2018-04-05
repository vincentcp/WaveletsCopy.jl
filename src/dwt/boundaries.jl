# BOUNDARY TYPES

abstract type WaveletBoundary end

# Periodic boundary condition
struct PeriodicBoundary <: WaveletBoundary
end

# Symmetric extension
struct SymmetricBoundary <: WaveletBoundary
end

# zero padding
struct ZeropaddingBoundary <: WaveletBoundary
end

# constant padding
struct CPBoundary{T} <: WaveletBoundary
    constant    ::  T
end

perbound = DWT.PeriodicBoundary()
symbound = DWT.SymmetricBoundary()
zerobound = DWT.ZeropaddingBoundary()
