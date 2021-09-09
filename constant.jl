mutable struct constant
    
    Nx::UInt32
    Ny::UInt32
    Nz::UInt32

    Lx::Float64
    Ly::Float64
    Lz::Float64
    
    realMaxIter::UInt32
    pseudoMaxIter::UInt32
    pseudoMaxResidual::Float64

    CFL::Float64
    
    Δt::Float64
end
