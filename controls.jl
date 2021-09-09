mutable struct controls
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
    Lco::Float64
    Uco::Float64
    
    time::Float64
    realIter::UInt32
    pseudoIter::UInt32
    residual::Float64

    p::UInt32
    u::UInt32
    v::UInt32
    w::UInt32
    T::UInt32
    Y₁::UInt32
    ρ::UInt32
    Hₜ::UInt32
    c::UInt32
    Y₂::UInt32
    α₁::UInt32
    α₂::UInt32
    ∂ρ∂p::UInt32
    ∂ρ∂T::UInt32
    ∂ρ∂Y₁::UInt32
    ∂Hₜ∂p::UInt32
    ∂Hₜ∂T::UInt32
    ∂Hₜ∂Y₁::UInt32
    Δτ::UInt32
    Vᵣ::UInt32
end
