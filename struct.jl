mutable struct Point
    x::Float64
    y::Float64
    z::Float64
    Point() = new(0.0,0.0,0.0)
end


mutable struct Face
    x::Float64
    y::Float64
    z::Float64
    owner::UInt32
    neighbour::UInt32
    n̂::Vector{Float64}
    ΔS::Float64
    points::Vector{Point}
    varₗ::Vector{Float64}
    varᵣ::Vector{Float64}
end


mutable struct Cell
    x::Float64
    y::Float64
    z::Float64
    faces::Vector{UInt32}
    points::Vector{UInt32}
    Ω::Float64
    var::Vector{Float64}
    Qᵐ::Vector{Float64}
    Qⁿ::Vector{Float64}
    Qⁿ⁻¹::Vector{Float64}
    Cell() = new(0.0,0.0,0.0,[],[],0.0,[],[],[],[])
end
