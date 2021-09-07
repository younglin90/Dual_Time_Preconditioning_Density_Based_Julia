module mesh

    points::Array{point}
    faces::Array{face}
    cells::Array{cell}

    struct point
        x::Float64
        y::Float64
        z::Float64
    end


    struct face
        points::Int32

        cellₗ::Int32
        cellᵣ::Int32

        pₗ::Float64
        uₗ::Float64
        vₗ::Float64
        wₗ::Float64
        Tₗ::Float64
        Yᵢₗ::Array{Float64}
        αᵢₗ::Array{Float64}
        
        pᵣ::Float64
        uᵣ::Float64
        vᵣ::Float64
        wᵣ::Float64
        Tᵣ::Float64
        Yᵢᵣ::Array{Float64}
        αᵢᵣ::Array{Float64}
    end

    struct cell
        points::Int32
        faces::Int32
        p::Float64
        u::Float64
        v::Float64
        w::Float64
        T::Float64
        Yᵢ::Array{Float64}
        αᵢ::Array{Float64}
        Qᵐ::Array{Float64}
        Qⁿ::Array{Float64}
        Qⁿ⁻¹::Array{Float64}
        Hₜ::Float64
        ρ::Float64
        c::Float64
        ∂ρ_∂p::Float64
        ∂Hₜ_∂p::Float64
        ∂ρ_∂T::Float64
        ∂Hₜ_∂T::Float64
        ∂ρ_∂Yᵢ::Array{Float64}
        ∂Hₜ_∂Yᵢ::Array{Float64}
    end

        
end
