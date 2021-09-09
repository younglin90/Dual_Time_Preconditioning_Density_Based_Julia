function construct_A_matrix_explicit!(
    👉::controls, 
    A::Array{Float64},
    cells::Vector{mesh.Cell}
)

    for i in 1:length(cells)

        ρ = cells[i].var[👉.ρ]
        u = cells[i].var[👉.u]
        v = cells[i].var[👉.v]
        w = cells[i].var[👉.w]
        T = cells[i].var[👉.T]
        Y₁ = cells[i].var[👉.Y₁]
        Hₜ = cells[i].var[👉.Hₜ]
        ∂ρ∂p = cells[i].var[👉.∂ρ∂p]
        ∂ρ∂T = cells[i].var[👉.∂ρ∂T]
        ∂Hₜ∂p = cells[i].var[👉.∂Hₜ∂p]
        ∂Hₜ∂T = cells[i].var[👉.∂Hₜ∂T]
        ∂ρ∂Y₁ = cells[i].var[👉.∂ρ∂Y₁]
        ∂Hₜ∂Y₁ = cells[i].var[👉.∂Hₜ∂Y₁]

        T = zeros(Float64, 6, 6)
        T[1, 1] = ∂ρ∂p
        T[2, 1] = ∂ρ∂p*u
        T[3, 1] = ∂ρ∂p*v
        T[4, 1] = ∂ρ∂p*w
        T[5, 1] = ∂ρ∂p*Hₜ + ρ*∂Hₜ∂p - 1.0
        
        T[2, 2] = ρ
        T[5, 2] = ρ*u
        
        T[3, 3] = ρ
        T[5, 3] = ρ*v
        
        T[4, 4] = ρ
        T[5, 4] = ρ*w
        
        T[1, 5] = ∂ρ∂T
        T[2, 5] = ∂ρ∂T*u
        T[3, 5] = ∂ρ∂T*v
        T[4, 5] = ∂ρ∂T*w
        T[5, 5] = ∂ρ∂T*Hₜ + ρ*∂Hₜ∂T

        T[6, 1] = ∂ρ∂p*Y₁
        T[6, 5] = ∂ρ∂T*Y₁

        T[6, 6] = ∂ρ∂Y₁*Y₁
        T[6, 6] += ρ

        T[1, 6] = ∂ρ∂Y₁
        T[2, 6] = ∂ρ∂Y₁*u
        T[3, 6] = ∂ρ∂Y₁*v
        T[4, 6] = ∂ρ∂Y₁*w
        T[5, 6] = ∂ρ∂Y₁*Hₜ + ρ*∂Hₜ∂Y₁

        P = deepcopy(T)
        β = 1.0/cells[i].var[👉.Vᵣ]^2.0
            + ∂ρ∂p - 1.0/cells[i].var[👉.c]^2.0
        P[1, 1] = β
        P[2, 1] = β*u
        P[3, 1] = β*v
        P[4, 1] = β*w
        P[5, 1] = β*Hₜ + ρ*∂Hₜ∂p - 1.0
        P[6, 1] = β*Y₁

        A[i, :, :] = P[:, :]./cells[i].var[👉.Δτ] + 1.5*T[:, :]./👉.Δt


    end


end



function construct_A_matrix_implicit!(
    👉::controls, 
    A::Array{Float64},
    cells::Vector{mesh.Cell},
    faces::Vector{mesh.Face}
)




end