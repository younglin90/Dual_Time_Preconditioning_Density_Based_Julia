function construct_A_matrix_explicit!(
    👉::controls, 
    cells::Vector{mesh.Cell},
    faces::Vector{mesh.Face},
    faces_internal::Vector{mesh.Face},
    faces_boundary::Vector{mesh.Face},
    faces_boundary_top::Vector{mesh.Face},
    faces_boundary_bottom::Vector{mesh.Face},
    faces_boundary_left::Vector{mesh.Face},
    faces_boundary_right::Vector{mesh.Face},
    B::Array{Float64}
)

    x = zeros(Float64,length(cells),6)

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

        A = zeros(Float64,6,6)
        A[:, :] = P[:, :]./cells[i].var[👉.Δτ] .*cells[i].Ω + 1.5*T[:, :]./👉.Δt .*cells[i].Ω

        x[i, :] = A[:, :]\B[i, :]
        
        cells[i].var[👉.p] += x[i, 1]
        cells[i].var[👉.u] += x[i, 2]
        cells[i].var[👉.v] += x[i, 3]
        cells[i].var[👉.w] += x[i, 4]
        cells[i].var[👉.T] += x[i, 5]
        cells[i].var[👉.Y₁] += x[i, 6]

        cells[i].var[👉.p] = max(cells[i].var[👉.p],1.e-200)
        cells[i].var[👉.T] = max(cells[i].var[👉.T],1.e-200)

    end

    

    return norm(x)


end



function construct_A_matrix_implicit!(
    👉::controls, 
    cells::Vector{mesh.Cell},
    faces::Vector{mesh.Face},
    faces_internal::Vector{mesh.Face},
    faces_boundary::Vector{mesh.Face},
    faces_boundary_top::Vector{mesh.Face},
    faces_boundary_bottom::Vector{mesh.Face},
    faces_boundary_left::Vector{mesh.Face},
    faces_boundary_right::Vector{mesh.Face},
    B_inp::Array{Float64}
)
    
    B_n = 5
    A_n = B_n * B_n

    A_rows = zeros(Int64, length(cells)*A_n)
    A_cols = zeros(Int64, length(cells)*A_n)
    A_vals = zeros(Float64, length(cells)*A_n)
    B = zeros(Float64, length(cells)*B_n)


    diagon = 1
    
    for cell in cells

        ijStart = B_n*(diagon-1)
        Astart = A_n*(diagon-1)
        i = Astart

        Ω = cell.Ω
        Δt = 👉.Δt
        
        ρ = cells[diagon].var[👉.ρ]
        u = cells[diagon].var[👉.u]
        v = cells[diagon].var[👉.v]
        w = cells[diagon].var[👉.w]
        T = cells[diagon].var[👉.T]
        Y₁ = cells[diagon].var[👉.Y₁]
        Hₜ = cells[diagon].var[👉.Hₜ]
        ∂ρ∂p = cells[diagon].var[👉.∂ρ∂p]
        ∂ρ∂T = cells[diagon].var[👉.∂ρ∂T]
        ∂Hₜ∂p = cells[diagon].var[👉.∂Hₜ∂p]
        ∂Hₜ∂T = cells[diagon].var[👉.∂Hₜ∂T]
        ∂ρ∂Y₁ = cells[diagon].var[👉.∂ρ∂Y₁]
        ∂Hₜ∂Y₁ = cells[diagon].var[👉.∂Hₜ∂Y₁]

        T = zeros(Float64, B_n, B_n)
        T[1, 1] = ∂ρ∂p
        T[2, 1] = ∂ρ∂p*u
        T[3, 1] = ∂ρ∂p*v
        T[4, 1] = ∂ρ∂p*Hₜ + ρ*∂Hₜ∂p - 1.0
        
        T[2, 2] = ρ
        T[4, 2] = ρ*u
        
        T[3, 3] = ρ
        T[4, 3] = ρ*v
        
        T[1, 4] = ∂ρ∂T
        T[2, 4] = ∂ρ∂T*u
        T[3, 4] = ∂ρ∂T*v
        T[4, 4] = ∂ρ∂T*Hₜ + ρ*∂Hₜ∂T

        T[5, 1] = ∂ρ∂p*Y₁
        T[5, 4] = ∂ρ∂T*Y₁

        T[5, 5] = ∂ρ∂Y₁*Y₁
        T[5, 5] += ρ

        T[1, 5] = ∂ρ∂Y₁
        T[2, 5] = ∂ρ∂Y₁*u
        T[3, 5] = ∂ρ∂Y₁*v
        T[4, 5] = ∂ρ∂Y₁*Hₜ + ρ*∂Hₜ∂Y₁

        P = deepcopy(T)
        β = 1.0/cells[diagon].var[👉.Vᵣ]^2.0
            + ∂ρ∂p - 1.0/cells[diagon].var[👉.c]^2.0
        P[1, 1] = β
        P[2, 1] = β*u
        P[3, 1] = β*v
        P[4, 1] = β*Hₜ + ρ*∂Hₜ∂p - 1.0
        P[5, 1] = β*Y₁

        
        A = zeros(Float64, B_n , B_n)
        for i_ in 1:B_n
            for j_ in 1:B_n
                A[i_, j_] = P[i_, j_]/cells[diagon].var[👉.Δτ]*Ω + 1.5*T[i_, j_]/👉.Δt*Ω
            end
        end
        #A[:, :] = P[:, :]./cells[diagon].var[👉.Δτ] + 1.5*T[:, :]./👉.Δt
        #A = A.*Ω

        A[3, 1] += (-∂ρ∂p*(-9.8)*Ω)
        A[3, 4] += (-∂ρ∂T*(-9.8)*Ω)
        A[3, 5] += (-∂ρ∂Y₁*(-9.8)*Ω)


        #println(pⁿ,uⁿ,vⁿ,Hₜⁿ)
        # continuity
        for i_ in 1:B_n
            for j_ in 1:B_n
                i += 1
                A_rows[i] = ijStart + i_; A_cols[i] = ijStart + j_
                A_vals[i] = A[i_, j_]
            end
        end

        B[ijStart + 1] = B_inp[diagon, 1]
        B[ijStart + 2] = B_inp[diagon, 2]
        B[ijStart + 3] = B_inp[diagon, 3]
        B[ijStart + 4] = B_inp[diagon, 5]
        B[ijStart + 5] = B_inp[diagon, 6]


        diagon += 1


    end




        
    ∂Δp∂x0 = zeros(Float64, length(cells), 3)
    for face in faces_internal
        pₙ = 0.5 * (cells[face.owner].var[👉.p] + cells[face.neighbour].var[👉.p])
        ∂Δp∂x0[face.owner, 1] += pₙ * face.n̂[1] * face.ΔS / cells[face.owner].Ω
        ∂Δp∂x0[face.owner, 2] += pₙ * face.n̂[2] * face.ΔS / cells[face.owner].Ω
        ∂Δp∂x0[face.owner, 3] += pₙ * face.n̂[3] * face.ΔS / cells[face.owner].Ω
        ∂Δp∂x0[face.neighbour, 1] -= pₙ * face.n̂[1] * face.ΔS / cells[face.neighbour].Ω
        ∂Δp∂x0[face.neighbour, 2] -= pₙ * face.n̂[2] * face.ΔS / cells[face.neighbour].Ω
        ∂Δp∂x0[face.neighbour, 3] -= pₙ * face.n̂[3] * face.ΔS / cells[face.neighbour].Ω
    end

    for face in faces_boundary
        pₙ = cells[face.owner].var[👉.p]
        ∂Δp∂x0[face.owner, 1] += pₙ * face.n̂[1] * face.ΔS / cells[face.owner].Ω
        ∂Δp∂x0[face.owner, 2] += pₙ * face.n̂[2] * face.ΔS / cells[face.owner].Ω
        ∂Δp∂x0[face.owner, 3] += pₙ * face.n̂[3] * face.ΔS / cells[face.owner].Ω
    end




    for face in faces_internal

        ΔS = face.ΔS

        ijStartₗ = B_n*(face.owner-1)
        ijStartᵣ = B_n*(face.neighbour-1)

        ρₗ = cells[face.owner].var[👉.ρ]
        ρᵣ = cells[face.neighbour].var[👉.ρ]
        pₗ = cells[face.owner].var[👉.p]
        pᵣ = cells[face.neighbour].var[👉.p]
        uₗ = cells[face.owner].var[👉.u]
        uᵣ = cells[face.neighbour].var[👉.u]
        vₗ = cells[face.owner].var[👉.v]
        vᵣ = cells[face.neighbour].var[👉.v]
        wₗ = cells[face.owner].var[👉.w]
        wᵣ = cells[face.neighbour].var[👉.w]
        Hₜₗ = cells[face.owner].var[👉.Hₜ]
        Hₜᵣ = cells[face.neighbour].var[👉.Hₜ]
        #μₗ = cells[face.owner].var[👉.μ]
        #μᵣ = cells[face.neighbour].var[👉.μ]
        ∂ρ∂pₗ = cells[face.owner].var[👉.∂ρ∂p]
        ∂ρ∂pᵣ = cells[face.neighbour].var[👉.∂ρ∂p]
        ∂ρ∂Tₗ = cells[face.owner].var[👉.∂ρ∂T]
        ∂ρ∂Tᵣ = cells[face.neighbour].var[👉.∂ρ∂T]
        ∂Hₜ∂pₗ = cells[face.owner].var[👉.∂Hₜ∂p]
        ∂Hₜ∂pᵣ = cells[face.neighbour].var[👉.∂Hₜ∂p]
        ∂Hₜ∂Tₗ = cells[face.owner].var[👉.∂Hₜ∂T]
        ∂Hₜ∂Tᵣ = cells[face.neighbour].var[👉.∂Hₜ∂T]
        Y₁ₗ = cells[face.owner].var[👉.Y₁]
        Y₁ᵣ = cells[face.neighbour].var[👉.Y₁]
        ∂ρ∂Y₁ₗ = cells[face.owner].var[👉.∂ρ∂Y₁]
        ∂ρ∂Y₁ᵣ = cells[face.neighbour].var[👉.∂ρ∂Y₁]
        ∂Hₜ∂Y₁ₗ = cells[face.owner].var[👉.∂Hₜ∂Y₁]
        ∂Hₜ∂Y₁ᵣ = cells[face.neighbour].var[👉.∂Hₜ∂Y₁]
        cₗ = cells[face.owner].var[👉.c]
        cᵣ = cells[face.neighbour].var[👉.c]

        Uₙₗ = uₗ * face.n̂[1] + vₗ * face.n̂[2]
        Uₙᵣ = uᵣ * face.n̂[1] + vᵣ * face.n̂[2]
        Uₙ = 0.5 * (Uₙₗ + Uₙᵣ)
        c̄ = 0.5*(cₗ+cᵣ)

        ρˢ = 1.0 / (0.5/ρₗ + 0.5/ρᵣ)
        
        # roe
        Kp = 0.5
        #Uₙ -= Kp * (pᵣ-pₗ) / c̄ / ρˢ

        

        centerₗ = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centerᵣ = [cells[face.neighbour].x, cells[face.neighbour].y, cells[face.neighbour].z]
        ΔLR = norm(centerᵣ - centerₗ)

        #ρˢ = 1.0 / (0.5/ρₗ + 0.5/ρᵣ)
        #d = 0.5 * (cells[face.owner].Ω / (Ap[face.owner]+1.e-250) + cells[face.neighbour].Ω / (Ap[face.neighbour]+1.e-250) )
        d̂ = 👉.Δt / ρˢ
        #d̂ = d / (1.0 + ρˢ / 👉.Δt * d)
        
        # Rhie-Chow
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 1] * face.n̂[1]
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 2] * face.n̂[2]
        Uₙ += d̂ * ρˢ * 0.5 / ρₗ * ∂Δp∂x0[face.owner, 3] * face.n̂[3]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 1] * face.n̂[1]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 2] * face.n̂[2]
        Uₙ += d̂ * ρˢ * 0.5 / ρᵣ * ∂Δp∂x0[face.neighbour, 3] * face.n̂[3]
        Uₙ -= d̂ * (pᵣ-pₗ) / ΔLR


        
        Wₗ = 0.5 * (1.0 + sign(Uₙ))
        Wᵣ = 1.0 - Wₗ

        ρₙ = Wₗ * ρₗ + Wᵣ * ρᵣ
        uₙ = Wₗ * uₗ + Wᵣ * uᵣ
        vₙ = Wₗ * vₗ + Wᵣ * vᵣ
        wₙ = 0.0#Wₗ * wₗ + Wᵣ * wᵣ
        Hₜₙ = Wₗ * Hₜₗ + Wᵣ * Hₜᵣ
        Y₁ₙ = Wₗ * Y₁ₗ + Wᵣ * Y₁ᵣ

        pₙ = 0.5 * (pₗ + pᵣ)

        
        iₗ = A_n*(face.owner-1)
        iᵣ = A_n*(face.neighbour-1)

        ∂U∂p_diff = d̂ / ΔLR #Kp / c̄ / ρˢ

        #------------------------
        # continuity
        # p'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * Uₙ * ΔS + ρₙ * ∂U∂p_diff * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ * Uₙ * ΔS - ρₙ * ∂U∂p_diff * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * Uₙ * ΔS - ρₙ * ∂U∂p_diff * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ * Uₙ * ΔS + ρₙ * ∂U∂p_diff * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1
        
        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[1] * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[1] * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[1] * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[1] * ΔS ))
        
        # v'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[2] * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[2] * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[2] * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[2] * ΔS ))
        
        # T'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -( Wₗ * ∂ρ∂Tₗ * Uₙ * ΔS ))
        
        # Y₁'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂Y₁ₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 1); push!(A_cols, ijStartᵣ + 5)
        push!(A_vals, ( Wᵣ * ∂ρ∂Y₁ᵣ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Y₁ᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 1); push!(A_cols, ijStartₗ + 5)
        push!(A_vals, -( Wₗ * ∂ρ∂Y₁ₗ * Uₙ * ΔS ))
        


        

        #------------------------
        # x-momentum

        # p'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS + ρₙ * ∂U∂p_diff * uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS - ρₙ * ∂U∂p_diff * uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS - ρₙ * ∂U∂p_diff * uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ * uₙ * Uₙ * ΔS + 0.5 * face.n̂[1] * ΔS + ρₙ * ∂U∂p_diff * uₙ * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[1] * uₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[1] * uₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[1] * uₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[1] * uₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS ))
        
        # v'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[2] * uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[2] * uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[2] * uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[2] * uₙ * ΔS ))

        # T'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * uₙ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ * uₙ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * uₙ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -( Wₗ * ∂ρ∂Tₗ * uₙ * Uₙ * ΔS ))

        # Y₁'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂Y₁ₗ * uₙ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 2); push!(A_cols, ijStartᵣ + 5)
        push!(A_vals, ( Wᵣ * ∂ρ∂Y₁ᵣ * uₙ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Y₁ᵣ * uₙ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 2); push!(A_cols, ijStartₗ + 5)
        push!(A_vals, -( Wₗ * ∂ρ∂Y₁ₗ * uₙ * Uₙ * ΔS ))



        #------------------------
        # y-momentum
        
        # p'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] +=  ( Wₗ * ∂ρ∂pₗ* vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS + ρₙ * ∂U∂p_diff * vₙ * ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ * vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS - ρₙ * ∂U∂p_diff * vₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ* vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS - ρₙ * ∂U∂p_diff * vₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ * vₙ * Uₙ * ΔS + 0.5 * face.n̂[2] * ΔS + ρₙ * ∂U∂p_diff * vₙ * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[1] * vₙ * ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[1] * vₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[1] * vₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[1] * vₙ * ΔS ))
        
        # v'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[2] * vₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[2] * vₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[2] * vₙ * ΔS + Wᵣ * ρᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[2] * vₙ * ΔS + Wₗ * ρₗ * Uₙ * ΔS ))

        # T'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * vₙ * Uₙ *ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ * vₙ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * vₙ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -(  Wₗ * ∂ρ∂Tₗ * vₙ * Uₙ *ΔS ))

        # Y₁'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂Y₁ₗ * vₙ * Uₙ *ΔS )
        push!(A_rows, ijStartₗ + 3); push!(A_cols, ijStartᵣ + 5)
        push!(A_vals, ( Wᵣ * ∂ρ∂Y₁ᵣ * vₙ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Y₁ᵣ * vₙ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 3); push!(A_cols, ijStartₗ + 5)
        push!(A_vals, -(  Wₗ * ∂ρ∂Y₁ₗ * vₙ * Uₙ *ΔS ))


        

        #------------------------
        # energy
        # p'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * Hₜₙ * Uₙ * ΔS + ρₙ * Wₗ * ∂Hₜ∂pₗ * Uₙ * ΔS + ρₙ * ∂U∂p_diff * Hₜₙ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ * Hₜₙ * Uₙ * ΔS + ρₙ * Wᵣ * ∂Hₜ∂pᵣ * Uₙ * ΔS - ρₙ * ∂U∂p_diff * Hₜₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * Hₜₙ * Uₙ * ΔS + ρₙ * Wᵣ * ∂Hₜ∂pᵣ * Uₙ * ΔS - ρₙ * ∂U∂p_diff * Hₜₙ * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ * Hₜₙ * Uₙ * ΔS + ρₙ * Wₗ * ∂Hₜ∂pₗ * Uₙ * ΔS + ρₙ * ∂U∂p_diff * Hₜₙ * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * Wₗ * uₗ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * Wᵣ * uᵣ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * Wᵣ * uᵣ * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * Wₗ * uₗ * ΔS ))

        # v'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * Wₗ * vₗ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * Wᵣ * vᵣ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * Wᵣ * vᵣ * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * Wₗ * vₗ * ΔS ))

        
        # T'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * Hₜₗ * Uₙ * ΔS + ρₙ * Wₗ * ∂Hₜ∂Tₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ * Hₜᵣ * Uₙ * ΔS + ρₙ * Wᵣ * ∂Hₜ∂Tᵣ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * Hₜᵣ * Uₙ * ΔS + ρₙ * Wᵣ * ∂Hₜ∂Tᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -( Wₗ * ∂ρ∂Tₗ * Hₜₗ * Uₙ * ΔS + ρₙ * Wₗ * ∂Hₜ∂Tₗ * Uₙ * ΔS ))

        
        # Y₁'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂Y₁ₗ * Hₜₗ * Uₙ * ΔS + ρₙ * Wₗ * ∂Hₜ∂Y₁ₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 4); push!(A_cols, ijStartᵣ + 5)
        push!(A_vals, ( Wᵣ * ∂ρ∂Y₁ᵣ * Hₜᵣ * Uₙ * ΔS + ρₙ * Wᵣ * ∂Hₜ∂Y₁ᵣ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Y₁ᵣ * Hₜᵣ * Uₙ * ΔS + ρₙ * Wᵣ * ∂Hₜ∂Y₁ᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 4); push!(A_cols, ijStartₗ + 5)
        push!(A_vals, -( Wₗ * ∂ρ∂Y₁ₗ * Hₜₗ * Uₙ * ΔS + ρₙ * Wₗ * ∂Hₜ∂Y₁ₗ * Uₙ * ΔS ))
        

        

        #------------------------
        # mass fraction
        # p'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂pₗ * Y₁ₙ * Uₙ * ΔS + ρₙ * ∂U∂p_diff * Y₁ₙ * ΔS )
        push!(A_rows, ijStartₗ + 5); push!(A_cols, ijStartᵣ + 1)
        push!(A_vals, ( Wᵣ * ∂ρ∂pᵣ * Y₁ₙ * Uₙ * ΔS - ρₙ * ∂U∂p_diff * Y₁ₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂pᵣ * Y₁ₙ * Uₙ * ΔS - ρₙ * ∂U∂p_diff * Y₁ₙ * ΔS )
        push!(A_rows, ijStartᵣ + 5); push!(A_cols, ijStartₗ + 1)
        push!(A_vals, -( Wₗ * ∂ρ∂pₗ * Y₁ₙ * Uₙ * ΔS + ρₙ * ∂U∂p_diff * Y₁ₙ * ΔS ))
        
        # u'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[1] * Y₁ₙ * ΔS )
        push!(A_rows, ijStartₗ + 5); push!(A_cols, ijStartᵣ + 2)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[1] * Y₁ₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[1] * Y₁ₙ * ΔS )
        push!(A_rows, ijStartᵣ + 5); push!(A_cols, ijStartₗ + 2)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[1] * Y₁ₙ * ΔS ))

        # v'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( ρₙ * 0.5 * face.n̂[2] * Y₁ₙ * ΔS )
        push!(A_rows, ijStartₗ + 5); push!(A_cols, ijStartᵣ + 3)
        push!(A_vals, ( ρₙ * 0.5 * face.n̂[2] * Y₁ₙ * ΔS ))
        
        A_vals[iᵣ] -= ( ρₙ * 0.5 * face.n̂[2] * Y₁ₙ * ΔS )
        push!(A_rows, ijStartᵣ + 5); push!(A_cols, ijStartₗ + 3)
        push!(A_vals, -( ρₙ * 0.5 * face.n̂[2] * Y₁ₙ * ΔS ))

        
        # T'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂Tₗ * Y₁ₙ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 5); push!(A_cols, ijStartᵣ + 4)
        push!(A_vals, ( Wᵣ * ∂ρ∂Tᵣ * Y₁ₙ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Tᵣ * Y₁ₙ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 5); push!(A_cols, ijStartₗ + 4)
        push!(A_vals, -( Wₗ * ∂ρ∂Tₗ * Y₁ₙ * Uₙ * ΔS ))

        
        # Y₁'
        iₗ += 1; iᵣ += 1

        A_vals[iₗ] += ( Wₗ * ∂ρ∂Y₁ₗ * Y₁ₙ * Uₙ * ΔS + ρₙ * Wₗ * Uₙ * ΔS )
        push!(A_rows, ijStartₗ + 5); push!(A_cols, ijStartᵣ + 5)
        push!(A_vals, ( Wᵣ * ∂ρ∂Y₁ᵣ * Y₁ₙ * Uₙ * ΔS + ρₙ * Wᵣ * Uₙ * ΔS ))
        
        A_vals[iᵣ] -= ( Wᵣ * ∂ρ∂Y₁ᵣ * Y₁ₙ * Uₙ * ΔS + ρₙ * Wᵣ * Uₙ * ΔS )
        push!(A_rows, ijStartᵣ + 5); push!(A_cols, ijStartₗ + 5)
        push!(A_vals, -( Wₗ * ∂ρ∂Y₁ₗ * Y₁ₙ * Uₙ * ΔS + ρₙ * Wₗ * Uₙ * ΔS ))
        

        # ----------------------------
#=
        # B
        B[ijStartₗ + 1] -= ( ρₙ * Uₙ * ΔS )
        B[ijStartᵣ + 1] += ( ρₙ * Uₙ * ΔS )
        # B
        B[ijStartₗ + 2] -= ( ρₙ * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )
        B[ijStartᵣ + 2] += ( ρₙ * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )
        # B
        B[ijStartₗ + 3] -= ( ρₙ * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )
        B[ijStartᵣ + 3] += ( ρₙ * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )
        # B
        B[ijStartₗ + 4] -= ( ρₙ * Hₜₙ * Uₙ * ΔS )
        B[ijStartᵣ + 4] += ( ρₙ * Hₜₙ * Uₙ * ΔS )
        # B
        B[ijStartₗ + 5] -= ( ρₙ * Y₁ₙ * Uₙ * ΔS )
        B[ijStartᵣ + 5] += ( ρₙ * Y₁ₙ * Uₙ * ΔS )
=#


    end



    # boundary faces
    #boundary = append(faces_boundary_top , faces_boundary_bottom , faces_boundary_left , faces_boundary_right )
    bc_wall = []
    append!( bc_wall, faces_boundary_top )
    append!( bc_wall, faces_boundary_bottom )
    append!( bc_wall, faces_boundary_left )
    append!( bc_wall, faces_boundary_right )

    bc_slipwall = []
    
    bc_subinlet = []
    
    bc_suboutlet = []
    
    bc_supoutlet = []

    for face in bc_wall
        
        ijStartₗ = B_n*(face.owner-1)

        i = A_n*(face.owner-1)

        ρₙ = cells[face.owner].var[👉.ρ]
        ∂ρ∂pₙ = cells[face.owner].var[👉.∂ρ∂p]
        ∂ρ∂Tₙ = cells[face.owner].var[👉.∂ρ∂T]
        ∂Hₜ∂pₙ = cells[face.owner].var[👉.∂Hₜ∂p]
        ∂Hₜ∂Tₙ = cells[face.owner].var[👉.∂Hₜ∂T]
        ∂Hₜ∂Y₁ₙ = cells[face.owner].var[👉.∂Hₜ∂Y₁]
        ∂ρ∂Y₁ₙ = cells[face.owner].var[👉.∂ρ∂Y₁]
        pₙ = cells[face.owner].var[👉.p]
        Hₜₙ = cells[face.owner].var[👉.Hₜ]
        Y₁ₙ = cells[face.owner].var[👉.Y₁]

        ΔS = face.ΔS

        uₙ = 0.0
        vₙ = 0.0
        wₙ = 0.0
        Uₙ = 0.0
        Tₙ = cells[face.owner].var[👉.T]
        α₁ₙ = cells[face.owner].var[👉.α₁]


        ρₙ, Hₜₙ, cₙ = faceEOS!(pₙ,uₙ,vₙ,wₙ,Tₙ,α₁ₙ)

        # continuity
        i += 1
        A_vals[i] += ∂ρ∂pₙ * Uₙ * ΔS
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * Uₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Y₁ₙ * Uₙ * ΔS

        
        # x-momentum
        i += 1
        A_vals[i] += ∂ρ∂pₙ * uₙ * Uₙ * ΔS
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * uₙ * Uₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Y₁ₙ * uₙ * Uₙ * ΔS

        
        # y-momentum
        i += 1
        A_vals[i] += ∂ρ∂pₙ * vₙ * Uₙ * ΔS
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * vₙ * Uₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Y₁ₙ * vₙ * Uₙ * ΔS


        # energy
        i += 1
        A_vals[i] += ∂ρ∂pₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂pₙ * ΔS
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂Tₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Y₁ₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂Y₁ₙ * ΔS


        # massfraction
        i += 1
        A_vals[i] += (∂ρ∂pₙ * Uₙ * Y₁ₙ * ΔS)# + ρₙ * Hₜₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * Uₙ * Y₁ₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Y₁ₙ * Uₙ * Y₁ₙ * ΔS + ρₙ * Uₙ * ΔS

#=
        B[ijStartₗ + 1] -= ( ρₙ * Uₙ * ΔS )
        B[ijStartₗ + 2] -= ( ρₙ * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )
        B[ijStartₗ + 3] -= ( ρₙ * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )
        B[ijStartₗ + 4] -= ( ρₙ * Hₜₙ * Uₙ * ΔS )
        B[ijStartₗ + 5] -= ( ρₙ * Y₁ₙ * Uₙ * ΔS )
=#
        

    end
 

    for face in bc_slipwall
        
        ijStartₗ = B_n*(face.owner-1)

        i = A_n*(face.owner-1)

        ρₙ = cells[face.owner].var[👉.ρ]
        ∂ρ∂pₙ = cells[face.owner].var[👉.∂ρ∂p]
        ∂ρ∂Tₙ = cells[face.owner].var[👉.∂ρ∂T]
        ∂Hₜ∂pₙ = cells[face.owner].var[👉.∂Hₜ∂p]
        ∂Hₜ∂Tₙ = cells[face.owner].var[👉.∂Hₜ∂T]
        ∂Hₜ∂Y₁ₙ = cells[face.owner].var[👉.∂Hₜ∂Y₁]
        ∂ρ∂Y₁ₙ = cells[face.owner].var[👉.∂ρ∂Y₁]
        pₙ = cells[face.owner].var[👉.p]
        Hₜₙ = cells[face.owner].var[👉.Hₜ]
        Y₁ₙ = cells[face.owner].var[👉.Y₁]

        ΔS = face.ΔS

        Uₙ = 0.0
        Uₙ += cells[face.owner].var[👉.u]*face.n̂[1]
        Uₙ += cells[face.owner].var[👉.v]*face.n̂[2]
        Uₙ += cells[face.owner].var[👉.w]*face.n̂[3]

        uₙ = cells[face.owner].var[👉.u] - Uₙ * face.n̂[1]
        vₙ = cells[face.owner].var[👉.v] - Uₙ * face.n̂[2]
        wₙ = cells[face.owner].var[👉.w] - Uₙ * face.n̂[3]

        Uₙ = uₙ * face.n̂[1] + vₙ * face.n̂[2] + wₙ * face.n̂[3]

        Tₙ = cells[face.owner].var[👉.T]
        α₁ₙ = cells[face.owner].var[👉.α₁]


        ρₙ, Hₜₙ, cₙ = faceEOS!(pₙ,uₙ,vₙ,wₙ,Tₙ,α₁ₙ)

        # continuity
        i += 1
        A_vals[i] += ∂ρ∂pₙ * Uₙ * ΔS
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * Uₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Y₁ₙ * Uₙ * ΔS

        
        # x-momentum
        i += 1
        A_vals[i] += ∂ρ∂pₙ * uₙ * Uₙ * ΔS
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * uₙ * Uₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Y₁ₙ * uₙ * Uₙ * ΔS

        
        # y-momentum
        i += 1
        A_vals[i] += ∂ρ∂pₙ * vₙ * Uₙ * ΔS
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * vₙ * Uₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Y₁ₙ * vₙ * Uₙ * ΔS


        # energy
        i += 1
        A_vals[i] += ∂ρ∂pₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂pₙ * ΔS
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂Tₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Y₁ₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂Y₁ₙ * ΔS


        # massfraction
        i += 1
        A_vals[i] += (∂ρ∂pₙ * Uₙ * Y₁ₙ * ΔS)# + ρₙ * Hₜₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * Uₙ * Y₁ₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Y₁ₙ * Uₙ * Y₁ₙ * ΔS + ρₙ * Uₙ * ΔS



        B[ijStartₗ + 1] -= ( ρₙ * Uₙ * ΔS )
        B[ijStartₗ + 2] -= ( ρₙ * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )
        B[ijStartₗ + 3] -= ( ρₙ * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )
        B[ijStartₗ + 4] -= ( ρₙ * Hₜₙ * Uₙ * ΔS )
        B[ijStartₗ + 5] -= ( ρₙ * Y₁ₙ * Uₙ * ΔS )
        

    end
 
    for face in bc_subinlet
        
        ijStartₗ = B_n*(face.owner-1)

        i = A_n*(face.owner-1)

        ρₙ = cells[face.owner].var[👉.ρ]
        ∂ρ∂pₙ = cells[face.owner].var[👉.∂ρ∂p]
        ∂ρ∂Tₙ = cells[face.owner].var[👉.∂ρ∂T]
        ∂Hₜ∂pₙ = cells[face.owner].var[👉.∂Hₜ∂p]
        ∂Hₜ∂Tₙ = cells[face.owner].var[👉.∂Hₜ∂T]
        ∂Hₜ∂Y₁ₙ = cells[face.owner].var[👉.∂Hₜ∂Y₁]
        ∂ρ∂Y₁ₙ = cells[face.owner].var[👉.∂ρ∂Y₁]
        pₙ = cells[face.owner].var[👉.p]
        Hₜₙ = cells[face.owner].var[👉.Hₜ]
        Y₁ₙ = cells[face.owner].var[👉.Y₁]

        ΔS = face.ΔS

        uₙ = 1.0
        #uₙ = 0.5 * ( 1.0 + cells[face.owner].var[👉.u] )
        vₙ = 0.0
        wₙ = 0.0
        Uₙ = uₙ*face.n̂[1] + vₙ*face.n̂[2]

        Tₙ = 300.0
        #Tₙ = 0.5 * ( 300.0 + cells[face.owner].var[👉.T] )
        Y₁ₙ = 1.0
        
        #pₙ = 101325.0

        ρₙ, Hₜₙ, cₙ = faceEOS!(pₙ,uₙ,vₙ,wₙ,Tₙ,Y₁ₙ)

        centerₗ = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centerᵣ = [face.x, face.y, face.z]
        ΔLR = 1.0 * norm(centerᵣ - centerₗ)

        # continuity
        i += 1
        A_vals[i] += ∂ρ∂pₙ * Uₙ * ΔS# + ρₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += 0.0#0.5 * (ρₙ * face.n̂[1] * ΔS)
        i += 1
        A_vals[i] += 0.0#0.5 * (ρₙ * face.n̂[2] * ΔS)
        i += 1
        A_vals[i] += 0.0#0.5 * (∂ρ∂Tₙ * Uₙ * ΔS)
        i += 1
        A_vals[i] += 0.0#0.5 * (∂ρ∂Tₙ * Uₙ * ΔS)

        
        # x-momentum
        i += 1
        A_vals[i] += ∂ρ∂pₙ * uₙ * Uₙ * ΔS# + ρₙ * uₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += 0.0#0.5 * (ρₙ * Uₙ * ΔS + ρₙ * uₙ * face.n̂[1] * ΔS)
        i += 1
        A_vals[i] += 0.0#0.5 * (ρₙ * uₙ * face.n̂[2] * ΔS)
        i += 1
        A_vals[i] += 0.0#0.5 * (∂ρ∂Tₙ * uₙ * Uₙ * ΔS)
        i += 1
        A_vals[i] += 0.0#0.5 * (∂ρ∂Tₙ * uₙ * Uₙ * ΔS)

        
        # y-momentum
        i += 1
        A_vals[i] += ∂ρ∂pₙ * vₙ * Uₙ * ΔS# + ρₙ * vₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += 0.0#0.5 * ρₙ * vₙ * face.n̂[1] * ΔS
        i += 1
        A_vals[i] += 0.0#0.5 * (ρₙ * Uₙ * ΔS + ρₙ * vₙ * face.n̂[2] * ΔS)
        i += 1
        A_vals[i] += 0.0#0.5 * (∂ρ∂Tₙ * vₙ * Uₙ * ΔS)
        i += 1
        A_vals[i] += 0.0#0.5 * (∂ρ∂Tₙ * vₙ * Uₙ * ΔS)


        # energy
        i += 1
        A_vals[i] += (∂ρ∂pₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂pₙ * ΔS)# + ρₙ * Hₜₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += 0.0#0.5 * (ρₙ * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * uₙ * ΔS)
        i += 1
        A_vals[i] += 0.0#0.5 * (ρₙ * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * vₙ * ΔS)
        i += 1
        A_vals[i] += 0.0#0.5 * (∂ρ∂Tₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂Tₙ * ΔS)
        i += 1
        A_vals[i] += 0.0#0.5 * (∂ρ∂Tₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂Tₙ * ΔS)


        # massfraction
        i += 1
        A_vals[i] += (∂ρ∂pₙ * Uₙ * Y₁ₙ * ΔS)# + ρₙ * Hₜₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += 0.0#
        i += 1
        A_vals[i] += 0.0#ρₙ * face.n̂[2] * Y₁ₙ * ΔS
        i += 1
        A_vals[i] += 0.0#∂ρ∂Tₙ * Uₙ * Y₁ₙ * ΔS
        i += 1
        A_vals[i] += 0.0#∂ρ∂Y₁ₙ * Uₙ * Y₁ₙ * ΔS + ρₙ * Uₙ * ΔS


        B[ijStartₗ + 1] -= ( ρₙ * Uₙ * ΔS )
        B[ijStartₗ + 2] -= ( ρₙ * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )
        B[ijStartₗ + 3] -= ( ρₙ * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )
        B[ijStartₗ + 4] -= ( ρₙ * Hₜₙ * Uₙ * ΔS )
        B[ijStartₗ + 5] -= ( ρₙ * Y₁ₙ * Uₙ * ΔS )
        

    end
 
    for face in bc_suboutlet
        
        ijStartₗ = B_n*(face.owner-1)

        i = A_n*(face.owner-1)

        ρₙ = cells[face.owner].var[👉.ρ]
        ∂ρ∂pₙ = cells[face.owner].var[👉.∂ρ∂p]
        ∂ρ∂Tₙ = cells[face.owner].var[👉.∂ρ∂T]
        ∂Hₜ∂pₙ = cells[face.owner].var[👉.∂Hₜ∂p]
        ∂Hₜ∂Tₙ = cells[face.owner].var[👉.∂Hₜ∂T]
        ∂Hₜ∂Y₁ₙ = cells[face.owner].var[👉.∂Hₜ∂Y₁]
        ∂ρ∂Y₁ₙ = cells[face.owner].var[👉.∂ρ∂Y₁]
        pₙ = cells[face.owner].var[👉.p]
        Hₜₙ = cells[face.owner].var[👉.Hₜ]
        Y₁ₙ = cells[face.owner].var[👉.Y₁]

        ΔS = face.ΔS

        uₙ = cells[face.owner].var[👉.u]
        vₙ = cells[face.owner].var[👉.v]
        wₙ = cells[face.owner].var[👉.w]
        Uₙ = uₙ*face.n̂[1] + vₙ*face.n̂[2]

        Tₙ = cells[face.owner].var[👉.T]
        α₁ₙ = cells[face.owner].var[👉.α₁]

        pₙ = 101325.0 #0.5 * ( 101325.0 + cells[face.owner].var[👉.p] )
        
        ρₙ, Hₜₙ, cₙ = faceEOS!(pₙ,uₙ,vₙ,wₙ,Tₙ,α₁ₙ)

        centerₗ = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centerᵣ = [face.x, face.y, face.z]
        ΔLR = 2.0 * norm(centerᵣ - centerₗ)

        # continuity
        i += 1
        A_vals[i] += 0.0#0.5 * (∂ρ∂pₙ * Uₙ * ΔS) + ρₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += ρₙ * face.n̂[1] * ΔS
        i += 1
        A_vals[i] += ρₙ * face.n̂[2] * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * Uₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Y₁ₙ * Uₙ * ΔS

        
        # x-momentum
        i += 1
        A_vals[i] += 0.0#0.5 * (∂ρ∂pₙ * uₙ * Uₙ * ΔS) + ρₙ * uₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += (ρₙ * Uₙ * ΔS + ρₙ * uₙ * face.n̂[1] * ΔS)
        i += 1
        A_vals[i] += ρₙ * uₙ * face.n̂[2] * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * uₙ * Uₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Y₁ₙ * uₙ * Uₙ * ΔS

        
        # y-momentum
        i += 1
        A_vals[i] += 0.0#0.5 * (∂ρ∂pₙ * vₙ * Uₙ * ΔS) + ρₙ * vₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += ρₙ * vₙ * face.n̂[1] * ΔS
        i += 1
        A_vals[i] += (ρₙ * Uₙ * ΔS + ρₙ * vₙ * face.n̂[2] * ΔS)
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * vₙ * Uₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Y₁ₙ * vₙ * Uₙ * ΔS


        # energy
        i += 1
        A_vals[i] += 0.0#0.5 * (∂ρ∂pₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂pₙ * ΔS) + ρₙ * Hₜₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += ρₙ * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * uₙ * ΔS
        i += 1
        A_vals[i] += ρₙ * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * vₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂Tₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Y₁ₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂Y₁ₙ * ΔS


        # massfraction
        i += 1
        A_vals[i] += 0.0
        i += 1
        A_vals[i] += ρₙ * face.n̂[1] * Y₁ₙ * ΔS
        i += 1
        A_vals[i] += ρₙ * face.n̂[2] * Y₁ₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * Uₙ * Y₁ₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Y₁ₙ * Uₙ * Y₁ₙ * ΔS + ρₙ * Uₙ * ΔS

        B[ijStartₗ + 1] -= ( ρₙ * Uₙ * ΔS )
        B[ijStartₗ + 2] -= ( ρₙ * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )
        B[ijStartₗ + 3] -= ( ρₙ * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )
        B[ijStartₗ + 4] -= ( ρₙ * Hₜₙ * Uₙ * ΔS )
        B[ijStartₗ + 5] -= ( ρₙ * Y₁ₙ * Uₙ * ΔS )
        

    end
 
    
    
    for face in bc_supoutlet
        
        ijStartₗ = B_n*(face.owner-1)

        i = A_n*(face.owner-1)

        ρₙ = cells[face.owner].var[👉.ρ]
        ∂ρ∂pₙ = cells[face.owner].var[👉.∂ρ∂p]
        ∂ρ∂Tₙ = cells[face.owner].var[👉.∂ρ∂T]
        ∂Hₜ∂pₙ = cells[face.owner].var[👉.∂Hₜ∂p]
        ∂Hₜ∂Tₙ = cells[face.owner].var[👉.∂Hₜ∂T]
        ∂Hₜ∂Y₁ₙ = cells[face.owner].var[👉.∂Hₜ∂Y₁]
        ∂ρ∂Y₁ₙ = cells[face.owner].var[👉.∂ρ∂Y₁]
        pₙ = cells[face.owner].var[👉.p]
        Hₜₙ = cells[face.owner].var[👉.Hₜ]
        Y₁ₙ = cells[face.owner].var[👉.Y₁]

        ΔS = face.ΔS

        uₙ = cells[face.owner].var[👉.u]
        vₙ = cells[face.owner].var[👉.v]
        wₙ = cells[face.owner].var[👉.w]
        Uₙ = uₙ*face.n̂[1] + vₙ*face.n̂[2]

        Tₙ = cells[face.owner].var[👉.T]
        α₁ₙ = cells[face.owner].var[👉.α₁]

        ρₙ, Hₜₙ, cₙ = faceEOS!(pₙ,uₙ,vₙ,wₙ,Tₙ,α₁ₙ)

        centerₗ = [cells[face.owner].x, cells[face.owner].y, cells[face.owner].z]
        centerᵣ = [face.x, face.y, face.z]
        ΔLR = 2.0 * norm(centerᵣ - centerₗ)

        # continuity
        i += 1
        A_vals[i] += (∂ρ∂pₙ * Uₙ * ΔS) #+ ρₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += ρₙ * face.n̂[1] * ΔS
        i += 1
        A_vals[i] += ρₙ * face.n̂[2] * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * Uₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Y₁ₙ * Uₙ * ΔS

        
        # x-momentum
        i += 1
        A_vals[i] += (∂ρ∂pₙ * uₙ * Uₙ * ΔS)# + ρₙ * uₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += (ρₙ * Uₙ * ΔS + ρₙ * uₙ * face.n̂[1] * ΔS)
        i += 1
        A_vals[i] += ρₙ * uₙ * face.n̂[2] * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * uₙ * Uₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Y₁ₙ * uₙ * Uₙ * ΔS

        
        # y-momentum
        i += 1
        A_vals[i] += (∂ρ∂pₙ * vₙ * Uₙ * ΔS)# + ρₙ * vₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += ρₙ * vₙ * face.n̂[1] * ΔS
        i += 1
        A_vals[i] += (ρₙ * Uₙ * ΔS + ρₙ * vₙ * face.n̂[2] * ΔS)
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * vₙ * Uₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Y₁ₙ * vₙ * Uₙ * ΔS


        # energy
        i += 1
        A_vals[i] += (∂ρ∂pₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂pₙ * ΔS)# + ρₙ * Hₜₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += ρₙ * face.n̂[1] * Hₜₙ * ΔS + ρₙ * Uₙ * uₙ * ΔS
        i += 1
        A_vals[i] += ρₙ * face.n̂[2] * Hₜₙ * ΔS + ρₙ * Uₙ * vₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂Tₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Y₁ₙ * Uₙ * Hₜₙ * ΔS + ρₙ * Uₙ * ∂Hₜ∂Y₁ₙ * ΔS


        # massfraction
        i += 1
        A_vals[i] += (∂ρ∂pₙ * Uₙ * Y₁ₙ * ΔS)# + ρₙ * Hₜₙ * 👉.Δt/ρₙ / ΔLR * ΔS
        i += 1
        A_vals[i] += ρₙ * face.n̂[1] * Y₁ₙ * ΔS
        i += 1
        A_vals[i] += ρₙ * face.n̂[2] * Y₁ₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Tₙ * Uₙ * Y₁ₙ * ΔS
        i += 1
        A_vals[i] += ∂ρ∂Y₁ₙ * Uₙ * Y₁ₙ * ΔS + ρₙ * Uₙ * ΔS


        B[ijStartₗ + 1] -= ( ρₙ * Uₙ * ΔS )
        B[ijStartₗ + 2] -= ( ρₙ * uₙ * Uₙ * ΔS + pₙ * face.n̂[1] * ΔS )
        B[ijStartₗ + 3] -= ( ρₙ * vₙ * Uₙ * ΔS + pₙ * face.n̂[2] * ΔS )
        B[ijStartₗ + 4] -= ( ρₙ * Hₜₙ * Uₙ * ΔS )
        B[ijStartₗ + 5] -= ( ρₙ * Y₁ₙ * Uₙ * ΔS )
        

    end
 


    A = sparse(A_rows,A_cols,A_vals)

    ps = MKLPardisoSolver()
    set_matrixtype!(ps, Pardiso.REAL_NONSYM)
    ΔQ = solve(ps, A, B)


    relax_p = 0.9
    relax_U = 0.9
    relax_T = 0.9
    relax_Y = 0.9


    diagon = 1
    for cell in cells

        ijStart = B_n*(diagon-1)
        Astart = A_n*(diagon-1)
        i = Astart

        
        cell.var[👉.p] += relax_p * ΔQ[ijStart + 1]
        cell.var[👉.u] += relax_U * ΔQ[ijStart + 2]
        cell.var[👉.v] += relax_U * ΔQ[ijStart + 3]
        cell.var[👉.T] += relax_T * ΔQ[ijStart + 4]
        cell.var[👉.Y₁] += relax_Y * ΔQ[ijStart + 5]

        cell.var[👉.p] = max(cell.var[👉.p],1.e-200)
        cell.var[👉.T] = max(cell.var[👉.T],1.e-200)
        

        diagon += 1
    end

    return norm(ΔQ)

end