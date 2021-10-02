include("./module.jl")
include("./controls.jl")
using .mesh

function structured_grid_uniform!(
    👉::controls,
    cell::Vector{mesh.Cell},
    face::Vector{mesh.Face},
    face_internal::Vector{mesh.Face},
    face_boundary::Vector{mesh.Face},
    face_boundary_top::Vector{mesh.Face},
    face_boundary_bottom::Vector{mesh.Face},
    face_boundary_left::Vector{mesh.Face},
    face_boundary_right::Vector{mesh.Face}
    )

    N = 👉.Nx*👉.Ny*👉.Nz

    Δx = 👉.Lx/👉.Nx
    Δy = 👉.Ly/👉.Ny 
    Δz = 👉.Lz/👉.Nz 

    x = 0.5*Δx:Δx:👉.Lx
    y = 0.5*Δy:Δy:👉.Ly
    z = 0.5*Δz:Δz:👉.Lz

    length(y)

    for i in 1:N 
        push!(cell, mesh.Cell())
    end

    #println(N)
    #println(length(cell))


    for i in 1:👉.Nx 
        for j in 1:👉.Ny
            k=1
            ijk = i + 👉.Nx*(j-1) + 👉.Nx*👉.Ny*(k-1)
            cell[ijk].x = x[i]
            cell[ijk].y = y[j]
            cell[ijk].z = 0.5*Δz #z[i]

            cell[ijk].Ω = Δx*Δy*Δz

            cell[ijk].Qᵐ = zeros(Float64,6)
            cell[ijk].Qⁿ = zeros(Float64,6)
            cell[ijk].Qⁿ⁻¹ = zeros(Float64,6)
        end
    end



    #println(length(cell))
    # internal faces
    for i in 2:👉.Nx
        for j in 2:👉.Ny
            k = 1
            ijk = i + 👉.Nx*(j-1) + 👉.Nx*👉.Ny*(k-1)
            imjk = (i-1) + 👉.Nx*(j-1) + 👉.Nx*👉.Ny*(k-1)
            ijmk = i + 👉.Nx*(j-2) + 👉.Nx*👉.Ny*(k-1)
            n̂ = [0.0, 1.0, 0.0]
            ΔS = Δx*Δz
            owner = ijmk
            neighbour = ijk
            push!(face, mesh.Face(cell[ijk].x, Δy*(j-1), 0.5*Δz, owner, neighbour, n̂, ΔS, [], [], []))
            n̂ = [1.0, 0.0, 0.0]
            ΔS = Δy*Δz
            owner = imjk
            neighbour = ijk
            push!(face, mesh.Face(Δx*(i-1), cell[ijk].y, 0.5*Δz, owner, neighbour, n̂, ΔS, [], [], []))
        end
    end

    for j in 2:👉.Ny
        i = 1
        k = 1
        ijk = i + 👉.Nx*(j-1) + 👉.Nx*👉.Ny*(k-1)
        imjk = (i-1) + 👉.Nx*(j-1) + 👉.Nx*👉.Ny*(k-1)
        ijmk = i + 👉.Nx*(j-2) + 👉.Nx*👉.Ny*(k-1)
        n̂ = [0.0, 1.0, 0.0]
        ΔS = Δx*Δz
        owner = ijmk
        neighbour = ijk
        push!(face, mesh.Face(cell[ijk].x, Δy*(j-1), 0.5*Δz, owner, neighbour, n̂, ΔS, [], [], []))
    end

    for i in 2:👉.Nx
        j = 1
        k = 1
        ijk = i + 👉.Nx*(j-1) + 👉.Nx*👉.Ny*(k-1)
        imjk = (i-1) + 👉.Nx*(j-1) + 👉.Nx*👉.Ny*(k-1)
        ijmk = i + 👉.Nx*(j-2) + 👉.Nx*👉.Ny*(k-1)
        n̂ = [1.0, 0.0, 0.0]
        ΔS = Δy*Δz
        owner = imjk
        neighbour = ijk
        push!(face, mesh.Face(Δx*(i-1), cell[ijk].y, 0.5*Δz, owner, neighbour, n̂, ΔS, [], [], []))
    end

    face_internal_num = length(face)
    for i in face
        push!(face_internal, i)
    end

    # boundary faces
    # Left
    for j in 1:👉.Ny
        i = 1
        k = 1
        ijk = i + 👉.Nx*(j-1) + 👉.Nx*👉.Ny*(k-1)
        imjk = (i-1) + 👉.Nx*(j-1) + 👉.Nx*👉.Ny*(k-1)
        ijmk = i + 👉.Nx*(j-2) + 👉.Nx*👉.Ny*(k-1)
        n̂ = [-1.0, 0.0, 0.0]
        ΔS = Δy*Δz
        owner = ijk
        neighbour = 0
        push!(face, mesh.Face(Δx*(i-1), cell[ijk].y, 0.5*Δz, owner, neighbour, n̂, ΔS, [], [], []))
    end

    face_total_left = length(face)
    for i in face_internal_num+1:face_total_left
        push!(face_boundary_left, face[i])
    end

    # Bottom
    for i in 1:👉.Nx
        j = 1
        k = 1
        ijk = i + 👉.Nx*(j-1) + 👉.Nx*👉.Ny*(k-1)
        imjk = (i-1) + 👉.Nx*(j-1) + 👉.Nx*👉.Ny*(k-1)
        ijmk = i + 👉.Nx*(j-2) + 👉.Nx*👉.Ny*(k-1)
        n̂ = [0.0, -1.0, 0.0]
        ΔS = Δx*Δz
        owner = ijk
        neighbour = 0
        push!(face, mesh.Face(cell[ijk].x, Δy*(j-1), 0.5*Δz, owner, neighbour, n̂, ΔS, [], [], []))
    end
    
    face_total_bottom = length(face)
    for i in face_total_left+1:face_total_bottom
        push!(face_boundary_bottom, face[i])
    end

    # Right
    for j in 1:👉.Ny
        i = 👉.Nx
        k = 1
        ijk = i + 👉.Nx*(j-1) + 👉.Nx*👉.Ny*(k-1)
        imjk = (i-1) + 👉.Nx*(j-1) + 👉.Nx*👉.Ny*(k-1)
        ijmk = i + 👉.Nx*(j-2) + 👉.Nx*👉.Ny*(k-1)
        n̂ = [1.0, 0.0, 0.0]
        ΔS = Δy*Δz
        owner = ijk
        neighbour = 0
        push!(face, mesh.Face(Δx*👉.Nx, cell[ijk].y, 0.5*Δz, owner, neighbour, n̂, ΔS, [], [], []))
    end
    
    face_total_right = length(face)
    for i in face_total_bottom+1:face_total_right
        push!(face_boundary_right, face[i])
    end

    # Top
    for i in 1:👉.Nx
        j = 👉.Ny
        k = 1
        ijk = i + 👉.Nx*(j-1) + 👉.Nx*👉.Ny*(k-1)
        imjk = (i-1) + 👉.Nx*(j-1) + 👉.Nx*👉.Ny*(k-1)
        ijmk = i + 👉.Nx*(j-2) + 👉.Nx*👉.Ny*(k-1)
        n̂ = [0.0, 1.0, 0.0]
        ΔS = Δx*Δz
        owner = ijk
        neighbour = 0
        push!(face, mesh.Face(cell[ijk].x, Δy*👉.Ny, 0.5*Δz, owner, neighbour, n̂, ΔS, [], [], []))
    end
    
    face_total_top = length(face)
    for i in face_total_right+1:face_total_top
        push!(face_boundary_top, face[i])
    end

    
    face_total = length(face)
    for i in face_internal_num+1:face_total
        push!(face_boundary, face[i])
    end
    

    for i in cell
        for j in 1:31
            push!(i.var,0.0)
        end
    end

    for i in face
        for j in 1:9
            push!(i.varₗ,0.0)
            push!(i.varᵣ,0.0)
        end
    end
    #=
    =#

    return nothing

end