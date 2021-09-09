function face_left_right!(
    👉::controls, 
    cells::Vector{mesh.Cell},
    faces::Vector{mesh.Face},
    faces_internal::Vector{mesh.Face},
    faces_boundary::Vector{mesh.Face},
    faces_boundary_top::Vector{mesh.Face},
    faces_boundary_bottom::Vector{mesh.Face},
    faces_boundary_left::Vector{mesh.Face},
    faces_boundary_right::Vector{mesh.Face}
)

    # internal faces
    for face in faces_internal

        face.varₗ[👉.p] = cells[face.owner].var[👉.p]
        face.varᵣ[👉.p] = cells[face.neighbour].var[👉.p]

        face.varₗ[👉.u] = cells[face.owner].var[👉.u]
        face.varᵣ[👉.u] = cells[face.neighbour].var[👉.u]
        
        face.varₗ[👉.v] = cells[face.owner].var[👉.v]
        face.varᵣ[👉.v] = cells[face.neighbour].var[👉.v]
        
        face.varₗ[👉.w] = cells[face.owner].var[👉.w]
        face.varᵣ[👉.w] = cells[face.neighbour].var[👉.w]
        
        face.varₗ[👉.T] = cells[face.owner].var[👉.T]
        face.varᵣ[👉.T] = cells[face.neighbour].var[👉.T]
        
        face.varₗ[👉.Y₁] = cells[face.owner].var[👉.Y₁]
        face.varᵣ[👉.Y₁] = cells[face.neighbour].var[👉.Y₁]

    end


    # Top
    #bc_sub_outlet!(👉, faces_boundary_top, cells, 101325.0)
    #bc_moving_wall!(👉, faces_boundary_top, cells, 1.0, 0.0, 0.0)
    bc_slip_wall!(👉, faces_boundary_top, cells)

    # Bottom
    #bc_wall!(👉, faces_boundary_bottom, cells)
    bc_slip_wall!(👉, faces_boundary_bottom, cells)

    # Left
    #bc_slip_wall!(👉, faces_boundary_left, cells)
    #bc_wall!(👉, faces_boundary_left, cells)
    bc_sup_outlet!(👉, faces_boundary_left, cells)

    # Right
    #bc_slip_wall!(👉, faces_boundary_right, cells)
    #bc_wall!(👉, faces_boundary_right, cells)
    bc_sup_outlet!(👉, faces_boundary_right, cells)


    # EOS
    for face in faces
        face.varₗ[👉.ρ], face.varₗ[👉.Hₜ], face.varₗ[👉.c] = 
        faceEOS!(face.varₗ[👉.p], 
        face.varₗ[👉.u], face.varₗ[👉.v], face.varₗ[👉.w], face.varₗ[👉.T], face.varₗ[👉.Y₁])
        
        face.varᵣ[👉.ρ], face.varᵣ[👉.Hₜ], face.varᵣ[👉.c] = 
        faceEOS!(face.varᵣ[👉.p], 
        face.varᵣ[👉.u], face.varᵣ[👉.v], face.varᵣ[👉.w], face.varᵣ[👉.T], face.varᵣ[👉.Y₁])

    end


end


function bc_wall!(
    👉::controls, 
    faces_boundary::Vector{mesh.Face},
    cells::Vector{mesh.Cell}
)
    for face in faces_boundary

        face.varₗ[👉.p] = cells[face.owner].var[👉.p]
        face.varᵣ[👉.p] = face.varₗ[👉.p]

        face.varₗ[👉.u] = 0.0
        face.varᵣ[👉.u] = 0.0
        
        face.varₗ[👉.v] = 0.0
        face.varᵣ[👉.v] = 0.0
        
        face.varₗ[👉.w] = 0.0
        face.varᵣ[👉.w] = 0.0
        
        face.varₗ[👉.T] = cells[face.owner].var[👉.T]
        face.varᵣ[👉.T] = face.varₗ[👉.T]
        
        face.varₗ[👉.Y₁] = cells[face.owner].var[👉.Y₁]
        face.varᵣ[👉.Y₁] = face.varₗ[👉.Y₁]

    end

end

function bc_slip_wall!(
    👉::controls, 
    faces_boundary::Vector{mesh.Face},
    cells::Vector{mesh.Cell}
)

    for face in faces_boundary

        face.varₗ[👉.p] = cells[face.owner].var[👉.p]
        face.varᵣ[👉.p] = face.varₗ[👉.p]

        Un = 0.0
        Un += cells[face.owner].var[👉.u]*face.n̂[1]
        Un += cells[face.owner].var[👉.v]*face.n̂[2]
        Un += cells[face.owner].var[👉.w]*face.n̂[3]

        invU = cells[face.owner].var[👉.u] - Un * face.n̂[1]
        invV = cells[face.owner].var[👉.v] - Un * face.n̂[2]
        invW = cells[face.owner].var[👉.w] - Un * face.n̂[3]

        face.varₗ[👉.u] = invU
        face.varᵣ[👉.u] = invU
        
        face.varₗ[👉.v] = invV
        face.varᵣ[👉.v] = invV
        
        face.varₗ[👉.w] = invW
        face.varᵣ[👉.w] = invW
        
        face.varₗ[👉.T] = cells[face.owner].var[👉.T]
        face.varᵣ[👉.T] = face.varₗ[👉.T]
        
        face.varₗ[👉.Y₁] = cells[face.owner].var[👉.Y₁]
        face.varᵣ[👉.Y₁] = face.varₗ[👉.Y₁]

    end

end


    
function bc_moving_wall!(
    👉::controls, 
    faces_boundary::Vector{mesh.Face},
    cells::Vector{mesh.Cell},
    u::Float64, v::Float64, w::Float64
)

    for face in faces_boundary

        face.varₗ[👉.p] = cells[face.owner].var[👉.p]
        face.varᵣ[👉.p] = face.varₗ[👉.p]

        face.varₗ[👉.u] = u
        face.varᵣ[👉.u] = u
        
        face.varₗ[👉.v] = v
        face.varᵣ[👉.v] = v
        
        face.varₗ[👉.w] = w
        face.varᵣ[👉.w] = w
        
        face.varₗ[👉.T] = cells[face.owner].var[👉.T]
        face.varᵣ[👉.T] = face.varₗ[👉.T]
        
        face.varₗ[👉.Y₁] = cells[face.owner].var[👉.Y₁]
        face.varᵣ[👉.Y₁] = face.varₗ[👉.Y₁]

    end

end


function bc_sub_outlet!(
    👉::controls, 
    faces_boundary::Vector{mesh.Face},
    cells::Vector{mesh.Cell},
    set_p::Float64
)

    for face in faces_boundary

        face.varₗ[👉.p] = set_p
        face.varᵣ[👉.p] = set_p

        face.varₗ[👉.u] = cells[face.owner].var[👉.u]
        face.varᵣ[👉.u] = face.varₗ[👉.u]
        
        face.varₗ[👉.v] = cells[face.owner].var[👉.v]
        face.varᵣ[👉.v] = face.varₗ[👉.v]
        
        face.varₗ[👉.w] = cells[face.owner].var[👉.w]
        face.varᵣ[👉.w] = face.varₗ[👉.w]
        
        face.varₗ[👉.T] = cells[face.owner].var[👉.T]
        face.varᵣ[👉.T] = face.varₗ[👉.T]
        
        face.varₗ[👉.Y₁] = cells[face.owner].var[👉.Y₁]
        face.varᵣ[👉.Y₁] = face.varₗ[👉.Y₁]

    end


end
    


function bc_sup_outlet!(
    👉::controls, 
    faces_boundary::Vector{mesh.Face},
    cells::Vector{mesh.Cell}
)

    for face in faces_boundary

        face.varₗ[👉.p] = cells[face.owner].var[👉.p]
        face.varᵣ[👉.p] = face.varₗ[👉.p]

        face.varₗ[👉.u] = cells[face.owner].var[👉.u]
        face.varᵣ[👉.u] = face.varₗ[👉.u]
        
        face.varₗ[👉.v] = cells[face.owner].var[👉.v]
        face.varᵣ[👉.v] = face.varₗ[👉.v]
        
        face.varₗ[👉.w] = cells[face.owner].var[👉.w]
        face.varᵣ[👉.w] = face.varₗ[👉.w]
        
        face.varₗ[👉.T] = cells[face.owner].var[👉.T]
        face.varᵣ[👉.T] = face.varₗ[👉.T]
        
        face.varₗ[👉.Y₁] = cells[face.owner].var[👉.Y₁]
        face.varᵣ[👉.Y₁] = face.varₗ[👉.Y₁]

    end


end
    

