include("./structured_grid_uniform.jl")
include("./constant.jl")
include("./controls.jl")
include("./timestep.jl")
include("./construct_A_matrix.jl")
include("./face_left_right.jl")
include("./flux.jl")
include("./linear_solver.jl")
include("./update_conservative.jl")
include("./update_primitive.jl")
include("./real_time_terms.jl")
include("./residual_norm.jl")
include("./EOS.jl")

using Plots

function main()


        #🎲
        #⬜
        #◽

    Nx = 160
    Ny = 1
    Nz = 1
    Lx = 1.0
    Ly = 0.1
    Lz = 0.1
    realMaxIter = 1000000
    pseudoMaxIter = 30
    pseudoMaxResidual = -50.0

    CFL = 0.5
    Δt = 1.e-6
    Lco = 1.0
    Uco = 1.0

    👉 = controls(
        Nx,Ny,Nz, Lx,Ly,Lz, 
        realMaxIter,pseudoMaxIter,pseudoMaxResidual, 
        CFL, Δt, Lco, Uco,
        0.0, 0, 0, 0.0, 
        1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20
    )

    cells = Vector{mesh.Cell}(undef, 0)
    faces = Vector{mesh.Face}(undef, 0)
    faces_internal = Vector{mesh.Face}(undef, 0)
    faces_boundary = Vector{mesh.Face}(undef, 0)
    faces_boundary_top = Vector{mesh.Face}(undef, 0)
    faces_boundary_bottom = Vector{mesh.Face}(undef, 0)
    faces_boundary_left = Vector{mesh.Face}(undef, 0)
    faces_boundary_right = Vector{mesh.Face}(undef, 0)

    structured_grid_uniform!(
        👉,
        cells,
        faces,
        faces_internal,
        faces_boundary,
        faces_boundary_top,
        faces_boundary_bottom,
        faces_boundary_left,
        faces_boundary_right
    )

    println(" Cell size = ",length(cells))
    println(" Face size = ",length(faces))
    println(" Face internal size = ",length(faces_internal))
    println(" Face boundary size = ",length(faces_boundary))
    println(" Face boundary top size = ",length(faces_boundary_top))
    println(" Face boundary bottom size = ",length(faces_boundary_bottom))
    println(" Face boundary left size = ",length(faces_boundary_left))
    println(" Face boundary right size = ",length(faces_boundary_right))

    # initialization
    for cell in cells
        cell.var[👉.p] = 101325.0
        cell.var[👉.u] = 0.0
        cell.var[👉.v] = 0.0
        cell.var[👉.w] = 0.0
        cell.var[👉.T] = 300.0
        cell.var[👉.Y₁] = 0.0
    end

    half_cell_num::Int32 = round(length(cells)/2)
    for i in 1:half_cell_num
        cells[i].var[👉.p] = 1.e6
        cells[i].var[👉.u] = -2500.0
        cells[i].var[👉.v] = 0.0
        cells[i].var[👉.w] = 0.0
        cells[i].var[👉.T] = 1000.0
        cells[i].var[👉.Y₁] = 0.0
    end
    
    for i in half_cell_num+1:length(cells)
        cells[i].var[👉.p] = 1.e6
        cells[i].var[👉.u] = 2500.0
        cells[i].var[👉.v] = 0.0
        cells[i].var[👉.w] = 0.0
        cells[i].var[👉.T] = 950.0
        cells[i].var[👉.Y₁] = 0.0
    end

    # EOS
    EOS!(👉, cells)


    # solver
    👉.time = 0.0
    👉.realIter = 1
    while(
        👉.realIter <= 👉.realMaxIter
    )
        println("real-time Step: $(👉.realIter) \t Time: $(round((👉.Δt),digits=3))")

        # Qⁿ, Qⁿ⁻¹
        if 👉.realIter == 1
            update_real_conservative0!(👉, cells)
        else
            update_real_conservative!(👉, cells)
        end

        👉.pseudoIter = 1
        👉.residual = 10000.0
        residual0 = 10000.0
        while(
            👉.pseudoIter <= 👉.pseudoMaxIter &&
            👉.residual >= 👉.pseudoMaxResidual
        )

            if 👉.pseudoIter == 1
                CFL = 0.01
            else
                CFL = 1.0
            end

            # time-step
            timestep!(👉, cells)

            # face left, Right
            face_left_right!(👉, cells, faces, faces_internal, faces_boundary, 
                faces_boundary_top, faces_boundary_bottom, faces_boundary_left, faces_boundary_right)

            # right hand side
            B = zeros(Float64, length(cells), 6)

            # flux
            flux!(👉, B, cells, faces_internal, faces_boundary)

            # Qᵐ
            update_pseudo_conservative!(👉, cells)

            # real time terms
            real_time_terms!(👉, B, cells)

            # sparse A matrix
            A = zeros(Float64, length(cells), 6, 6)
            #construct_A_matrix_implicit!(👉, A, cells, faces)
            construct_A_matrix_explicit!(👉, A, cells)

            # Ax=B
            x = zeros(Float64, length(cells), 6)
            #linear_solver_implicit!(A, x, B)
            linear_solver_explicit!(A, x, B)

            # residual norm
            👉.residual = residual_norm!(x, cells)
            if 👉.pseudoIter == 1
                residual0 = 👉.residual
            end
            
            # update primitive
            update_primitive!(👉, x, cells)
            
            # EOS
            EOS!(👉, cells)
            
        
            println("- pseudo-time Step: $(👉.pseudoIter) \t",
            "log₁₀|ΔR|₂: $(round((👉.residual-residual0),digits=8))")

            gr()
            X = zeros(Float64, length(cells), 1)
            Y = zeros(Float64, length(cells), 6)
            for i in 1:length(cells)
                X[i] = cells[i].x
                Y[i,1] = cells[i].var[👉.p]
                Y[i,2] = cells[i].var[👉.u]
                Y[i,3] = cells[i].var[👉.T]
                Y[i,4] = cells[i].var[👉.Y₁]
                Y[i,5] = cells[i].var[👉.ρ]
                Y[i,6] = cells[i].var[👉.c]
            end
            plt = plot(X,Y,layout = (3, 2),label = ["p" "u" "T" "Y₁" "ρ" "c"] )
            gui(); sleep(0.5)


            👉.pseudoIter += 1

        end

        👉.realIter += 1
        👉.time += 👉.Δt

    end





    #Δt = CFL * Δx
end



# calculation main
main()


