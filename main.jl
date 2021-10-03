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
#using PlotlyJS
using LinearAlgebra
using SparseArrays
using IterativeSolvers
#using AlgebraicMultigrid

using Pardiso


function plotting1D(
    Nx, Ny, 
    👉::controls,
    cells::Vector{mesh.Cell}
)

    X = zeros(Float64, length(cells), 1)
    Y = zeros(Float64, length(cells), 8)
    for i in 1:length(cells)
        X[i] = cells[i].x
        Y[i,1] = cells[i].var[👉.p]
        Y[i,2] = cells[i].var[👉.u]
        Y[i,3] = cells[i].var[👉.v]
        Y[i,4] = cells[i].var[👉.T]
        Y[i,5] = cells[i].var[👉.α₁]
        Y[i,6] = cells[i].var[👉.ρ]
        Y[i,7] = cells[i].var[👉.c]
        Y[i,8] = cells[i].var[👉.Hₜ]
        
    end
    plot(X,Y,layout = grid(3, 3), label = ["p" "u" "v" "T" "α₁" "ρ" "c" "Hₜ"] )

    gui()

end


function plotting2D(
    Nx, Ny, 
    👉::controls,
    cells::Vector{mesh.Cell}
)

    #plt = plot(X,VAR,layout = 
    #grid(3, 2),
    #label = ["p" "u" "T" "Y₁" "ρ" "c"] )
    #plot(plt)
    #contourf!(X,Y,VAR)

    X = zeros(Float64, Nx)
    Y = zeros(Float64, Ny)
    VAR1 = zeros(Float64, Nx, Ny)
    VAR2 = zeros(Float64, Nx, Ny)
    VAR3 = zeros(Float64, Nx, Ny)
    VAR4 = zeros(Float64, Nx, Ny)
    VAR5 = zeros(Float64, Nx, Ny)
    VAR6 = zeros(Float64, Nx, Ny)
    for i in 1:Nx
        for j in 1:Ny
            k=1
            ijk = i + Nx*(j-1) + Nx*Ny*(k-1)
            X[i] = cells[ijk].x
            Y[j] = cells[ijk].y
            VAR1[i,j] = cells[ijk].var[👉.p]
            VAR2[i,j] = cells[ijk].var[👉.ρ]
            #VAR2[i,j] = cells[ijk].var[👉.α₁]
            VAR3[i,j] = cells[ijk].var[👉.u]
            VAR4[i,j] = cells[ijk].var[👉.v]
            VAR5[i,j] = cells[ijk].var[👉.w]
            VAR6[i,j] = cells[ijk].var[👉.T]

        end
    end

    #plotlyjs()
    #X = 0.5*Δx:Δx:👉.Lx
    #Y = 0.5*Δy:Δy:👉.Ly
    #X = repeat(reshape(x, 1, :), length(y), 1)
    #Y = repeat(y, 1, length(x))
    #plot(contour(X, Y, VAR2, fill = true))
    plot(
        heatmap(X, Y, VAR1', c = :bluesreds),
        heatmap(X, Y, VAR2', c = :bluesreds),
        heatmap(X, Y, VAR3', c = :bluesreds),
        heatmap(X, Y, VAR4', c = :bluesreds),
        heatmap(X, Y, VAR5', c = :bluesreds),
        heatmap(X, Y, VAR6', c = :bluesreds);
        layout = grid(3, 2)
    )

    gui()
#=
    plot(contour(
        x=0.5*Δx:Δx:👉.Lx,#X, # horizontal axis
        y=0.5*Δy:Δy:👉.Ly,#Y, # vertical axis
        z=VAR2'#VAR[:,5]'
    ))
=#



end



function main()


        #🎲
        #⬜
        #◽

    Nx = 50
    Ny = 50
    Nz = 1
    Lx = 1.0
    Ly = 1.0
    Lz = 0.1
    realMaxIter = 1000000
    pseudoMaxIter = 30
    pseudoMaxResidual = -10000000.0

    CFL = 1.0
    Δt = 1.e-3
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



    # dam break
    for cell in cells
        cell.var[👉.p] = 101325.0
        cell.var[👉.u] = 0.0
        cell.var[👉.v] = 0.0
        cell.var[👉.w] = 0.0
        cell.var[👉.T] = 300.0
        cell.var[👉.Y₁] = 0.0
        cell.var[👉.α₁] = 0.0

        if cell.x < 0.4 && cell.y < 0.4
            cell.var[👉.Y₁] = 1.0
            cell.var[👉.α₁] = 1.0
        end
    end




#=
    # 1D cavitation
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
=#

#=

    # 1D high pressure water & low pressure air
    for cell in cells
        if cell.x < 0.7
            cell.var[👉.p] = 1.e9
            cell.var[👉.u] = 0.0
            cell.var[👉.v] = 0.0
            cell.var[👉.w] = 0.0
            cell.var[👉.T] = 300.0
            cell.var[👉.Y₁] = 1.0
            cell.var[👉.α₁] = 1.0
        else
            cell.var[👉.p] = 1.e5
            cell.var[👉.u] = 0.0
            cell.var[👉.v] = 0.0
            cell.var[👉.w] = 0.0
            cell.var[👉.T] = 6.968
            cell.var[👉.Y₁] = 0.0
            cell.var[👉.α₁] = 0.0
        end
    end

=#

#=
    # 1D high pressure water & low pressure air
    for cell in cells
        if cell.x < 0.5
            cell.var[👉.p] = 1.0
            cell.var[👉.u] = 0.0
            cell.var[👉.v] = 0.0
            cell.var[👉.w] = 0.0
            cell.var[👉.T] = 0.003484
            cell.var[👉.Y₁] = 0.0
            cell.var[👉.α₁] = 0.0
        else
            cell.var[👉.p] = 0.1
            cell.var[👉.u] = 0.0
            cell.var[👉.v] = 0.0
            cell.var[👉.w] = 0.0
            cell.var[👉.T] = 0.002787
            cell.var[👉.Y₁] = 0.0
            cell.var[👉.α₁] = 0.0
        end
    end
=#


#=
    # One-dimensional helium-bubble in air
    for cell in cells
        
        if cell.x < 0.3
            cell.var[👉.p] = 1.245e5
            cell.var[👉.u] = 55.33
            cell.var[👉.v] = 0.0
            cell.var[👉.w] = 0.0
            cell.var[👉.T] = 319.48
            cell.var[👉.Y₁] = 0.0
            cell.var[👉.α₁] = 0.0
        else
            cell.var[👉.p] = 1.e5
            cell.var[👉.u] = 0.0
            cell.var[👉.v] = 0.0
            cell.var[👉.w] = 0.0
            cell.var[👉.T] = 300.0
            cell.var[👉.Y₁] = 0.0
            cell.var[👉.α₁] = 0.0
        end

        if 0.5 < cell.x < 0.7
            cell.var[👉.Y₁] = 1.0
            cell.var[👉.α₁] = 1.0
        end
    end
=#

    # EOS
    EOS!(👉, cells)


    # solver
    👉.time = 0.0
    👉.realIter = 1
    total_iter = 1.0
    
    plt2 = plot([],[])

    while(
        👉.realIter <= 👉.realMaxIter
    )
        println("real-time Step: $(👉.realIter) \t Time: $(👉.time)")

        # Qⁿ, Qⁿ⁻¹
        if 👉.realIter == 1
            update_real_conservative0!(👉, cells)
        else
            update_real_conservative!(👉, cells)
        end

        if 👉.realIter < 3
            👉.CFL = 0.01
            👉.pseudoMaxIter = 15
        else
            👉.CFL = 0.01
            👉.pseudoMaxIter = 15
        end

        👉.pseudoIter = 1
        👉.residual = 10000.0
        residual0 = 10000.0
        while(
            👉.pseudoIter ≤ 👉.pseudoMaxIter &&
            👉.residual-residual0 ≥ 👉.pseudoMaxResidual
        )
        #=
            if 👉.pseudoIter == 1
                👉.CFL = 0.0000001
            else
                👉.CFL = 0.1
            end
        =#
            # time-step
            timestep!(👉, cells, Lx/Nx, Ly/Ny, Lz/Nz)

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
            #construct_A_matrix_explicit!(👉, A, cells)

            👉.residual = 
            #construct_A_matrix_implicit!(
            construct_A_matrix_explicit!(
                👉, 
                cells,
                faces,
                faces_internal,
                faces_boundary,
                faces_boundary_top,
                faces_boundary_bottom,
                faces_boundary_left,
                faces_boundary_right,
                B
            )
                

            #=
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
            =#
            
            # EOS
            EOS!(👉, cells)
            
        
            #println("- pseudo-time Step: $(👉.pseudoIter) \t",
            #"log₁₀|ΔR|₂: $(round((👉.residual-residual0),digits=8))")
            println("- pseudo-time Step: $(👉.pseudoIter) \t",
            "log₁₀|ΔR|₂: $(round((👉.residual),digits=8))")

            plotting2D(Nx, Ny, 👉, cells)

            👉.pseudoIter += 1
            total_iter += 1.0

        end

        👉.realIter += 1
        👉.time += 👉.Δt

    end





    #Δt = CFL * Δx
end



# calculation main
main()


