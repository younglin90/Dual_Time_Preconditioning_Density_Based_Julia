function linear_solver_explicit!(
    A::Array{Float64},
    x::Array{Float64},
    B::Array{Float64}
)

    n = length(A[:, 1, 1])
    for i in 1:n
        #println(A[i, :, :])

        #T = A[i, :, :]
        invA = inv(A[i, :, :])
        #println(invA)


        x[i, :] = invA[:, :]*B[i, :]
    end


end





function linear_solver_implicit!(
    👉::controls, 
    cell::Vector{mesh.Cell},
    face::Vector{mesh.Face}
)


end