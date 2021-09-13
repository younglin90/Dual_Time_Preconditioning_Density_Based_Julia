function linear_solver_explicit!(
    A::Array{Float64},
    x::Array{Float64},
    B::Array{Float64}
)

#@time begin
    n = length(A[:, 1, 1])
    for i in 1:n
        #invA = inv(A[i, :, :])
        #x[i, :] = invA[:, :]*B[i, :]
        
        x[i, :] = A[i, :, :]\B[i, :]
    end
#end

end





function linear_solver_implicit!(
    👉::controls, 
    cell::Vector{mesh.Cell},
    face::Vector{mesh.Face}
)


end