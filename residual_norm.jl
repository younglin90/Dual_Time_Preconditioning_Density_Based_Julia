function residual_norm!(
    x::Array{Float64},
    cells::Vector{mesh.Cell}
)
    RMS::Float64 = 0.0
    maxQ::Float64 = -100000.0
    for i in 1:length(cells)
        for j in 1:6
            RMS += x[i, j]^2
            maxQ = max(maxQ, abs(x[i, j]))
            #=
            maxQ = max(maxQ, abs(cells[i].var[1]))
            maxQ = max(maxQ, abs(cells[i].var[2]))
            maxQ = max(maxQ, abs(cells[i].var[3]))
            maxQ = max(maxQ, abs(cells[i].var[4]))
            maxQ = max(maxQ, abs(cells[i].var[5]))
            maxQ = max(maxQ, abs(cells[i].var[6]))
            =#
        end
    end
    return log10(√(RMS/6.0)/(maxQ+1.e-200))

end
