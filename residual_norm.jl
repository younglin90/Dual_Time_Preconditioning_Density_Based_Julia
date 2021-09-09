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
        end
    end
    return log10(RMS/(maxQ+1.e-200))

end
