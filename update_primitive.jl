function update_primitive!(
    👉::controls, 
    x::Array{Float64},
    cells::Vector{mesh.Cell}
)

    for i in 1:length(cells)
        cells[i].var[👉.p] += x[i, 1]
        cells[i].var[👉.u] += x[i, 2]
        cells[i].var[👉.v] += x[i, 3]
        cells[i].var[👉.w] += x[i, 4]
        cells[i].var[👉.T] += x[i, 5]
        cells[i].var[👉.Y₁] += x[i, 6]
    end

end
