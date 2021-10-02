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

        cells[i].var[👉.p] = max(cells[i].var[👉.p],1.e-200)
        cells[i].var[👉.T] = max(cells[i].var[👉.T],1.e-200)
    end

end
