function update_pseudo_conservative!(
    👉::controls, 
    cells::Vector{mesh.Cell}
)
    for cell in cells
        cell.Qᵐ[1] = cell.var[👉.ρ]
        cell.Qᵐ[2] = cell.var[👉.ρ]*cell.var[👉.u]
        cell.Qᵐ[3] = cell.var[👉.ρ]*cell.var[👉.v]
        cell.Qᵐ[4] = cell.var[👉.ρ]*cell.var[👉.w]
        cell.Qᵐ[5] = cell.var[👉.ρ]*cell.var[👉.Hₜ]-cell.var[👉.p]
        cell.Qᵐ[6] = cell.var[👉.ρ]*cell.var[👉.Y₁]
    end

end



function update_real_conservative!(
    👉::controls, 
    cells::Vector{mesh.Cell}
)

    for cell in cells
        
        cell.Qⁿ⁻¹[:] = cell.Qⁿ[:]

        cell.Qⁿ[1] = cell.var[👉.ρ]
        cell.Qⁿ[2] = cell.var[👉.ρ]*cell.var[👉.u]
        cell.Qⁿ[3] = cell.var[👉.ρ]*cell.var[👉.v]
        cell.Qⁿ[4] = cell.var[👉.ρ]*cell.var[👉.w]
        cell.Qⁿ[5] = cell.var[👉.ρ]*cell.var[👉.Hₜ]-cell.var[👉.p]
        cell.Qⁿ[6] = cell.var[👉.ρ]*cell.var[👉.Y₁]
    end

end




function update_real_conservative0!(
    👉::controls, 
    cells::Vector{mesh.Cell}
)

    for cell in cells
        
        cell.Qⁿ⁻¹[1] = cell.var[👉.ρ]
        cell.Qⁿ⁻¹[2] = cell.var[👉.ρ]*cell.var[👉.u]
        cell.Qⁿ⁻¹[3] = cell.var[👉.ρ]*cell.var[👉.v]
        cell.Qⁿ⁻¹[4] = cell.var[👉.ρ]*cell.var[👉.w]
        cell.Qⁿ⁻¹[5] = cell.var[👉.ρ]*cell.var[👉.Hₜ]-cell.var[👉.p]
        cell.Qⁿ⁻¹[6] = cell.var[👉.ρ]*cell.var[👉.Y₁]

        cell.Qⁿ[:] = cell.Qⁿ⁻¹[:]
    end

end
