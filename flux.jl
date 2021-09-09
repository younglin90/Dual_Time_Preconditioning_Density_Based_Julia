function flux!(
    👉::controls, 
    RHS::Array{Float64},
    cells::Vector{mesh.Cell},
    faces_internal::Vector{mesh.Face},
    faces_boundary::Vector{mesh.Face}
)

    for face in faces_internal
        flux = zeros(Float64, 6, 1)
        flux = AUSMPWP_N(
            face.varₗ[👉.p],face.varₗ[👉.u],face.varₗ[👉.v],face.varₗ[👉.w],
            face.varₗ[👉.T],face.varₗ[👉.Y₁],face.varₗ[👉.ρ],face.varₗ[👉.Hₜ],face.varₗ[👉.c],
            face.varᵣ[👉.p],face.varᵣ[👉.u],face.varᵣ[👉.v],face.varᵣ[👉.w],
            face.varᵣ[👉.T],face.varᵣ[👉.Y₁],face.varᵣ[👉.ρ],face.varᵣ[👉.Hₜ],face.varᵣ[👉.c],
            👉.Lco,👉.Uco,👉.Δt,
            1.0, 1.0,
            face.n̂[1],face.n̂[2],face.n̂[3]
            )

        RHS[face.owner, :] -= flux[:]*face.ΔS/cells[face.owner].Ω
        RHS[face.neighbour, :] += flux[:]*face.ΔS/cells[face.neighbour].Ω

    end

    for face in faces_boundary

        flux = zeros(Float64, 6, 1)
        flux = AUSMPWP_N(
            face.varₗ[👉.p],face.varₗ[👉.u],face.varₗ[👉.v],face.varₗ[👉.w],
            face.varₗ[👉.T],face.varₗ[👉.Y₁],face.varₗ[👉.ρ],face.varₗ[👉.Hₜ],face.varₗ[👉.c],
            face.varᵣ[👉.p],face.varᵣ[👉.u],face.varᵣ[👉.v],face.varᵣ[👉.w],
            face.varᵣ[👉.T],face.varᵣ[👉.Y₁],face.varᵣ[👉.ρ],face.varᵣ[👉.Hₜ],face.varᵣ[👉.c],
            👉.Lco,👉.Uco,👉.Δt,
            1.0, 1.0,
            face.n̂[1],face.n̂[2],face.n̂[3]
            )

        RHS[face.owner, :] -= flux[:]*face.ΔS/cells[face.owner].Ω

    end


end




function AUSMPWP_N(
    pₗ,uₗ,vₗ,wₗ,Tₗ,Y₁ₗ,ρₗ,Hₜₗ,cₗ,
    pᵣ,uᵣ,vᵣ,wᵣ,Tᵣ,Y₁ᵣ,ρᵣ,Hₜᵣ,cᵣ,
    Lco,Uco,Δt,
    w₁,w₂,
    nx,ny,nz
)


	# properties of Left
	Uₙₗ = uₗ*nx + vₗ*ny + wₗ*nz
	# properties of Right
	Uₙᵣ = uᵣ*nx + vᵣ*ny + wᵣ*nz

	RT = √(ρᵣ/ρₗ)
    ρ̄ = (ρₗ+RT*ρᵣ)/(1.0+RT)
    ū = (uₗ+RT*uᵣ)/(1.0+RT)
    v̄ = (vₗ+RT*vᵣ)/(1.0+RT)
    w̄ = (wₗ+RT*wᵣ)/(1.0+RT)
	# Low-Diffusion Flux-Splitting Methods for Real Fluid Flows at All Speeds, 1999
	c̄ = 0.5*(cₗ + cᵣ)
	Mco = Uco/c̄
	M̄ₖ = √(ū^2 + v̄^2 + w̄^2)/c̄
	Mun = Lco/π/Δt/c̄
	# Scaling ============
	θₜ = max(M̄ₖ,Mco)
	θₚ = min(1.0,max(θₜ,Mun))
	θᵤ = min(1.0,θₜ)
	ϕₚ = θₚ*(2.0-θₚ)
	ϕᵤ = θᵤ*(2.0-θᵤ)
	aSDST = 0.1875*(-4.0+5.0*ϕᵤ*ϕᵤ)
	Mₗ = Uₙₗ/c̄
	Mᵣ = Uₙᵣ/c̄
	# calculate M+ and P+ for left state
	Mₗ⁺ = M_func(Mₗ,1.0,0.125)
	p⁺ = pre_func(Mₗ,1.0,aSDST)
	# calculate M- and P- for left state
	Mᵣ⁻ = M_func(Mᵣ,-1.0,0.125)
	p⁻ = pre_func(Mᵣ,-1.0,aSDST)
	# ======= carefully, sensitive ===========
	Mₗᵣ = Mₗ⁺ + Mᵣ⁻
	ρₗᵣ = ( Mₗᵣ > 0.0 ? ρₗ : ρᵣ )
	w₁ = 1.0 - w₁^3
	w₂ = 1.0 - w₂^2
	w = max(w₁,w₂)
	fₗ = pₗ/(ρ̄ *c̄^2)*(1.0-w)*ρ̄ /ρₗᵣ   /ϕₚ
	fᵣ = pᵣ/(ρ̄ *c̄^2)*(1.0-w)*ρ̄ /ρₗᵣ   /ϕₚ
	# =======================================
	M̄⁺ = 0.0
    M̄⁻ = 0.0
	if Mₗᵣ >= 0.0
        M̄⁺ = Mₗ⁺ + Mᵣ⁻ * ((1.0-w)*(1.0+fᵣ)-fₗ)
        M̄⁻ = Mᵣ⁻ * w * (1.0+fᵣ)
	else
        M̄⁺ = Mₗ⁺ * w * (1.0+fₗ)
        M̄⁻ = Mᵣ⁻ + Mₗ⁺ * ((1.0-w)*(1.0+fₗ)-fᵣ)
    end
	ṁₗ = c̄ * ρₗ * M̄⁺
	ṁᵣ = c̄ * ρᵣ * M̄⁻
	Kᵤ = 0.5
	pᵤ = -2.0 * Kᵤ * p⁺ * p⁻ * ρ̄  * c̄ * (Uₙᵣ-Uₙₗ)
	pₗᵣ = p⁺*pₗ + p⁻*pᵣ +  pᵤ * ϕᵤ
	
	# comp. convective flux
    flux = zeros(Float64,6,1)
	flux[1] = ṁₗ + ṁᵣ
	flux[2] = ṁₗ*uₗ + ṁᵣ*uᵣ + pₗᵣ*nx
	flux[3] = ṁₗ*vₗ + ṁᵣ*vᵣ + pₗᵣ*ny
	flux[4] = ṁₗ*wₗ + ṁᵣ*wᵣ + pₗᵣ*nz
	flux[5] = ṁₗ*Hₜₗ + ṁᵣ*Hₜᵣ
	flux[6] = ṁₗ*Y₁ₗ + ṁᵣ*Y₁ᵣ

    return flux


end



function M_func(M::Float64, op::Float64, α::Float64)
    mu = 0.0
	if abs(M) > 1.0 
		mu = 0.5*(M + op*abs(M))
	else
		mu = op*0.25*(M + op)^2.0 + op*α*(M*M-1.0)^2.0
    end
	
	return mu
end

function pre_func(M::Float64, op::Float64, α::Float64)
    mu = 0.0
	if abs(M) > 1.0
		mu = 0.5*(1.0 + op*sign(M) )
	else
		mu = 0.25*(M + op)^2.0*(2.0-op*M) + op*α*M*(M*M-1.0)^2.0
    end
	
	return mu;
end

