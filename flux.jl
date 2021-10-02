function flux!(
    👉::controls, 
    RHS::Array{Float64},
    cells::Vector{mesh.Cell},
    faces_internal::Vector{mesh.Face},
    faces_boundary::Vector{mesh.Face}
)


	
	shock_sensor = shock_discontinuity_sensing_term!(
		👉,cells,faces_internal,faces_boundary
	)



    for face in faces_internal

		cpi2 = min(shock_sensor[face.owner],shock_sensor[face.neighbour])
		w₁ = cpi2
		pₗ = face.varₗ[👉.p]
		ρₗ = face.varₗ[👉.ρ]
		cₗ = face.varₗ[👉.c]
		pᵣ = face.varᵣ[👉.p]
		ρᵣ = face.varᵣ[👉.ρ]
		cᵣ = face.varᵣ[👉.c]
		preLs = abs(pₗ) + 0.1 * ρₗ*cₗ*cₗ
		preRs = abs(pᵣ) + 0.1 * ρᵣ*cᵣ*cᵣ
		w₂ = min(preLs/preRs,preRs/preLs)

        flux = zeros(Float64, 6, 1)
        flux = YYL_Riemann(
            face.varₗ[👉.p],face.varₗ[👉.u],face.varₗ[👉.v],face.varₗ[👉.w],
            face.varₗ[👉.T],face.varₗ[👉.Y₁],face.varₗ[👉.ρ],face.varₗ[👉.Hₜ],face.varₗ[👉.c],
            face.varᵣ[👉.p],face.varᵣ[👉.u],face.varᵣ[👉.v],face.varᵣ[👉.w],
            face.varᵣ[👉.T],face.varᵣ[👉.Y₁],face.varᵣ[👉.ρ],face.varᵣ[👉.Hₜ],face.varᵣ[👉.c],
            👉.Lco,👉.Uco,👉.Δt,
            w₁, w₂, cpi2,
            face.n̂[1],face.n̂[2],face.n̂[3]
            )


        RHS[face.owner, :] -= flux[:]*face.ΔS/cells[face.owner].Ω #* 👉.Δt
        RHS[face.neighbour, :] += flux[:]*face.ΔS/cells[face.neighbour].Ω #* 👉.Δt

    end


    for face in faces_boundary

        flux = zeros(Float64, 6, 1)
        flux = YYL_Riemann(
            face.varₗ[👉.p],face.varₗ[👉.u],face.varₗ[👉.v],face.varₗ[👉.w],
            face.varₗ[👉.T],face.varₗ[👉.Y₁],face.varₗ[👉.ρ],face.varₗ[👉.Hₜ],face.varₗ[👉.c],
            face.varᵣ[👉.p],face.varᵣ[👉.u],face.varᵣ[👉.v],face.varᵣ[👉.w],
            face.varᵣ[👉.T],face.varᵣ[👉.Y₁],face.varᵣ[👉.ρ],face.varᵣ[👉.Hₜ],face.varᵣ[👉.c],
            👉.Lco,👉.Uco,👉.Δt,
            1.0, 1.0, 1.0,
            face.n̂[1],face.n̂[2],face.n̂[3]
            )

        RHS[face.owner, :] -= flux[:]*face.ΔS/cells[face.owner].Ω #* 👉.Δt

    end

	i = 1
	for cell in cells

		RHS[i,3] += cell.var[👉.ρ]*(-9.8)

		i += 1
	end
	

end

function shock_discontinuity_sensing_term!(
    👉::controls, 
    cells::Vector{mesh.Cell},
    faces_internal::Vector{mesh.Face},
    faces_boundary::Vector{mesh.Face}
)


	udf_coeff_shock = 0.1

	
    cell_cpi2 = zeros(Float64,length(cells))

	for face in faces_internal

		pL = face.varₗ[👉.p]
		rhoL = face.varₗ[👉.ρ]
		cL = face.varₗ[👉.c]

		#/==========

		pR = face.varᵣ[👉.p]
		rhoR = face.varᵣ[👉.ρ]
		cR = face.varᵣ[👉.c]

		#/===========
		#RT = sqrt(rhoR/rhoL)
		pLs = pL + udf_coeff_shock * min( rhoL*cL*cL, rhoR*cR*cR ) 
		pRs = pR + udf_coeff_shock * min( rhoL*cL*cL, rhoR*cR*cR ) 

		cpi2 = min( pLs/pRs, pRs/pLs )

		cell_cpi2[face.owner] = min( cell_cpi2[face.owner] , cpi2 )
		cell_cpi2[face.neighbour] = min( cell_cpi2[face.neighbour] , cpi2 )

	end

	return cell_cpi2



end

function YYL_Riemann(
    pₗ::Float64,uₗ::Float64,vₗ::Float64,wₗ::Float64,
	Tₗ::Float64,Y₁ₗ::Float64,ρₗ::Float64,Hₜₗ::Float64,cₗ::Float64,
    pᵣ::Float64,uᵣ::Float64,vᵣ::Float64,wᵣ::Float64,
	Tᵣ::Float64,Y₁ᵣ::Float64,ρᵣ::Float64,Hₜᵣ::Float64,cᵣ::Float64,
    Lco::Float64,Uco::Float64,Δt::Float64,
    w₁::Float64,w₂::Float64,cpi2::Float64,
    nx::Float64,ny::Float64,nz::Float64
)

	Uₙₗ = uₗ*nx + vₗ*ny + wₗ*nz
	Uₙᵣ = uᵣ*nx + vᵣ*ny + wᵣ*nz

	c̄ = 0.5*(cₗ + cᵣ)
	#M̄ₖ = √(ū^2 + v̄^2 + w̄^2)/c̄
	Mₗ = Uₙₗ/c̄
	Mᵣ = Uₙᵣ/c̄
	# calculate M+ and P+ for left state
	Mₗ⁺ = M_func(Mₗ,1.0,0.125)
	p⁺ = pre_func(Mₗ,1.0,0.1875)
	# calculate M- and P- for left state
	Mᵣ⁻ = M_func(Mᵣ,-1.0,0.125)
	p⁻ = pre_func(Mᵣ,-1.0,0.1875)
	M̄ = (ρₗ*abs(Mₗ)+ρᵣ*abs(Mᵣ)) / (ρₗ+ρᵣ)
	KLR = sqrt(0.5*(uₗ^2+vₗ^2+wₗ^2+uᵣ^2+vᵣ^2+wᵣ^2))
	g = 1.0 + max(min(Mₗ,0.0),-1.0)*min(max(Mᵣ,0.0),1.0)
	D_L = Mₗ + (1.0-g)*abs(Mₗ)
	D_R = Mᵣ - (1.0-g)*abs(Mᵣ)
	D_rho = M̄*g

	preLs = abs(pₗ) + 0.1 * ρₗ*cₗ*cₗ
	preRs = abs(pᵣ) + 0.1 * ρᵣ*cᵣ*cᵣ
	w = 1.0 - min(preLs/preRs,preRs/preLs)^2.0

	ps = p⁺*pₗ+p⁻*pᵣ
	pll = 0.0
	if 3.0/4.0 <= min(pₗ/pᵣ,pᵣ/pₗ) && 1.0 > min(pₗ/pᵣ,pᵣ/pₗ) 
		pll=4.0*min(pₗ/pᵣ,pᵣ/pₗ)-3.0
	end
	fL = 0.0
	if abs(Mₗ) <= 1.0 
		fL = (pₗ/ps-1.0)*pll*abs(Mₗ⁺)*min(1.0,( (abs(Uₙₗ)/c̄) )^0.25)
	end
		
	fR = 0.0
	if abs(Mᵣ) <= 1.0
		fR = (pₗ/ps-1.0)*pll*abs(Mᵣ⁻)*min(1.0,( (abs(Uₙᵣ)/c̄) )^0.25)
	end

	MPP = Mₗ⁺+Mᵣ⁻
	MLP_AUSM = 0.5*(MPP+abs(MPP))
	MRM_AUSM = 0.5*(MPP-abs(MPP))

	MLP_SLAU = 0.5*(D_L+D_rho)
	MRM_SLAU = 0.5*(D_R-D_rho)

	fa1 = w

	MLPL = fa1*MLP_AUSM + (1.0-fa1)*MLP_SLAU
	MRMR = fa1*MRM_AUSM + (1.0-fa1)*MRM_SLAU

	ṁ = ρₗ*c̄*MLPL + ρᵣ*c̄*MRMR - 
		0.5*
		(1.0-0.5*(1.0-cos(3.141592*min(1.0,max(abs(Mₗ),abs(Mᵣ))))))*
		(1.0-0.5*(1.0+cos(3.141592*min(abs(pₗ/pᵣ),abs(pᵣ/pₗ)))))*
		(pᵣ-pₗ)/c̄

	Mₗ⁺ = 0.5*(Mₗ+abs(Mₗ))
	Mᵣ⁻ = 0.5*(Mᵣ-abs(Mᵣ))
	ṁₗ = 0.0
	ṁᵣ = 0.0
    if ṁ >= 0.0
		ṁₗ = ṁ - (ρₗ*c̄*Mᵣ⁻)*( w*(1.0+fR)-fR+fL )
		ṁᵣ = (ρᵣ*c̄*Mᵣ⁻)*( w*(1.0+fR) )
    else
		ṁₗ = (ρₗ*c̄*Mₗ⁺)*( w*(1.0+fL) )
		ṁᵣ = ṁ - (ρᵣ*c̄*Mₗ⁺)*( w*(1.0+fL)-fL+fR )
	end
	

	fa2 = w
	gam = 0.6

	pₗᵣ = 0.5*(pₗ+pᵣ) - 
			fa2*(KLR/c̄)*0.5*p⁺*p⁻*0.5*(pₗ+pᵣ)/c̄*(Uₙᵣ-Uₙₗ) + 
			max(0.2,gam)*(KLR/c̄)*0.5*(pₗ+pᵣ)*(p⁺+p⁻-1.0) - 
			0.5*(p⁺-p⁻)*(pᵣ-pₗ)

	#pₗᵣ = 0.5*(pₗ+pᵣ) + KLR*0.5*(ρₗ+ρᵣ)*c̄*(p⁺+p⁻-1.0) - 0.5*(p⁺-p⁻)*(pᵣ-pₗ)

	#pₗᵣ = pₗ*p⁺ + pᵣ*p⁻
	#pₗᵣ = 0.5*(pₗ + pᵣ)
	#=
	ṁ = ρₗ*c̄*MLP_SLAU + ρᵣ*c̄*MRM_SLAU - 
			0.5*(1.0-min(1.0,1.0/c̄*KLR))^2.0 *(pᵣ-pₗ)/c̄
	if ṁ >= 0.0
		ṁₗ = ṁ #- (ρₗ*c̄*Mᵣ⁻)*( w*(1.0+fR)-fR+fL )
		ṁᵣ = 0.0# (ρᵣ*c̄*Mᵣ⁻)*( w*(1.0+fR) )
	else
		ṁₗ = 0.0#(ρₗ*c̄*Mₗ⁺)*( w*(1.0+fL) )
		ṁᵣ = ṁ #- (ρᵣ*c̄*Mₗ⁺)*( w*(1.0+fL)-fL+fR )
	end
	=#
#=
	UU = 0.5*(Uₙₗ+Uₙᵣ)
	if UU > 0.0
		ṁₗ = ρₗ * UU
		ṁᵣ = 0.0
	else
		ṁₗ = 0.0
		ṁᵣ = ρᵣ * UU
	end
=#
	#rhohat = 0.5*(ρₗ+ρᵣ)
    #gam2 = 0.5*(1.0-tanh(15.0*(min(abs(pₗ/pᵣ),abs(pᵣ/pₗ)))+0.0))
    #gam = max( 0.5, gam2 ) 
	#pₗᵣ = 0.5*(pₗ+pᵣ) + KLR*rhohat*c̄*(p⁺+p⁻-1.0) - 
	#		0.5*(p⁺-p⁻)*(pᵣ-pₗ)

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


function SLAU2_HR(
    pₗ::Float64,uₗ::Float64,vₗ::Float64,wₗ::Float64,
	Tₗ::Float64,Y₁ₗ::Float64,ρₗ::Float64,Hₜₗ::Float64,cₗ::Float64,
    pᵣ::Float64,uᵣ::Float64,vᵣ::Float64,wᵣ::Float64,
	Tᵣ::Float64,Y₁ᵣ::Float64,ρᵣ::Float64,Hₜᵣ::Float64,cᵣ::Float64,
    Lco::Float64,Uco::Float64,Δt::Float64,
    w₁::Float64,w₂::Float64,cpi2::Float64,
    nx::Float64,ny::Float64,nz::Float64
)


	# properties of Left
	Uₙₗ = uₗ*nx + vₗ*ny + wₗ*nz
	# properties of Right
	Uₙᵣ = uᵣ*nx + vᵣ*ny + wᵣ*nz

	c̄ = 0.5*(cₗ + cᵣ)
	#M̄ₖ = √(ū^2 + v̄^2 + w̄^2)/c̄
	Mₗ = Uₙₗ/c̄
	Mᵣ = Uₙᵣ/c̄
	# calculate M+ and P+ for left state
	Mₗ⁺ = M_func(Mₗ,1.0,0.125)
	p⁺ = pre_func(Mₗ,1.0,0.0)
	# calculate M- and P- for left state
	Mᵣ⁻ = M_func(Mᵣ,-1.0,0.125)
	p⁻ = pre_func(Mᵣ,-1.0,0.0)
	# ======= carefully, sensitive ===========
	Mₗᵣ = Mₗ⁺ + Mᵣ⁻
	# =======================================
	M̄⁺ = 0.0
    M̄⁻ = 0.0
	if Mₗᵣ >= 0.0
        M̄⁺ = Mₗᵣ
        M̄⁻ = 0.0
	else
        M̄⁺ = 0.0
        M̄⁻ = Mₗᵣ
    end

	#KLR = sqrt(0.5*(uₗ^2+vₗ^2+wₗ^2+uᵣ^2+vᵣ^2+wᵣ^2))
	#g = -max(min(Mₗ,0.0),-1.0)*min(max(Mᵣ,0.0),1.0)
	#Mdash = min(1.0,KLR/c̄)
	#Vn = (ρₗ*abs(Uₙₗ)+ρᵣ*abs(Uₙᵣ)) / (ρₗ+ρᵣ)
	#Vnp = (1.0-g)*Vn + g*abs(Uₙₗ)
	#Vnm = (1.0-g)*Vn + g*abs(Uₙᵣ)
	#mdot = 0.5*(ρₗ*(Uₙₗ+Vnp)+ρᵣ*(Uₙᵣ-Vnm)-(1.0-Mdash)^2/c̄*(pᵣ-pₗ))

	ṁₗ = c̄ * ρₗ * M̄⁺
	ṁᵣ = c̄ * ρᵣ * M̄⁻
	#ṁₗ = 0.5*(mdot+abs(mdot))
	#ṁᵣ = 0.5*(mdot-abs(mdot))
	#Kᵤ = 0.5
	#pᵤ = -2.0 * Kᵤ * p⁺ * p⁻ * ρ̄  * c̄ * (Uₙᵣ-Uₙₗ)
	#pₗᵣ = p⁺*pₗ + p⁻*pᵣ +  pᵤ * ϕᵤ

	#pₗᵣ = 0.5*(pₗ+pᵣ) - 0.5*KLR/c̄*0.5*p⁺*p⁻*0.5*(pₗ+pᵣ)/c̄*(Uₙᵣ-Uₙₗ) + 
	#0.5*(KLR/c̄)*0.5*(pₗ+pᵣ)*(p⁺+p⁻-1.0) - 0.5*(p⁺-p⁻)*(pᵣ-pₗ)
	UU = 0.5*(Uₙₗ+Uₙᵣ)
	if UU > 0.0
		ṁₗ = ρₗ * UU
		ṁᵣ = 0.0
	else
		ṁₗ = 0.0
		ṁᵣ = ρᵣ * UU
	end

	#pₗᵣ = pₗ*p⁺ + pᵣ*p⁻

	pₗᵣ = 0.5*(pₗ + pᵣ)

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



function AUSMPWP_N(
    pₗ,uₗ,vₗ,wₗ,Tₗ,Y₁ₗ,ρₗ,Hₜₗ,cₗ,
    pᵣ,uᵣ,vᵣ,wᵣ,Tᵣ,Y₁ᵣ,ρᵣ,Hₜᵣ,cᵣ,
    Lco,Uco,Δt,
    w₁,w₂,cpi2,
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



function Roe(
    pₗ,uₗ,vₗ,wₗ,Tₗ,Y₁ₗ,ρₗ,Hₜₗ,cₗ,
    pᵣ,uᵣ,vᵣ,wᵣ,Tᵣ,Y₁ᵣ,ρᵣ,Hₜᵣ,cᵣ,
    Lco,Uco,Δt,
    w₁,w₂,cpi2,
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
    Ȳ = (Y₁ₗ+RT*Y₁ᵣ)/(1.0+RT)
	
	RT = 1.0
    ρ̂ = (ρₗ+RT*ρᵣ)/(1.0+RT)
    ū = (uₗ+RT*uᵣ)/(1.0+RT)
    v̄ = (vₗ+RT*vᵣ)/(1.0+RT)
    w̄ = (wₗ+RT*wᵣ)/(1.0+RT)
    Ĥ = (Hₜₗ+RT*Hₜᵣ)/(1.0+RT)
    Ȳ = (Y₁ₗ+RT*Y₁ᵣ)/(1.0+RT)
	Ûₙ = 0.5 * (Uₙₗ + Uₙᵣ)
	# Low-Diffusion Flux-Splitting Methods for Real Fluid Flows at All Speeds, 1999
	c̄ = 0.5*(cₗ + cᵣ)

	fluxAdvL = zeros(Float64,6)
	fluxAdvL[1] = ρₗ*Uₙₗ
	fluxAdvL[2] = ρₗ*Uₙₗ*uₗ + pₗ*nx
	fluxAdvL[3] = ρₗ*Uₙₗ*vₗ + pₗ*ny
	fluxAdvL[4] = ρₗ*Uₙₗ*wₗ + pₗ*nz
	fluxAdvL[5] = ρₗ*Uₙₗ*Hₜₗ
	fluxAdvL[6] = ρₗ*Uₙₗ*Y₁ₗ
	
	fluxAdvR = zeros(Float64,6)
	fluxAdvR[1] = ρᵣ*Uₙᵣ
	fluxAdvR[2] = ρᵣ*Uₙᵣ*uᵣ + pᵣ*nx
	fluxAdvR[3] = ρᵣ*Uₙᵣ*vᵣ + pᵣ*ny
	fluxAdvR[4] = ρᵣ*Uₙᵣ*wᵣ + pᵣ*nz
	fluxAdvR[5] = ρᵣ*Uₙᵣ*Hₜᵣ
	fluxAdvR[6] = ρᵣ*Uₙᵣ*Y₁ᵣ
	
	
	fluxDissDW = zeros(Float64,6)
	tmpUn = abs(Ûₙ)
	fluxDissDW[1] = tmpUn*(ρᵣ-ρₗ)
	fluxDissDW[2] = tmpUn*(ρᵣ*uᵣ-ρₗ*uₗ)
	fluxDissDW[3] = tmpUn*(ρᵣ*vᵣ-ρₗ*vₗ)
	fluxDissDW[4] = tmpUn*(ρᵣ*wᵣ-ρₗ*wₗ)
	fluxDissDW[5] = tmpUn*(ρᵣ*(Hₜᵣ-pᵣ/ρᵣ)-ρₗ*(Hₜₗ-pₗ/ρₗ))
	fluxDissDW[6] = tmpUn*(ρᵣ*Y₁ᵣ-ρₗ*Y₁ₗ)
	
	
	fluxDissDP = zeros(Float64,6)
	tmpUn = Ûₙ/c̄*(pᵣ-pₗ)
	fluxDissDP[1] = 0.0
	fluxDissDP[2] = tmpUn*nx
	fluxDissDP[3] = tmpUn*ny
	fluxDissDP[4] = tmpUn*nz
	fluxDissDP[5] = tmpUn*Ûₙ
	fluxDissDP[6] = 0.0
	
	fluxDissDU = zeros(Float64,6)
	tmpUn = (c̄-abs(Ûₙ))/c̄/c̄/ρ̂ *(pᵣ-pₗ) + Ûₙ/c̄*(Uₙᵣ-Uₙₗ) 
	fluxDissDU[1] = tmpUn*ρ̂ 
	fluxDissDU[2] = tmpUn*ρ̂ *ū
	fluxDissDU[3] = tmpUn*ρ̂ *v̄
	fluxDissDU[4] = tmpUn*ρ̂ *w̄
	fluxDissDU[5] = tmpUn*ρ̂ *Ĥ
	fluxDissDU[6] = tmpUn*ρ̂ *Ȳ
	

	# comp. convective flux
    flux = zeros(Float64,6)
	flux[:] = 0.5*(fluxAdvL[:]+fluxAdvR[:]) - 
				0.5*(fluxDissDW[:]+fluxDissDP[:]+fluxDissDU[:])


    return flux


end

