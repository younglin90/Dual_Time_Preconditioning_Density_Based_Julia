#= include("./module.jl")
include("./controls.jl")
using .mesh =#

function EOS!(
    👉::controls, cell::Vector{mesh.Cell}
)


    for i in cell
        ρᵢ = zeros(Float64, 2, 1)
        ρᵢ = zeros(Float64, 2, 1)
        Hₜᵢ = zeros(Float64, 2, 1)
        cᵢ = zeros(Float64, 2, 1)
        ∂ρ∂pᵢ = zeros(Float64, 2, 1)
        ∂ρ∂Tᵢ = zeros(Float64, 2, 1)
        ∂Hₜ∂pᵢ = zeros(Float64, 2, 1)
        ∂Hₜ∂Tᵢ = zeros(Float64, 2, 1)

        i.var[👉.Y₁] = max(min(i.var[👉.Y₁],1.0),0.0)
        i.var[👉.Y₂] = 1.0 - i.var[👉.Y₁]

        j = 1
        (ρᵢ[j], Hₜᵢ[j], cᵢ[j], ∂ρ∂pᵢ[j], ∂ρ∂Tᵢ[j], ∂Hₜ∂pᵢ[j], ∂Hₜ∂Tᵢ[j]) = 
        eosNASG(i.var[👉.p], i.var[👉.u], i.var[👉.v], i.var[👉.w], i.var[👉.T])
        j = 2
        (ρᵢ[j], Hₜᵢ[j], cᵢ[j], ∂ρ∂pᵢ[j], ∂ρ∂Tᵢ[j], ∂Hₜ∂pᵢ[j], ∂Hₜ∂Tᵢ[j]) = 
        eosIdeal(i.var[👉.p], i.var[👉.u], i.var[👉.v], i.var[👉.w], i.var[👉.T])

        ρ = 1.0/(i.var[👉.Y₁]/ρᵢ[1]+i.var[👉.Y₂]/ρᵢ[2])

        i.var[👉.α₁] = ρ*i.var[👉.Y₁]/ρᵢ[1]
        i.var[👉.α₁] = max(min(i.var[👉.α₁],1.0),0.0)
        i.var[👉.α₂] = 1.0 - i.var[👉.α₁]

        i.var[👉.ρ] = ρ
        i.var[👉.Hₜ] = i.var[👉.Y₁]*Hₜᵢ[1] + i.var[👉.Y₂]*Hₜᵢ[2]

        #i.var[👉.∂ρ∂p] = i.var[👉.α₁]*∂ρ∂pᵢ[1] + i.var[👉.α₂]*∂ρ∂pᵢ[2]
        #i.var[👉.∂ρ∂T] = i.var[👉.α₁]*∂ρ∂Tᵢ[1] + i.var[👉.α₂]*∂ρ∂Tᵢ[2]
        i.var[👉.∂ρ∂p] = 
			ρ^2*(i.var[👉.Y₁]/ρᵢ[1]^2*∂ρ∂pᵢ[1] + 
				 i.var[👉.Y₂]/ρᵢ[2]^2*∂ρ∂pᵢ[2])
        i.var[👉.∂ρ∂T] = 
			ρ^2*(i.var[👉.Y₁]/ρᵢ[1]^2*∂ρ∂Tᵢ[1] + 
				 i.var[👉.Y₂]/ρᵢ[2]^2*∂ρ∂Tᵢ[2])

        i.var[👉.∂Hₜ∂p] = i.var[👉.Y₁]*∂Hₜ∂pᵢ[1] + i.var[👉.Y₂]*∂Hₜ∂pᵢ[2]
        i.var[👉.∂Hₜ∂T] = i.var[👉.Y₁]*∂Hₜ∂Tᵢ[1] + i.var[👉.Y₂]*∂Hₜ∂Tᵢ[2]

        i.var[👉.∂ρ∂Y₁] = -ρ^2.0*(1.0/ρᵢ[1]-1.0/ρᵢ[2])
        i.var[👉.∂Hₜ∂Y₁] = Hₜᵢ[1]-Hₜᵢ[2]

        i.var[👉.c] = i.var[👉.∂ρ∂p] + 1.0/ρ*i.var[👉.∂ρ∂T]/i.var[👉.∂Hₜ∂T]*(1.0-ρ*i.var[👉.∂Hₜ∂p])
        i.var[👉.c] = √(1.0/i.var[👉.c])


    end


end



function faceEOS!(
	p::Float64,u::Float64,v::Float64,w::Float64,T::Float64,Y₁::Float64
)


	ρᵢ = zeros(Float64, 2, 1)
	ρᵢ = zeros(Float64, 2, 1)
	Hₜᵢ = zeros(Float64, 2, 1)
	cᵢ = zeros(Float64, 2, 1)
	∂ρ∂pᵢ = zeros(Float64, 2, 1)
	∂ρ∂Tᵢ = zeros(Float64, 2, 1)
	∂Hₜ∂pᵢ = zeros(Float64, 2, 1)
	∂Hₜ∂Tᵢ = zeros(Float64, 2, 1)

	Y₁ = max(min(Y₁,1.0),0.0)
	Y₂ = 1.0 - Y₁

	j = 1
	(ρᵢ[j], Hₜᵢ[j], cᵢ[j], ∂ρ∂pᵢ[j], ∂ρ∂Tᵢ[j], ∂Hₜ∂pᵢ[j], ∂Hₜ∂Tᵢ[j]) = 
	eosNASG(p, u, v, w, T)
	j = 2
	(ρᵢ[j], Hₜᵢ[j], cᵢ[j], ∂ρ∂pᵢ[j], ∂ρ∂Tᵢ[j], ∂Hₜ∂pᵢ[j], ∂Hₜ∂Tᵢ[j]) = 
	eosIdeal(p, u, v, w, T)

	ρ = 1.0/(Y₁/ρᵢ[1]+Y₂/ρᵢ[2])

	α₁ = ρ*Y₁/ρᵢ[1]
	α₁ = max(min(α₁,1.0),0.0)
	α₂ = 1.0 - α₁

	Hₜ = Y₁*Hₜᵢ[1] + Y₂*Hₜᵢ[2]

	#∂ρ∂p = α₁*∂ρ∂pᵢ[1] + α₂*∂ρ∂pᵢ[2]
	#∂ρ∂T = α₁*∂ρ∂Tᵢ[1] + α₂*∂ρ∂Tᵢ[2]
    ∂ρ∂p = ρ^2*(Y₁/ρᵢ[1]^2*∂ρ∂pᵢ[1] + Y₂/ρᵢ[2]^2*∂ρ∂pᵢ[2])
    ∂ρ∂T = ρ^2*(Y₁/ρᵢ[1]^2*∂ρ∂Tᵢ[1] + Y₂/ρᵢ[2]^2*∂ρ∂Tᵢ[2])

	∂Hₜ∂p = Y₁*∂Hₜ∂Tᵢ[1] + Y₂*∂Hₜ∂Tᵢ[2]
	∂Hₜ∂T = Y₁*∂Hₜ∂Tᵢ[1] + Y₂*∂Hₜ∂Tᵢ[2]

	#∂ρ∂Y₁ = -ρ*ρ*(1.0/ρᵢ[1]-1.0/ρᵢ[2])
	#∂Hₜ∂Y₁ = Hₜᵢ[1]-Hₜᵢ[2]

	c = ∂ρ∂p + 1.0/ρ*∂ρ∂T/∂Hₜ∂T*(1.0-ρ*∂Hₜ∂p)
	c = √(1.0/c)

	return ρ, Hₜ, c


end




function eosIdeal(
    p, u, v, w, T
)
	u² = u^2+v^2+w^2
	
	cᵥ = 707.0
    γ = 1.4
	
	cₚ = γ*cᵥ
		
	# density of each phase
	ρ = 1.0/( (γ-1.0)*cᵥ*T/p )
	c = √( γ/ρ*p )

	# d(rho)/d(p)
	∂ρ∂p = ρ/p
	
	# d(rho)/d(T)
	∂ρ∂T = -ρ/T

	# d(h)/d(p)
	∂Hₜ∂p = 0.0
	# d(h)/d(T)
	∂Hₜ∂T = γ*cᵥ

	# internal energy of each phase
	# internal_energy = (p+gam*pinf)/(gam-1.0)*(1.0/rhoi-bNASG)+q

	# enthalpy of each phase
	h = γ*cᵥ*T

		   
	# eti = internal_energy + 0.5*usqrt
	Hₜ = h + 0.5*u²

	# cvi = cv
	# cpi = cp	

    return ρ, Hₜ, c, ∂ρ∂p, ∂ρ∂T, ∂Hₜ∂p, ∂Hₜ∂T
end		

function eosNASG(
    p, u, v, w, T
)
	u² = u^2+v^2+w^2
	
    p∞ = 621780000.0
	cᵥ = 3610.0
    γ = 1.19
    b = 6.7212e-4
	q = -1177788.0
	
	cₚ = γ*cᵥ
		
	# density of each phase
	ρ = 1.0/( (γ-1.0)*cᵥ*T/(p+p∞)+b )
	c = √( γ/(ρ*ρ)*(p+p∞)/(1.0/ρ-b) )

	# d(rho)/d(p)
	∂ρ∂p = ρ*ρ*(1.0/ρ-b)/(p+p∞)
	
	# d(rho)/d(T)
	∂ρ∂T = -ρ*ρ*(1.0/ρ-b)/T

	# d(h)/d(p)
	∂Hₜ∂p = b
	# d(h)/d(T)
	∂Hₜ∂T = γ*cᵥ

	# internal energy of each phase
	# internal_energy = (p+gam*pinf)/(gam-1.0)*(1.0/rhoi-bNASG)+q

	# enthalpy of each phase
	h = γ*cᵥ*T + b*p + q

		   
	# eti = internal_energy + 0.5*usqrt
	Hₜ = h + 0.5*u²

	# cvi = cv
	# cpi = cp	

    return ρ, Hₜ, c, ∂ρ∂p, ∂ρ∂T, ∂Hₜ∂p, ∂Hₜ∂T
end