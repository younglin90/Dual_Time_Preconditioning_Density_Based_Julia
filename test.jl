using Plots




function eosIdeal(
    γ, cᵥ,
    p, u, v, w, T
)
	u² = u^2+v^2+w^2
	#= 
	cᵥ = 1401.0 #707.0
    γ = 1.330 #1.4 =#
	
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
	#= 
    p∞ = 621780000.0
	cᵥ = 3610.0
    γ = 1.19
    b = 6.7212e-4
	q = -1177788.0 =#
    
    p∞ = 7.028e8
	cᵥ = 3610.0
    γ = 1.19
    b = 0.00067212
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


function eosNASG2(
    p, u, v, w, T
)
	u² = u^2+v^2+w^2
	#= 
    p∞ = 621780000.0
	cᵥ = 3610.0
    γ = 1.19
    b = 6.7212e-4
	q = -1177788.0 =#
    
    p∞ = 0.0
	cᵥ = 955.0
    γ = 1.47
    b = 0.0
	q = 2077616.0
	
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



Xinp = 0.0:0.001:1.0


X = zeros(Float64, length(Xinp))
Y = zeros(Float64, length(Xinp))

for i in 1:length(X)
    ρᵢ = zeros(Float64, 2, 1)
    ρᵢ = zeros(Float64, 2, 1)
    Hₜᵢ = zeros(Float64, 2, 1)
    cᵢ = zeros(Float64, 2, 1)
    ∂ρ∂pᵢ = zeros(Float64, 2, 1)
    ∂ρ∂Tᵢ = zeros(Float64, 2, 1)
    ∂Hₜ∂pᵢ = zeros(Float64, 2, 1)
    ∂Hₜ∂Tᵢ = zeros(Float64, 2, 1)
#= 
    Y₁ = Xinp[i] #0.0

    Y₁ = max(min(Y₁,1.0),0.0)
    Y₂ = 1.0 - Y₁
 =#
    

    α₁ = Xinp[i]
    α₁ = max(min(α₁,1.0),0.0)
    α₂ = 1.0 - α₁


    j = 1
    (ρᵢ[j], Hₜᵢ[j], cᵢ[j], ∂ρ∂pᵢ[j], ∂ρ∂Tᵢ[j], ∂Hₜ∂pᵢ[j], ∂Hₜ∂Tᵢ[j]) = 
    #eosIdeal(1.4, 707.0, 101325.0, 0.0, 0.0, 0.0, 300.0)
    eosNASG(101325.0, 0.0, 0.0, 0.0, 300.0)
    j = 2
    (ρᵢ[j], Hₜᵢ[j], cᵢ[j], ∂ρ∂pᵢ[j], ∂ρ∂Tᵢ[j], ∂Hₜ∂pᵢ[j], ∂Hₜ∂Tᵢ[j]) = 
    #eosIdeal(1.33, 1401.0, 101325.0, 0.0, 0.0, 0.0, 300.0)
    eosIdeal(1.4, 707.0, 101325.0, 0.0, 0.0, 0.0, 300.0)
    #eosNASG2(101325.0, 0.0, 0.0, 0.0, 300.0)

#= 
    ρ = 1.0/(Y₁/ρᵢ[1]+Y₂/ρᵢ[2])
    α₁ = ρ*Y₁/ρᵢ[1]
    α₁ = max(min(α₁,1.0),0.0)
    α₂ = 1.0 - α₁
 =#

    ρ = α₁*ρᵢ[1]+α₂*ρᵢ[2]
    Y₁ = ρᵢ[1]*α₁/ρ
    Y₂ = 1.0 - Y₁

    Hₜ = Y₁*Hₜᵢ[1] + Y₂*Hₜᵢ[2]

   # ∂ρ∂p = α₁*∂ρ∂pᵢ[1] + α₂*∂ρ∂pᵢ[2]
    #∂ρ∂T = α₁*∂ρ∂Tᵢ[1] + α₂*∂ρ∂Tᵢ[2]
    
    ∂ρ∂p = ρ^2*(Y₁/ρᵢ[1]^2*∂ρ∂pᵢ[1] + Y₂/ρᵢ[2]^2*∂ρ∂pᵢ[2])
    ∂ρ∂T = ρ^2*(Y₁/ρᵢ[1]^2*∂ρ∂Tᵢ[1] + Y₂/ρᵢ[2]^2*∂ρ∂Tᵢ[2])

    ∂Hₜ∂p = Y₁*∂Hₜ∂pᵢ[1] + Y₂*∂Hₜ∂pᵢ[2]
    ∂Hₜ∂T = Y₁*∂Hₜ∂Tᵢ[1] + Y₂*∂Hₜ∂Tᵢ[2]

    ∂ρ∂Y₁ = -ρ^2.0*(1.0/ρᵢ[1]-1.0/ρᵢ[2])
    ∂Hₜ∂Y₁ = Hₜᵢ[1]-Hₜᵢ[2]

    c = √(ρ*∂Hₜ∂T/(ρ*∂ρ∂p*∂Hₜ∂T+(1.0-ρ*∂Hₜ∂p)*∂ρ∂T))


    X[i] = α₁
    Y[i] = c


    #invGam = α₁/(1.19-1) + (1.0-α₁)/(1.4-1.0)
    #gamInvGam = α₁*1.19*7.028e8/(1.19-1.0)+(1.0-α₁)*1.4*0.0/(1.4-1.0)
    #c2 = √(101325.0*(invGam+1.0)+gamInvGam)/(invGam*ρ)

    #c2 = √(1.0/(1.0/1500.0^2*(α₁^2+α₁*(1.0-α₁)*1000.0/1.0)+1.0/483.0^2*((1.0-α₁)^2+α₁*(1.0-α₁)*1.0/1000.0)))

    c2 = √(1.0/(ρ*(α₁/ρᵢ[1]/1500.0^2+α₂/ρᵢ[2]/487.0^2)))

    #Y[i] = c2

    #println(ρ," ",Y₁," ",Y₂," ",α₁," ",α₂," ",c," ",c2)

end

plot(X,Y,yaxis=:log)

#plot(X,Y,yticks = [0 10 100 1000],yaxis=:log)
#= 
println(Y₁)
println(Y₂)
println(α₁)
println(α₂)
println(ρ)
println(∂ρ∂p)
println(∂ρ∂T)
println(∂Hₜ∂T)
println(∂Hₜ∂p)
println(c)
 =#

