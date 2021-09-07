export mesh




Nx = 5
Ny = 5

Lx = 1.0
Ly = 1.0
Δx  = Lx/Nx
Δy  = Ly/Ny

x_initial = 0.5*Δx
y_initial = 0.5*Δy
x = x_initial:Δx:Lx
y = y_initial:Δy:Ly

CFL = 0.5
Δt = CFL * Δx
