"""
Construct a uniform flow at Re=200
"""

# Set the flow parameters
Re = 200
U∞ = (1.0, 0.0)

# Set the domain grid and time step size
nx = 152; ny = 102
Ly = 2.0
Δx = Ly/(ny-2)
Δt = min(0.5*Δx,0.5*Δx^2*Re)
w₀ = Nodes(Dual,(nx,ny))
xg, yg = coordinates(w₀,dx=Δx)
t = 0.0
tf = 2Δt

# Set up initial conditions
w₀ .= 0.0
fx, fy = Float64[], Float64[]
thist, uhist = [], []
