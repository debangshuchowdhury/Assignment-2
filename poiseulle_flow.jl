"""
Implementation of 2D poiseulle flow using SIMPLE. Using forward staggering for both u and v.
Assumptions:
- Steady state
- Incompressible
- No body forces
- Fully Developed
- Planar 
Boundary Conditions:
- No slip at the walls (u = 0)
- Pressure is specified only at the inlet
- Velocity specified at the inlet (u = uref)
- Neumann BC at the outlet (du/dx = 0)
- No penetration at the walls (v = 0)
"""


using LinearAlgebra, Plots, Statistics
# pythonplot()


# Define computational domain
begin
    Ly = 1.0
    Lx = 5.0
    Ny = 10
    Nx = convert(Int32, Ny * (Lx / Ly))
    dx = Lx / Nx
    dy = Ly / Ny
    x = collect(round.(LinRange(-dx / 2, Lx + dx / 2, Nx + 2), digits=6))
    y = collect(round.(LinRange(-dy / 2, Ly + dy / 2, Ny + 2), digits=6))
    println(size(x))
    println(size(y))
end

# Plotting computational domain to check
begin
    p = plot([0, 0], [0, Ly], color=:black, label=:false)
    plot!([0, Lx], [0, 0], color=:black, label=:false)
    plot!([Lx, Lx], [0, Ly], color=:black, label=:false)
    plot!([0, Lx], [Ly, Ly], color=:black, label=:false)
    for yc in y
        for xc in x
            if xc > 0 && yc > 0 && xc < Lx && yc < Ly
                scatter!([xc], [yc], color=:green, markersize=:1, label=false)
            else
                scatter!([xc], [yc], color=:red, markersize=:1, label=false)
            end
        end
    end
    # scatter!(x, y, markersize=:1, label=:false, color=:black)
    plot!(aspect_ratio=1)
    plot!(ylims=[-0.5, Ly + 0.5], xlims=[-0.5, Lx + 0.5])
    display(p)
end


function apply_BC!(u, v, P, U_ref, P_ref)
    """
    Applies relevant boundary conditions to the variables
    """
    # Inlet
    u[0, :] = U_ref
    P[0, :] = P_ref

    # Oulet
    u[end, :] = u[end-1, :]


    # No-slip
    u[:, 0] = -u[:, 1]
    u[:, end] = -u[:, end-1]

    # No penetration
    v[:, 0] = 0
    v[:, end-1] = 0 # because of forward staggering
    v[:, end] = 0
    P[:, end] = P[:, end-1]
    P[:, 0] = P[:, 1]
end


"""
   Performs Gauss-Seidel iteration for the y-momentum equation and MUTATES v
   Arguments:
   - u: x-velocity field
   - v: y-velocity field
   - P: pressure field
   - dx: grid spacing in x direction
   - dy: grid spacing in y direction
   - Re: Reynolds number
   """
function y_momentum_GS!(u, v, P, dx, dy, Re, Nx, Ny)
    av_ij = zeros(Nx + 2, Ny + 2)
    av_ip1j = zeros(Nx + 2, Ny + 2)
    av_im1j = zeros(Nx + 2, Ny + 2)
    av_ijp1 = zeros(Nx + 2, Ny + 2)
    av_ijm1 = zeros(Nx + 2, Ny + 2)

    for j in 2:Ny+1
        for i in 2:Nx+1
            av_ip1j[i, j] = +0.25 * dy * (u[i, j] + u[i, j+1]) - (dy / (dx * Re))
            av_im1j[i, j] = -0.25 * dy * (u[i-1, j] + u[i-1, j+1]) - dy / (dx * Re)
            av_ijp1[i, j] = +0.25 * dx * (v[i, j+1] + v[i, j]) - dx / (dy * Re)
            av_ijm1[i, j] = -0.25 * dx * (v[i, j] + v[i, j-1]) - dx / (dy * Re)

            av_ij[i, j] += -0.25 * dx * (v[i, j-1] + v[i, j])
            av_ij[i, j] += dx / (dy * Re)
            av_ij[i, j] += 0.25 * dy * (u[i, j] + u[i, j+1])
            av_ij[i, j] += dy / (dx * Re)
            av_ij[i, j] += 0.25 * dx * (v[i, j+1] + v[i, j])
            av_ij[i, j] += dx / (dy * Re)
            av_ij[i, j] += -0.25 * dy * (u[i-1, j] + u[i-1, j+1])
        end
    end

    for k in 1:100
        for j in 2:Ny+1
            for i in 2:Nx+1
                v[i, j] = (dx * (P[i, j] - P[i, j+1]) - av_ip1j[i, j] * v[i+1, j]
                           -
                           av_im1j[i, j] * v[i-1, j] - av_ijp1[i, j] * v[i, j+1]
                           -
                           av_ijm1[i, j] * v[i, j-1]) / av_ij[i, j]
            end
        end
    end

end