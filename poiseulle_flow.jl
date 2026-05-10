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
    u[1, :] = U_ref
    P[1, :] = P_ref

    # Oulet
    u[end, :] = u[end-1, :]


    # No-slip
    u[:, 1] = -u[:, 2]
    u[:, end] = -u[:, end-1]

    # No penetration
    v[:, 1] = 0
    v[:, end-1] = 0 # because of forward staggering
    v[:, end] = 0
    P[:, end] = P[:, end-1]
    P[:, 1] = P[:, 2]
end


function v_coefficients(u, v, dx, dy, Re, Nx, Ny)
    """
    Recovers coefficients of the y-velocity field
    Arguments:
    - u: x-velocity field
    - v: y-velocity field
    - dx: grid spacing in x direction
    - dy: grid spacing in y direction
    - Re: Reynolds number
    Returns:
    - av_ij: coefficient of v[i, j]
    - av_ip1j: coefficient of v[i+1, j]
    - av_im1j: coefficient of v[i-1, j]
    - av_ijp1: coefficient of v[i, j+1]
    - av_ijm1: coefficient of v[i, j-1]
    """
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
            av_ij[i, j] += dy / (dx * Re)
        end
    end
    return av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1
end

function y_momentum_GS!(v, P, dx, Nx, Ny, av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1)
    """
    Performs Gauss-Seidel iteration for the y-momentum equation and MUTATES v
    """
    # av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1 = v_coefficients(u, v, dx, dy, Re, Nx, Ny)

    for k in 1:100
        for j in 2:Ny
            for i in 2:Nx+1
                v[i, j] = (dx * (P[i, j] - P[i, j+1]) - av_ip1j[i, j] * v[i+1, j]
                           -
                           av_im1j[i, j] * v[i-1, j] - av_ijp1[i, j] * v[i, j+1]
                           -
                           av_ijm1[i, j] * v[i, j-1]) / av_ij[i, j]
            end
        end
        v[:, 1] = 0
        v[:, end-1] = 0 # because of forward staggering
        v[:, end] = 0
    end

end

function u_coefficients(u, v, dx, dy, Re, Nx, Ny)
    """
    Recovers coefficients of the x-velocity field
    Arguments:
    - u: x-velocity field
    - v: y-velocity field
    - dx: grid spacing in x direction
    - dy: grid spacing in y direction
    - Re: Reynolds number
    - Nx: number of grid points in x direction
    - Ny: number of grid points in y direction
    Returns:
    - au_ij: coefficient of u[i, j]
    - au_ip1j: coefficient of u[i+1, j]
    - au_im1j: coefficient of u[i-1, j]
    - au_ijp1: coefficient of u[i, j+1]
    - au_ijm1: coefficient of u[i, j-1]
    """
    au_ij = zeros(Nx + 2, Ny + 2)
    au_ip1j = zeros(Nx + 2, Ny + 2)
    au_im1j = zeros(Nx + 2, Ny + 2)
    au_ijp1 = zeros(Nx + 2, Ny + 2)
    au_ijm1 = zeros(Nx + 2, Ny + 2)

    for j in 2:Ny+1
        for i in 2:Nx+1
            au_ip1j[i, j] = +0.25 * dy * (u[i, j] + u[i+1, j]) - dy / (dx * Re)
            au_im1j[i, j] = -0.25 * dy * (u[i, j] + u[i-1, j]) - dy / (dx * Re)
            au_ijp1[i, j] = +0.25 * dx * (v[i, j] + v[i+1, j]) - dx / (dy * Re)
            au_ijm1[i, j] = -0.25 * dx * (v[i, j-1] + v[i+1, j-1]) - dx / (dy * Re)

            au_ij[i, j] += +0.25 * dy * (u[i, j] + u[i+1, j]) + dy / (dx * Re)
            au_ij[i, j] += -0.25 * dy * (u[i, j] + u[i-1, j]) + dy / (dx * Re)
            au_ij[i, j] += +0.25 * dx * (v[i, j] + v[i+1, j]) + dx / (dy * Re)
            au_ij[i, j] += -0.25 * dx * (v[i, j-1] + v[i+1, j-1]) + dx / (dy * Re)
        end
    end
    return au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1
end

function x_momentum_GS!(u, P, dy, Nx, Ny, au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1)
    """
    Performs Gauss-Seidel iteration for the x-momentum equation and MUTATES u

    """

    # au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1 = u_coefficients(u, v, dx, dy, Re, Nx, Ny)

    for k in 1:100
        for j in 2:Ny+1
            for i in 2:Nx+1
                u[i, j] = (dy * (P[i, j] - P[i+1, j]) - au_ip1j[i, j] * u[i+1, j]
                           -
                           au_im1j[i, j] * u[i-1, j] - au_ijp1[i, j] * u[i, j+1]
                           -
                           au_ijm1[i, j] * u[i, j-1]) / au_ij[i, j]
            end
        end
        # Inlet
        u[0, :] = U_ref
        # Oulet
        u[end, :] = u[end-1, :]
        # No-slip
        u[:, 0] = -u[:, 1]
        u[:, end] = -u[:, end-1]
    end
end

function pressure_correction(u, v, P, dx, dy, Nx, Ny, auij, avij)
    """
    Solves the pressure correction equation and returns the pressure correction field
    """

    # au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1 = u_coefficients(u, v, dx, dy, Re, Nx, Ny)
    # av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1 = v_coefficients(u, v, dx, dy, Re, Nx, Ny)
    p_prime = zeros(Nx + 2, Ny + 2)
    scale = norm(P)

    for k in 1:10000
        for j in 2:Ny+1
            for i in 2:Nx+1
                a_e = dy^2 / auij[i, j]
                a_w = dy^2 / auij[i-1, j]
                a_n = dx^2 / avij[i, j]
                a_s = dx^2 / avij[i, j-1]
                a_p = a_e + a_w + a_n + a_s
                b = dy * (u[i-1, j] - u[i, j]) + dx * (v[i, j-1] - v[i, j])
                p_prime[i, j] = (a_e * p_prime[i+1, j] + a_w * p_prime[i-1, j] + a_n * p_prime[i, j+1] + a_s * p_prime[i, j-1] + b) / a_p
            end
        end
        p_prime[1, :] = -p_prime[2, :]
        if norm(p_prime) / scale < 1e-6
            return p_prime
        end
    end
    return p_prime
end