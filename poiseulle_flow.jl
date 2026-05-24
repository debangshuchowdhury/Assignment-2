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


using LinearAlgebra, Plots, Statistics, Pkg
Pkg.activate("solvers/")
using solvers
# pythonplot()


# Define computational domain
begin
    Ly = 1.0
    Lx = 4.0
    Ny = 5
    Nx = convert(Int32, Ny * (Lx / Ly))
    dx = Lx / Nx
    dy = Ly / Ny
    x = collect(round.(LinRange(-dx / 2, Lx + dx / 2, Nx + 2), digits=6))
    y = collect(round.(LinRange(-dy / 2, Ly + dy / 2, Ny + 2), digits=6))
    println(size(x))
    println(size(y))
end

# Plotting computational domain to check
# begin
#     p = plot([0, 0], [0, Ly], color=:black, label=:false)
#     plot!([0, Lx], [0, 0], color=:black, label=:false)
#     plot!([Lx, Lx], [0, Ly], color=:black, label=:false)
#     plot!([0, Lx], [Ly, Ly], color=:black, label=:false)
#     for yc in y
#         for xc in x
#             if xc > 0 && yc > 0 && xc < Lx && yc < Ly
#                 scatter!([xc], [yc], color=:green, markersize=:1, label=false)
#             else
#                 scatter!([xc], [yc], color=:red, markersize=:1, label=false)
#             end
#         end
#     end
#     # scatter!(x, y, markersize=:1, label=:false, color=:black)
#     plot!(aspect_ratio=1)
#     plot!(ylims=[-0.5, Ly + 0.5], xlims=[-0.5, Lx + 0.5])
#     display(p)
# end

function apply_BC!(u, v, P, P_ref, u_in)
    """
    Applies relevant boundary conditions to the variables
    """
    # Inlet
    u[1, :] = u_in #(Re * PGrad) .* (0.5 .* y) .* (1 .- y)
    P[1, :] .= P_ref

    # Oulet
    # u[end, :] = u[end-1, :]
    u[end, :] = u_in



    # No-slip
    u[:, 1] = -u[:, 2]
    u[:, end] = -u[:, end-1]

    # No penetration
    v[:, 1] .= 0
    v[:, end-1] .= 0 # because of forward staggering
    v[:, end] .= 0
    # P[:, end] = P[:, end-1]
    # P[:, 1] = P[:, 2]
end



function nancheck(a, k)
    if any(isnan, a)
        throw(ErrorException("NaN in this array. (Iteration $k)"))
    end
end

# Computational loop
begin
    alpha = 1
    Re = 1
    PGrad = 0.1
    u = ones(Nx + 2, Ny + 2)
    v = zeros(Nx + 2, Ny + 2)
    P = zeros(Nx + 2, Ny + 2)
    Uref = 1.0
    U_inlet = (Re * PGrad) .* (0.5 .* y) .* (1 .- y)
    P_ref = 1.0
    # apply_BC!(u, v, P, Uref, P_ref)
    apply_BC!(u, v, P, P_ref, U_inlet)
    println("u max = $(maximum(abs.(u)))")

    # clim = maximum(abs, P)
    hmm = contourf(x, y, P',
        xlabel="x",
        ylabel="y",
        title="Pressure (initial condition)",
        color=:viridis,
        levels=20,          # number of contour levels; adjust as needed
        linewidth=0,
        aspect_ratio=1,
        # clims=(-clim, clim)
    )
    display(hmm)
    # sleep(5)

    for k in 1:10000
        av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1 = solvers.v_coefficients(u, v, dx, dy, Re, Nx, Ny)
        au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1 = solvers.u_coefficients(u, v, dx, dy, Re, Nx, Ny)
        # y_momentum_GS!(v, P, dx, Nx, Ny, av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1)
        # x_momentum_GS!(u, P, dy, Nx, Ny, au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1, Uref)
        # p_prime = pressure_correction(u, v, P, dx, dy, Nx, Ny, au_ij, av_ij)
        vstar = solvers.y_momentum_GS_returning(v, P, dx, Nx, Ny, av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1)
        ustar = solvers.x_momentum_GS_returning(u, P, dy, Nx, Ny, au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1, U_inlet)
        p_prime = solvers.pressure_correction(ustar, vstar, P, dx, dy, Nx, Ny, au_ij, av_ij)

        P += alpha * p_prime

        nancheck(ustar, k)
        nancheck(vstar, k)
        nancheck(p_prime, k)
        nancheck(P, k)

        uc = zeros(Nx + 2, Ny + 2)
        vc = zeros(Nx + 2, Ny + 2)
        for j in 2:Ny+1
            for i in 2:Nx+1
                uc[i, j] = (dy / au_ij[i, j]) * (p_prime[i, j] - p_prime[i+1, j])
                vc[i, j] = (dx / av_ij[i, j]) * (p_prime[i, j] - p_prime[i, j+1])
            end
        end

        nancheck(uc, k)
        nancheck(vc, k)

        u = ustar + alpha * uc
        v = vstar + alpha * vc

        apply_BC!(u, v, P, P_ref, U_inlet)

        error = maximum(abs.(uc)) / Uref
        if error <= 1e-10
            print("Converged in $k iterations. error = $error")

            break
        end


        println("Iteration $k: Error=$error")
        println(", pprime max = ", maximum(abs.(p_prime)))
        println("ustar max = $(maximum(abs.(ustar)))")
        println("--------------------------------\n")

        # clim = maximum(abs, P)
        contourf!(x, y, P',
            xlabel="x",
            ylabel="y",
            title="Pressure (Iteration $k)",
            color=:viridis,
            levels=20,          # number of contour levels; adjust as needed
            linewidth=0,
            aspect_ratio=1,
            # clims=(-clim, clim)
        )
        display(hmm)
        # sleep(1)

        if maximum(P) > 100
            print("iteration $k")
            throw(ErrorException("P is diverging."))
        end

    end
end

begin
    contourf(x, y, P',
        xlabel="x",
        ylabel="y",
        title="Pressure",
        color=:viridis,
        levels=20,          # number of contour levels; adjust as needed
        linewidth=0,
        aspect_ratio=1)

end
