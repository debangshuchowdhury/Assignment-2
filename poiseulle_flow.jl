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

function apply_BC!(u, v, P, P_ref, u_in, P_out, grad, dx)
    """
    Applies relevant boundary conditions to the variables
    """
    # Inlet
    # u[1, :] = u_in #(Re * PGrad) .* (0.5 .* y) .* (1 .- y)
    # u[2, :] = u[1, :]
    v[1, :] .= 0
    # P[1, :] .= P_ref
    P[2, :] = P[3, :] .- grad * dx
    P[1, :] = P[2, :] .- grad * dx
    # u[2, :] = u[3, :]
    u[1, :] = u[2, :]

    # Oulet
    # u[end, :] = u[end-1, :]
    # u[end, :] = u_in
    # u[end-1, :] = u[end, :]
    # u[end-1, :] = u[end-2, :]
    u[end, :] = u[end-1, :]
    P[end, :] .= 0
    v[end, :] .= 0



    # No-slip
    u[:, 1] = -u[:, 2]
    u[:, end] = -u[:, end-1]

    # No penetration
    v[:, 1] .= 0
    v[:, end-1] .= 0 # because of forward staggering
    v[:, end] .= 0
    P[:, end] = P[:, end-1]
    P[:, 2] = P[:, 3]
    P[:, 1] = P[:, 2]
end


function nancheck(a, k)
    if any(isnan, a)
        throw(ErrorException("NaN in this array. (Iteration $k)"))
    end
end

# Computational loop
begin
    alpha_p = 0.5
    alpha_u = 1
    alpha_v = 0.001
    Re = 10
    PGrad = -0.1
    u = ones(Nx + 2, Ny + 2)
    v = ones(Nx + 2, Ny + 2)
    P = zeros(Nx + 2, Ny + 2)
    # for i in range(Nx + 1, 1, step=-1)
    #     P[i, :] = P[i+1, :] .- PGrad * dx
    # end
    # Uref = 1.0
    U_inlet = (Re * -PGrad) .* (0.5 .* y) .* (1 .- y)
    U_inlet[1] = -U_inlet[2]
    U_inlet[end] = -U_inlet[end-1]
    # for l in 1:Nx+2
    #     u[l, :] = U_inlet
    # end
    P_ref = 1.0
    P_outlet = 0 #P_ref - PGrad * (Lx)
    ttt = plot(range(1, Nx + 2), P[:, 3])
    display(ttt)
    println("P outlet = $P_outlet")
    # apply_BC!(u, v, P, Uref, P_ref)
    apply_BC!(u, v, P, P_ref, U_inlet, P_outlet, PGrad, dx)
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
    # clf!()
    println("numx = $(Nx+2), numy = $(Ny+2)")
    println("MIN P = $(minimum(P))")
    println("MAX P = $(maximum(P))")
    # sleep(5)
    # p_prime = zeros(Nx + 2, Ny + 2)

    for k in 1:10000
        av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1 = solvers.v_coefficients(u, v, dx, dy, Re, Nx, Ny)
        au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1 = solvers.u_coefficients(u, v, dx, dy, Re, Nx, Ny)
        vstar = solvers.y_momentum_GS_returning(v, P, dx, Nx, Ny, av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1, alpha_v)
        # vstar = zeros(Nx + 2, Ny + 2)
        ustar = solvers.x_momentum_GS_returning(u, P, dy, Nx, Ny, au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1, U_inlet, alpha_u)
        # ustar = deepcopy(u)
        # solvers.x_momentum_upwind!(ustar, v, P, dx, dy, Nx, Ny, Re, U_inlet)
        p_prime = solvers.pressure_correction(ustar, vstar, P, dx, dy, Nx, Ny, au_ij, av_ij, 1.5)

        P += alpha_p * p_prime

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

        # nancheck(uc, k)
        # nancheck(vc, k)

        u = ustar + uc
        v = vstar + vc

        apply_BC!(u, v, P, P_ref, U_inlet, P_outlet, PGrad, dx)

        error = norm(uc) / norm(u)

        # clim = maximum(abs, P)
        if k % 10 == 0
            println("Iteration $k: Error=$error")
            println(", pprime max = ", maximum(abs.(p_prime)))
            println("ustar max = $(maximum(abs.(ustar)))")
            println("vstar max = $(maximum(abs.(vstar)))")
            println("--------------------------------\n")
            lal = contourf(x, y, P',
                xlabel="x",
                ylabel="y",
                title="Pressure (Iteration $k)",
                color=:viridis,
                levels=20,          # number of contour levels; adjust as needed
                linewidth=0,
                aspect_ratio=1,
                # clims=(-clim, clim)
            )
            display(lal)
            # sleep(1)
        end

        if maximum(P) > 100
            print("iteration $k: max P = $(maximum(P))")
            throw(ErrorException("P is diverging."))
        end
        if error <= 1e-4
            print("Converged in $k iterations. error = $error")

            break
        end
    end
end

begin
    tttt = plot(range(1, Nx + 2), P[:, 6])
    # plot!(ylims=[0.0122, 0.0130])
    # plot!(U_inlet, range(1, Ny + 2), linestyle=:dashdot)
    display(tttt)
end

begin
    asdf = plot(U_inlet, range(1, Ny + 2), label="real", color=:black)
    for i in range(1, Nx + 2)
        plot!(u[i, :], range(1, Ny + 2), label="$i")
        display(asdf)
        sleep(0.05)
        # if i > 5
        #     break
        # end
    end
end

contourf(x, y, v',
    xlabel="x",
    ylabel="y",
    # title="x velocity (Iteration $k)",
    color=:viridis,
    levels=20,          # number of contour levels; adjust as needed
    linewidth=0,
    aspect_ratio=1,
    # clims=(-clim, clim)
)