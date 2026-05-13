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
    u[end, :] = u[end-1, :]


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
    Re = 100
    PGrad = 1
    u = ones(Nx + 2, Ny + 2)
    v = zeros(Nx + 2, Ny + 2)
    P = zeros(Nx + 2, Ny + 2)
    Uref = 1.0
    U_inlet = (Re * PGrad) .* (0.5 .* y) .* (1 .- y)
    P_ref = 1.0
    # apply_BC!(u, v, P, Uref, P_ref)
    apply_BC!(u, v, P, P_ref, U_inlet)

    hmm = contourf(x, y, P',
        xlabel="x",
        ylabel="y",
        title="Pressure",
        color=:viridis,
        levels=20,          # number of contour levels; adjust as needed
        linewidth=0,
        aspect_ratio=1)
    display(hmm)
    sleep(5)

    for k in 1:10000
        av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1 = v_coefficients(u, v, dx, dy, Re, Nx, Ny)
        au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1 = u_coefficients(u, v, dx, dy, Re, Nx, Ny)
        # y_momentum_GS!(v, P, dx, Nx, Ny, av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1)
        # x_momentum_GS!(u, P, dy, Nx, Ny, au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1, Uref)
        # p_prime = pressure_correction(u, v, P, dx, dy, Nx, Ny, au_ij, av_ij)
        vstar = y_momentum_GS_returning(v, P, dx, Nx, Ny, av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1)
        ustar = x_momentum_GS_returning(u, P, dy, Nx, Ny, au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1, Uref)
        p_prime = pressure_correction(ustar, vstar, P, dx, dy, Nx, Ny, au_ij, av_ij)

        print("pprime max = ", maximum(abs.(p_prime)))
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

        if k % 100 == 0
        end
        println("Iteration $k: Error=$error")
        println("--------------------------------\n")

        contourf!(x, y, P',
            xlabel="x",
            ylabel="y",
            title="Pressure",
            color=:viridis,
            levels=20,          # number of contour levels; adjust as needed
            linewidth=0,
            aspect_ratio=1)
        display(hmm)
        sleep(1)

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

# function v_coefficients(u, v, dx, dy, Re, Nx, Ny)
#     """
#     Recovers coefficients of the y-velocity field
#     Arguments:
#     - u: x-velocity field
#     - v: y-velocity field
#     - dx: grid spacing in x direction
#     - dy: grid spacing in y direction
#     - Re: Reynolds number
#     Returns:
#     - av_ij: coefficient of v[i, j]
#     - av_ip1j: coefficient of v[i+1, j]
#     - av_im1j: coefficient of v[i-1, j]
#     - av_ijp1: coefficient of v[i, j+1]
#     - av_ijm1: coefficient of v[i, j-1]
#     """
#     av_ij = zeros(Nx + 2, Ny + 2)
#     av_ip1j = zeros(Nx + 2, Ny + 2)
#     av_im1j = zeros(Nx + 2, Ny + 2)
#     av_ijp1 = zeros(Nx + 2, Ny + 2)
#     av_ijm1 = zeros(Nx + 2, Ny + 2)

#     for j in 2:Ny+1
#         for i in 2:Nx+1
#             av_ip1j[i, j] = +0.25 * dy * (u[i, j] + u[i, j+1]) - (dy / (dx * Re))
#             av_im1j[i, j] = -0.25 * dy * (u[i-1, j] + u[i-1, j+1]) - dy / (dx * Re)
#             av_ijp1[i, j] = +0.25 * dx * (v[i, j+1] + v[i, j]) - dx / (dy * Re)
#             av_ijm1[i, j] = -0.25 * dx * (v[i, j] + v[i, j-1]) - dx / (dy * Re)

#             av_ij[i, j] += -0.25 * dx * (v[i, j-1] + v[i, j])
#             av_ij[i, j] += dx / (dy * Re)
#             av_ij[i, j] += 0.25 * dy * (u[i, j] + u[i, j+1])
#             av_ij[i, j] += dy / (dx * Re)
#             av_ij[i, j] += 0.25 * dx * (v[i, j+1] + v[i, j])
#             av_ij[i, j] += dx / (dy * Re)
#             av_ij[i, j] += -0.25 * dy * (u[i-1, j] + u[i-1, j+1])
#             av_ij[i, j] += dy / (dx * Re)
#         end
#     end
#     return av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1
# end

# function y_momentum_GS!(v, P, dx, Nx, Ny, av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1)
#     """
#     Performs Gauss-Seidel iteration for the y-momentum equation and MUTATES v
#     """
#     # av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1 = v_coefficients(u, v, dx, dy, Re, Nx, Ny)

#     for k in 1:100
#         for j in 2:Ny
#             for i in 2:Nx+1
#                 v[i, j] = (dx * (P[i, j] - P[i, j+1]) - av_ip1j[i, j] * v[i+1, j]
#                            -
#                            av_im1j[i, j] * v[i-1, j] - av_ijp1[i, j] * v[i, j+1]
#                            -
#                            av_ijm1[i, j] * v[i, j-1]) / av_ij[i, j]
#             end
#         end
#         v[:, 1] .= 0
#         v[:, end-1] .= 0 # because of forward staggering
#         v[:, end] .= 0
#     end

# end

# function y_momentum_GS_returning(v, P, dx, Nx, Ny, av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1)
#     V = zeros(Nx + 2, Ny + 2)
#     for k in 1:100
#         for j in 2:Ny
#             for i in 2:Nx+1
#                 V[i, j] = (dx * (P[i, j] - P[i, j+1]) - av_ip1j[i, j] * v[i+1, j]
#                            -
#                            av_im1j[i, j] * v[i-1, j] - av_ijp1[i, j] * v[i, j+1]
#                            -
#                            av_ijm1[i, j] * v[i, j-1]) / av_ij[i, j]
#                 if V[i, j] == 0
#                 end
#             end
#         end
#         V[:, 1] .= 0
#         V[:, end-1] .= 0 # because of forward staggering
#         V[:, end] .= 0
#     end
#     return V
# end

# function u_coefficients(u, v, dx, dy, Re, Nx, Ny)
#     """
#     Recovers coefficients of the x-velocity field
#     Arguments:
#     - u: x-velocity field
#     - v: y-velocity field
#     - dx: grid spacing in x direction
#     - dy: grid spacing in y direction
#     - Re: Reynolds number
#     - Nx: number of grid points in x direction
#     - Ny: number of grid points in y direction
#     Returns:
#     - au_ij: coefficient of u[i, j]
#     - au_ip1j: coefficient of u[i+1, j]
#     - au_im1j: coefficient of u[i-1, j]
#     - au_ijp1: coefficient of u[i, j+1]
#     - au_ijm1: coefficient of u[i, j-1]
#     """
#     au_ij = zeros(Nx + 2, Ny + 2)
#     au_ip1j = zeros(Nx + 2, Ny + 2)
#     au_im1j = zeros(Nx + 2, Ny + 2)
#     au_ijp1 = zeros(Nx + 2, Ny + 2)
#     au_ijm1 = zeros(Nx + 2, Ny + 2)

#     for j in 2:Ny+1
#         for i in 2:Nx+1
#             au_ip1j[i, j] = +0.25 * dy * (u[i, j] + u[i+1, j]) - dy / (dx * Re)
#             au_im1j[i, j] = -0.25 * dy * (u[i, j] + u[i-1, j]) - dy / (dx * Re)
#             au_ijp1[i, j] = +0.25 * dx * (v[i, j] + v[i+1, j]) - dx / (dy * Re)
#             au_ijm1[i, j] = -0.25 * dx * (v[i, j-1] + v[i+1, j-1]) - dx / (dy * Re)

#             au_ij[i, j] += 0.25 * dy * (u[i, j] + u[i+1, j]) + dy / (dx * Re)
#             au_ij[i, j] += (-0.25) * dy * (u[i, j] + u[i-1, j]) + dy / (dx * Re)
#             au_ij[i, j] += 0.25 * dx * (v[i, j] + v[i+1, j]) + dx / (dy * Re)
#             au_ij[i, j] += (-0.25) * dx * (v[i, j-1] + v[i+1, j-1]) + dx / (dy * Re)

#             if isnan(au_ij[i, j])
#                 throw(ErrorException("what? -joe biden"))
#             end
#         end
#     end
#     return au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1
# end

# function x_momentum_GS!(u, P, dy, Nx, Ny, au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1, Uref)
#     """
#     Performs Gauss-Seidel iteration for the x-momentum equation and MUTATES u

#     """

#     # au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1 = u_coefficients(u, v, dx, dy, Re, Nx, Ny)

#     for k in 1:100
#         for j in 2:Ny+1
#             for i in 2:Nx+1
#                 u[i, j] = (dy * (P[i, j] - P[i+1, j]) - au_ip1j[i, j] * u[i+1, j]
#                            -
#                            au_im1j[i, j] * u[i-1, j] - au_ijp1[i, j] * u[i, j+1]
#                            -
#                            au_ijm1[i, j] * u[i, j-1]) / au_ij[i, j]
#                 if isnan(u[i, j])
#                     println("i=$i, j=$j")
#                     throw(ErrorException("nan in u here"))
#                 end
#             end
#         end
#         # Inlet
#         u[1, :] .= Uref
#         # Oulet
#         u[end, :] = u[end-1, :]
#         # No-slip
#         u[:, 1] = -u[:, 2]
#         u[:, end] = -u[:, end-1]
#     end
# end

# function x_momentum_GS_returning(u, P, dy, Nx, Ny, au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1, Uref)
#     U = zeros(Nx + 2, Ny + 2)
#     for k in 1:100
#         for j in 2:Ny+1
#             for i in 2:Nx+1
#                 U[i, j] = (dy * (P[i, j] - P[i+1, j]) - au_ip1j[i, j] * u[i+1, j]
#                            -
#                            au_im1j[i, j] * u[i-1, j] - au_ijp1[i, j] * u[i, j+1]
#                            -
#                            au_ijm1[i, j] * u[i, j-1]) / au_ij[i, j]
#                 if isnan(U[i, j])
#                     println("i=$i, j=$j")
#                     throw(ErrorException("nan in u here"))
#                 end
#             end
#         end
#         # Inlet
#         U[1, :] .= Uref
#         # Oulet
#         U[end, :] = U[end-1, :]
#         # No-slip
#         U[:, 1] = -U[:, 2]
#         U[:, end] = -U[:, end-1]
#     end
#     return U
# end

# function pressure_correction(u, v, P, dx, dy, Nx, Ny, auij, avij)
#     """
#     Solves the pressure correction equation and returns the pressure correction field
#     """

#     p_prime = zeros(Float64, Nx + 2, Ny + 2)
#     scale = norm(P)

#     for k in 1:1000
#         for j in 3:Ny+1
#             for i in 3:Nx+1
#                 a_e = -dy^2 / auij[i, j]
#                 a_w = -dy^2 / auij[i-1, j]
#                 a_n = -dx^2 / avij[i, j]
#                 a_s = -dx^2 / avij[i, j-1]
#                 a_p = a_e + a_w + a_n + a_s
#                 b = dy * (u[i-1, j] - u[i, j]) + dx * (v[i, j-1] - v[i, j])
#                 test = (a_e * p_prime[i+1, j] + a_w * p_prime[i-1, j] + a_n * p_prime[i, j+1] + a_s * p_prime[i, j-1] + b) / a_p
#                 p_prime[i, j] = 0.9999999 * p_prime[i, j] + 0.0000001 * test
#                 if isnan(p_prime[i, j])
#                     println("a_e=$a_e, a_w=$a_w, a_n=$a_n, a_s=$a_s, a_p=$a_p")
#                     println("i=$i, j=$j")
#                     # println("u[i-1, j]=$(u[i-1, j]), u[i, j]=$(u[i, j]), v[i, j-1]=$(v[i, j-1]), v[i, j]=$(v[i, j])")
#                     println("b=$b")
#                     println("p_prime[i+1, j]=$(p_prime[i+1, j]), p_prime[i-1, j]=$(p_prime[i-1, j]), p_prime[i, j+1]=$(p_prime[i, j+1]), p_prime[i, j-1]=$(p_prime[i, j-1])")
#                     throw(ErrorException("NaN in pressure correction"))
#                 end
#             end
#         end
#         p_prime[1, :] = -p_prime[2, :]
#         if norm(p_prime) / scale < 1e-6
#             return p_prime
#         end
#     end
#     return p_prime
# end
