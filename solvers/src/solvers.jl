module solvers
# export poisson_solver

using LinearAlgebra


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

#     for j in 2:(Ny+1)
#         for i in 2:(Nx+1)
#             av_ip1j[i, j] = (+0.25 * dy * (u[i, j] + u[i, j+1]) - (dy / (dx * Re)))
#             av_im1j[i, j] = (-0.25 * dy * (u[i-1, j] + u[i-1, j+1]) - dy / (dx * Re))
#             av_ijp1[i, j] = (+0.25 * dx * (v[i, j+1] + v[i, j]) - dx / (dy * Re))
#             av_ijm1[i, j] = (-0.25 * dx * (v[i, j] + v[i, j-1]) - dx / (dy * Re))

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

# function y_momentum_GS_returning(v, P, dx, Nx, Ny, av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1, omega=1)
#     V = deepcopy(v)
#     for k in 1:2000
#         for j in 2:(Ny+1)
#             for i in 2:(Nx+1)
#                 test = (dx * (P[i, j] - P[i, j+1]) - av_ip1j[i, j] * V[i+1, j]
#                         -
#                         av_im1j[i, j] * V[i-1, j] - av_ijp1[i, j] * V[i, j+1]
#                         -
#                         av_ijm1[i, j] * V[i, j-1]) / av_ij[i, j]
#                 V[i, j] = (1 - omega) * V[i, j] + omega * test
#                 if isnan(test)
#                     println("i=$i, j=$j")
#                     throw(ErrorException("nan in v here"))
#                 end
#             end
#         end
#         # V[:, 1] .= 0
#         # V[:, end-1] .= 0 # because of forward staggering
#         # # V[:, end] .= 0
#         # v[1, :] = v[2, :]
#         # v[end, :] = v[end-1, :]
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

#     for j in 2:(Ny+1)
#         for i in 2:(Nx+1)
#             au_ip1j[i, j] = (+0.25 * dy * (u[i, j] + u[i+1, j]) - dy / (dx * Re))
#             au_im1j[i, j] = (-0.25 * dy * (u[i, j] + u[i-1, j]) - dy / (dx * Re))
#             au_ijp1[i, j] = (+0.25 * dx * (v[i, j] + v[i+1, j]) - dx / (dy * Re))
#             au_ijm1[i, j] = (-0.25 * dx * (v[i, j-1] + v[i+1, j-1]) - dx / (dy * Re))

#             au_ij[i, j] += 0.25 * dy * (u[i, j] + u[i+1, j]) + dy / (dx * Re)
#             au_ij[i, j] += (-0.25) * dy * (u[i, j] + u[i-1, j]) + dy / (dx * Re)
#             au_ij[i, j] += 0.25 * dx * (v[i, j] + v[i+1, j]) + dx / (dy * Re)
#             au_ij[i, j] += (-0.25) * dx * (v[i, j-1] + v[i+1, j-1]) + dx / (dy * Re)
#             # au_ij[i,j] = au_ip1j[i,j] + au_im1j[i,j] + au_ijp1[i,j] + au_ijm1[i,j] + 4*()

#             if isnan(au_ij[i, j])
#                 throw(ErrorException("what? -joe biden"))
#             end
#         end
#     end
#     return au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1
# end

# function x_momentum_GS_returning(u, P, dy, Nx, Ny, au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1, u_in, omega=1)
#     U = deepcopy(u)
#     for k in 1:2000
#         res = 0
#         for j in 2:(Ny+1)
#             for i in 2:(Nx+1)
#                 test = (dy * (P[i, j] - P[i+1, j]) - au_ip1j[i, j] * U[i+1, j]
#                         -
#                         au_im1j[i, j] * U[i-1, j] - au_ijp1[i, j] * U[i, j+1]
#                         -
#                         au_ijm1[i, j] * U[i, j-1]) / au_ij[i, j]
#                 pred = (1 - omega) * U[i, j] + omega * test
#                 res = maximum([abs((U[i, j] - pred))^2, res])
#                 U[i, j] = pred
#                 if isnan(test)
#                     println("i=$i, j=$j")
#                     throw(ErrorException("nan in u here"))
#                 end
#                 if pred > 100 * maximum(u_in)
#                     println("predi = $pred")
#                     throw(ErrorException("diverging"))
#                 end
#             end
#         end
#         # # Inlet
#         # U[1, :] = u_in
#         # # U[2, :] = U[3, :]
#         # # U[1, :] = U[2, :]

#         # # Oulet
#         # # U[end-1, :] = u_in    #U[end-1, :]
#         # # U[end-1, :] = U[end-2, :]
#         # U[end, :] = U[end-1, :]
#         # # U[end-2, :] = U[end-1, :]

#         # # No-slip
#         # U[:, 1] = -U[:, 2]
#         # U[:, end] = -U[:, end-1]

#         if res <= 1e-6
#             # println("xmom gs completed in $k iterations")
#             return U
#         end

#     end

#     return U
# end

# function pressure_correction(u, v, P, dx, dy, Nx, Ny, auij, avij, omega=1.5, p_prime=nothing)
#     """
#     Solves the pressure correction equation and returns the pressure correction field
#     """
#     if isnothing(p_prime)
#         p_prime = zeros(Float64, Nx + 2, Ny + 2)
#     end
#     # scale = norm(P)

#     for k in 1:10000
#         res = 0
#         for j in 3:(Ny+1)
#             for i in 3:(Nx+1)
#                 a_e = dy^2 / auij[i, j]
#                 a_w = dy^2 / auij[i-1, j]
#                 a_n = dx^2 / avij[i, j]
#                 a_s = dx^2 / avij[i, j-1]
#                 a_p = (a_e + a_w + a_n + a_s)
#                 b = dy * (u[i-1, j] - u[i, j]) + dx * (v[i, j-1] - v[i, j])
#                 test = (a_e * p_prime[i+1, j] + a_w * p_prime[i-1, j] + a_n * p_prime[i, j+1] + a_s * p_prime[i, j-1] + b) / a_p
#                 if isnan(p_prime[i, j])
#                     println("a_e=$a_e, a_w=$a_w, a_n=$a_n, a_s=$a_s, a_p=$a_p")
#                     println("i=$i, j=$j")
#                     # println("u[i-1, j]=$(u[i-1, j]), u[i, j]=$(u[i, j]), v[i, j-1]=$(v[i, j-1]), v[i, j]=$(v[i, j])")
#                     println("b=$b")
#                     println("p_prime[i+1, j]=$(p_prime[i+1, j]), p_prime[i-1, j]=$(p_prime[i-1, j]), p_prime[i, j+1]=$(p_prime[i, j+1]), p_prime[i, j-1]=$(p_prime[i, j-1])")
#                     throw(ErrorException("NaN in pressure correction"))
#                 end
#                 new = (1 - omega) * p_prime[i, j] + omega * test
#                 res += (p_prime[i, j] - res)^2
#                 p_prime[i, j] = new
#             end
#         end
#         # p_prime[1, :] .= 0
#         # p_prime[2, :] .= 0
#         # p_prime[end-1, :] .= 0
#         # p_prime[end, :] .= 0
#         if res < 1e-6
#             # println("breakout")
#             return p_prime
#         end
#     end
#     # println("pressure correction did not converge (err = $(norm(p_prime))")
#     return p_prime
# end

# function u_coefficients_upwind(u, v, dx, dy, Nx, Ny, Re)
#     """
#     Recovers coefficients of upwwind scheme for x momentum equation
#     """
#     au_ij = zeros(Nx + 2, Ny + 2)
#     au_ip1j = zeros(Nx + 2, Ny + 2)
#     au_im1j = zeros(Nx + 2, Ny + 2)
#     au_ijp1 = zeros(Nx + 2, Ny + 2)
#     au_ijm1 = zeros(Nx + 2, Ny + 2)

#     for j in 2:(Ny+1)
#         for i in 2:(Nx+1)
#             adv_e = 0.5 * (u[i+1, j] + u[i, j])
#             adv_w = 0.5 * (u[i-1, j] + u[i, j])
#             adv_n = 0.5 * (v[i+1, j] + v[i, j])
#             adv_s = 0.5 * (v[i, j-1] + v[i+1, j-1])
#             au_ip1j[i, j] = maximum([0, -adv_e]) * dy
#             au_im1j[i, j] = maximum([0, adv_w]) * dy
#             au_ijp1[i, j] = maximum([0, -adv_n]) * dx
#             au_ijm1[i, j] = maximum([0, adv_s]) * dx
#             au_ij[i, j] = au_ip1j[i, j] + au_im1j[i, j] + au_ijp1[i, j] + au_ijm1[i, j]
#         end
#     end
#     au_ip1j = au_ip1j .- (1 / Re) * (dy / dx)
#     au_im1j = au_im1j .- (1 / Re) * (dy / dx)
#     au_ijp1 = au_ijp1 .- (1 / Re) * (dx / dy)
#     au_ijm1 = au_ijm1 .- (1 / Re) * (dx / dy)
#     au_ij = au_ij .+ (2 / Re) * (dx / dy + dy / dx)

#     return au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1
# end

# function x_momentum_upwind(u, P, dy, Nx, Ny, au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1, omega=1, niter=2000)
#     """
#     upwind x momentum
#     """
#     U = deepcopy(u)
#     for k in 1:niter
#         # res = 0
#         for j in 2:(Ny+1)
#             for i in 2:(Nx+1)
#                 test = (dy * (P[i, j] - P[i+1, j]) - au_ip1j[i, j] * U[i+1, j]
#                         -
#                         au_im1j[i, j] * U[i-1, j] - au_ijp1[i, j] * U[i, j+1]
#                         -
#                         au_ijm1[i, j] * U[i, j-1]) / au_ij[i, j]
#                 U[i, j] = (1 - omega) * U[i, j] + omega * test
#             end
#         end
#     end
#     return U
# end

# function v_coefficients_upwind(u, v, dx, dy, Nx, Ny, Re)
#     """
#     Recovers coefficients of upwind scheme for y momentum equation
#     """
#     av_ij = zeros(Nx + 2, Ny + 2)
#     av_ip1j = zeros(Nx + 2, Ny + 2)
#     av_im1j = zeros(Nx + 2, Ny + 2)
#     av_ijp1 = zeros(Nx + 2, Ny + 2)
#     av_ijm1 = zeros(Nx + 2, Ny + 2)

#     for j in 2:(Ny+1)
#         for i in 2:(Nx+1)
#             adv_e = 0.5 * (u[i, j] + u[i, j+1])
#             adv_w = 0.5 * (u[i-1, j] + u[i-1, j+1])
#             adv_n = 0.5 * (v[i, j] + v[i, j+1])
#             adv_s = 0.5 * (v[i, j-1] + v[i, j])
#             av_ip1j[i, j] = maximum([0, -adv_e]) * dy
#             av_im1j[i, j] = maximum([0, adv_w]) * dy
#             av_ijp1[i, j] = maximum([0, -adv_n]) * dx
#             av_ijm1[i, j] = maximum([0, adv_s]) * dx
#             av_ij[i, j] = av_ip1j[i, j] + av_im1j[i, j] + av_ijp1[i, j] + av_ijm1[i, j]
#         end
#     end
#     av_ip1j = av_ip1j .- (1 / Re) * (dy / dx)
#     av_im1j = av_im1j .- (1 / Re) * (dy / dx)
#     av_ijp1 = av_ijp1 .- (1 / Re) * (dx / dy)
#     av_ijm1 = av_ijm1 .- (1 / Re) * (dx / dy)
#     av_ij = av_ij .+ (2 / Re) * (dx / dy + dy / dx)

#     return av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1
# end

# function y_momentum_upwind(v, P, dx, Nx, Ny, av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1, omega=1, niter=2000)
#     """
#     upwind y momentum
#     """
#     V = deepcopy(v)
#     for k in 1:niter
#         for j in 2:(Ny+1)
#             for i in 2:(Nx+1)
#                 test = (dx * (P[i, j] - P[i, j+1]) - av_ip1j[i, j] * V[i+1, j]
#                         -
#                         av_im1j[i, j] * V[i-1, j] - av_ijp1[i, j] * V[i, j+1]
#                         -
#                         av_ijm1[i, j] * V[i, j-1]) / av_ij[i, j]
#                 V[i, j] = (1 - omega) * V[i, j] + omega * test
#             end
#         end
#     end
#     return V
# end

# function x_momentum_GS_returning_step(u, P, dy, Nx, Ny, au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1, stepindx, stepindy, omega=1, niter=2000)
#     U = deepcopy(u)
#     for k in 1:niter
#         res = 0
#         for j in 2:(Ny+1)
#             for i in 2:(Nx+1)
#                 if i <= stepindx && j <= stepindy
#                     continue
#                 end
#                 test = (dy * (P[i, j] - P[i+1, j]) - au_ip1j[i, j] * U[i+1, j]
#                         -
#                         au_im1j[i, j] * U[i-1, j] - au_ijp1[i, j] * U[i, j+1]
#                         -
#                         au_ijm1[i, j] * U[i, j-1]) / au_ij[i, j]
#                 pred = (1 - omega) * U[i, j] + omega * test
#                 res = maximum([abs((U[i, j] - pred))^2, res])
#                 U[i, j] = pred
#                 if isnan(test)
#                     println("i=$i, j=$j")
#                     throw(ErrorException("nan in u here"))
#                 end
#             end
#         end

#         if res <= 1e-6
#             # println("xmom gs completed in $k iterations")
#             return U
#         end

#     end

#     return U
# end

# function y_momentum_GS_returning_step(v, P, dx, Nx, Ny, av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1, stepindx, stepindy, omega=1, niter=2000)
#     V = deepcopy(v)
#     for k in 1:niter
#         for j in 2:(Ny+1)
#             for i in 2:(Nx+1)
#                 if i <= stepindx && j <= stepindy
#                     continue
#                 end
#                 test = (dx * (P[i, j] - P[i, j+1]) - av_ip1j[i, j] * V[i+1, j]
#                         -
#                         av_im1j[i, j] * V[i-1, j] - av_ijp1[i, j] * V[i, j+1]
#                         -
#                         av_ijm1[i, j] * V[i, j-1]) / av_ij[i, j]
#                 V[i, j] = (1 - omega) * V[i, j] + omega * test
#                 if isnan(test)
#                     println("i=$i, j=$j")
#                     throw(ErrorException("nan in v here"))
#                 end
#             end
#         end
#     end
#     return V
# end

function u_coefficients_central_step(u, v, dx, dy, Nx, Ny, Re, stepindx, stepindy)
    """
    Recoveres coefficients of central fvm scheme for x momentum equation
    """
    au_ij = zeros(Nx + 2, Ny + 2)
    au_ip1j = zeros(Nx + 2, Ny + 2)
    au_im1j = zeros(Nx + 2, Ny + 2)
    au_ijp1 = zeros(Nx + 2, Ny + 2)
    au_ijm1 = zeros(Nx + 2, Ny + 2)

    for j in 2:(Ny+1)
        for i in 2:(Nx+1)
            if i<=stepindx && j<=stepindy
                continue
            end
            if i == stepindx+1 && j<=stepindy
                kw = 1
                au_im1j[i, j] = 0
            else
                kw = 1
                au_im1j[i, j] = 0.25 * dy * (u[i, j] + u[i-1, j])
            end
            if j == stepindy + 1 && i <= stepindx
                ks = 1
                au_ijm1[i, j] = 0
            else
                ks = 1
                au_ijm1[i, j] = 0.25 * dx * (v[i, j-1] + v[i+1, j-1])
            end
            au_ip1j[i, j] = -0.25 * dy * (u[i, j] + u[i+1, j])
            au_ijp1[i, j] = -0.25 * dx * (v[i, j] + v[i+1, j])

            au_ij[i, j] = - au_im1j[i, j] - au_ijm1[i, j] - au_ip1j[i, j] - au_ijp1[i, j]
            au_ij[i, j] += (1+kw)*(dy/dx)/Re + (1 + ks)*(dx/dy)/Re

            au_ip1j[i, j] += (1 / Re) * (dy / dx)
            au_im1j[i, j] += (1 / Re) * (dy / dx) * kw
            au_ijp1[i, j] += (1 / Re) * (dx / dy)
            au_ijm1[i, j] += (1 / Re) * (dx / dy) * ks
        end
    end

    return au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1
end

function v_coefficients_central_step(u, v, dx, dy, Nx, Ny, Re, stepindx, stepindy)
    """
    Recovers the coefficients of central fvm scheme for y momentum equation
    """
    av_ij = zeros(Nx + 2, Ny + 2)
    av_ip1j = zeros(Nx + 2, Ny + 2)
    av_im1j = zeros(Nx + 2, Ny + 2)
    av_ijp1 = zeros(Nx + 2, Ny + 2)
    av_ijm1 = zeros(Nx + 2, Ny + 2)

    for j in 2:(Ny+1)
        for i in 2:(Nx+1)
            if i <= stepindx && j <= stepindy
                continue
            end
            if i == stepindx + 1 && j <= stepindy
                kw = 1
                av_im1j[i, j] = 0
            else
                kw = 1
                av_im1j[i, j] = 0.25 * dy * (u[i-1, j] + u[i-1, j+1])
            end
            if j == stepindy + 1 && i <= stepindx
                ks = 1
                av_ijm1[i, j] = 0
            else
                ks = 1
                av_ijm1[i, j] = 0.25 * dx * (v[i, j] + v[i, j-1])
            end
            av_ip1j[i, j] = -0.25 * dy * (u[i, j] + u[i, j+1])
            av_ijp1[i, j] = -0.25 * dx * (v[i, j+1] + v[i, j])

            av_ij[i, j] = -av_ip1j[i, j] - av_im1j[i, j] - av_ijp1[i, j] - av_ijm1[i, j]
            av_ij[i, j] += (1 + ks) * (dx / dy) / Re + (1 + kw) * (dy / dx) / Re
            av_ip1j[i, j] += (1 / Re) * (dy / dx)
            av_im1j[i, j] += (1 / Re) * (dy / dx) * kw
            av_ijp1[i, j] += (1 / Re) * (dx / dy)
            av_ijm1[i, j] += (1 / Re) * (dx / dy) * ks
        end
    end

    return av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1
end

function u_coefficients_upwind_step(u, v, dx, dy, Nx, Ny, Re, stepindx, stepindy)
    """
    Recovers coefficients of upwwind scheme for x momentum equation
    """
    au_ij = zeros(Nx + 2, Ny + 2)
    au_ip1j = zeros(Nx + 2, Ny + 2)
    au_im1j = zeros(Nx + 2, Ny + 2)
    au_ijp1 = zeros(Nx + 2, Ny + 2)
    au_ijm1 = zeros(Nx + 2, Ny + 2)

    for j in 2:(Ny+1)
        for i in 2:(Nx+1)
            if i <= stepindx && j <= stepindy
                continue
            end
            if i == stepindx + 1 && j <= stepindy
                au_im1j[i, j] = 0
                kw = 1
            else
                adv_w = 0.5 * (u[i-1, j] + u[i, j])
                au_im1j[i, j] = maximum([0, adv_w]) * dy
                kw = 1
            end
            if j == stepindy + 1 && i <= stepindx
                au_ijm1[i, j] = 0
                ks = 1
            else
                adv_s = 0.5 * (v[i, j-1] + v[i+1, j-1])
                au_ijm1[i, j] = maximum([0, adv_s]) * dx
                ks = 1
            end
            adv_e = 0.5 * (u[i+1, j] + u[i, j])
            adv_n = 0.5 * (v[i+1, j] + v[i, j])
            au_ip1j[i, j] = maximum([0, -adv_e]) * dy
            au_ijp1[i, j] = maximum([0, -adv_n]) * dx
            au_ij[i, j] = au_ip1j[i, j] + au_im1j[i, j] + au_ijp1[i, j] + au_ijm1[i, j] + (1 + ks) * (dx / dy) / Re + (1 + kw) * (dy / dx) / Re
            au_ip1j[i, j] += (1 / Re) * (dy / dx)
            au_im1j[i, j] += (1 / Re) * (dy / dx) * kw
            au_ijp1[i, j] += (1 / Re) * (dx / dy)
            au_ijm1[i, j] += (1 / Re) * (dx / dy) * ks
        end
    end
    # au_ip1j = au_ip1j .- (1 / Re) * (dy / dx)
    # au_im1j = au_im1j .- (1 / Re) * (dy / dx)
    # au_ijp1 = au_ijp1 .- (1 / Re) * (dx / dy)
    # au_ijm1 = au_ijm1 .- (1 / Re) * (dx / dy)
    # au_ij = au_ij .+ (2 / Re) * (dx / dy + dy / dx)

    return au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1
end

function v_coefficients_upwind_step(u, v, dx, dy, Nx, Ny, Re, stepindx, stepindy)
    """
    Recovers coefficients of upwind scheme for y momentum equation
    """
    av_ij = zeros(Nx + 2, Ny + 2)
    av_ip1j = zeros(Nx + 2, Ny + 2)
    av_im1j = zeros(Nx + 2, Ny + 2)
    av_ijp1 = zeros(Nx + 2, Ny + 2)
    av_ijm1 = zeros(Nx + 2, Ny + 2)

    for j in 2:(Ny+1)
        for i in 2:(Nx+1)
            if i <= stepindx && j <= stepindy
                continue
            end
            if i == stepindx + 1 && j <= stepindy
                av_im1j[i, j] = 0
                kw = 1
            else
                adv_w = 0.5 * (u[i-1, j] + u[i-1, j+1])
                av_im1j[i, j] = maximum([0, adv_w]) * dy
                kw = 1
            end
            if j == stepindy + 1 && i <= stepindx
                av_ijm1[i, j] = 0
                ks = 1
            else
                adv_s = 0.5 * (v[i, j-1] + v[i, j])
                av_ijm1[i, j] = maximum([0, adv_s]) * dx
                ks = 1
            end
            adv_e = 0.5 * (u[i, j] + u[i, j+1])
            adv_n = 0.5 * (v[i, j] + v[i, j+1])
            av_ip1j[i, j] = maximum([0, -adv_e]) * dy
            av_ijp1[i, j] = maximum([0, -adv_n]) * dx

            av_ij[i, j] = av_ip1j[i, j] + av_im1j[i, j] + av_ijp1[i, j] + av_ijm1[i, j] + (1 + ks) * (dx / dy) / Re + (1 + kw) * (dy / dx) / Re
            av_ip1j[i, j] += (1 / Re) * (dy / dx)
            av_im1j[i, j] += (1 / Re) * (dy / dx) * kw
            av_ijp1[i, j] += (1 / Re) * (dx / dy)
            av_ijm1[i, j] += (1 / Re) * (dx / dy) * ks
        end
    end
    # av_ip1j = av_ip1j .- (1 / Re) * (dy / dx)
    # av_im1j = av_im1j .- (1 / Re) * (dy / dx)
    # av_ijp1 = av_ijp1 .- (1 / Re) * (dx / dy)
    # av_ijm1 = av_ijm1 .- (1 / Re) * (dx / dy)
    # av_ij = av_ij .+ (2 / Re) * (dx / dy + dy / dx)

    return av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1
end

function pressure_correction_step(u, v, dx, dy, Nx, Ny, auij, avij, stepindx, stepindy, omega=1.5, niter=10000)
    """
    Solves the pressure correction equation and returns the pressure correction field
    """

    p_prime = zeros(Float64, Nx + 2, Ny + 2)
    bs = 0
    for k in 1:niter
        res = 0
        bs = 0
        for j in 2:(Ny+1)
            for i in 2:(Nx+1)
                if i <= stepindx && j <= stepindy
                    continue
                end
                # a_e = auij[i, j] == 0.0 ? 0.0 : dy^2 / auij[i, j]
                a_e = dy^2 / auij[i, j]
                # a_w = dy^2 / auij[i-1, j]
                a_w = auij[i-1, j] == 0.0 ? 0.0 : dy^2 / auij[i-1, j]
                a_n = dx^2 / avij[i, j]
                # a_s = dx^2 / avij[i, j-1]
                a_s = avij[i, j-1] == 0.0 ? 0.0 : dx^2 / avij[i, j-1]
                a_p = (a_e + a_w + a_n + a_s)
                b = dy * (u[i-1, j] - u[i, j]) + dx * (v[i, j-1] - v[i, j])
                bs += b * b
                test = (a_e * p_prime[i+1, j] + a_w * p_prime[i-1, j] + a_n * p_prime[i, j+1] + a_s * p_prime[i, j-1] + b) / a_p
                if isnan(p_prime[i, j])
                    println("a_e=$a_e, a_w=$a_w, a_n=$a_n, a_s=$a_s, a_p=$a_p")
                    println("i=$i, j=$j")
                    # println("u[i-1, j]=$(u[i-1, j]), u[i, j]=$(u[i, j]), v[i, j-1]=$(v[i, j-1]), v[i, j]=$(v[i, j])")
                    println("b=$b")
                    println("p_prime[i+1, j]=$(p_prime[i+1, j]), p_prime[i-1, j]=$(p_prime[i-1, j]), p_prime[i, j+1]=$(p_prime[i, j+1]), p_prime[i, j-1]=$(p_prime[i, j-1])")
                    throw(ErrorException("NaN in pressure correction"))
                end
                new = (1 - omega) * p_prime[i, j] + omega * test
                res += (p_prime[i, j] - res)^2
                p_prime[i, j] = new
            end
        end
        p_prime[:, 1] = p_prime[:, 2]
        p_prime[1, :] = p_prime[2, :]
        p_prime[:, end] = p_prime[:, end-1]
        p_prime[1:stepindx, stepindy] = p_prime[1:stepindx, stepindy+1]
        p_prime[stepindx, 1:(stepindy-1)] = p_prime[stepindx+1, 1:(stepindy-1)]

        # p_prime[end, :] = p_prime[end-1, :]

    end
    # println("pressure correction did not converge (err = $(norm(p_prime))")
    return p_prime, bs
end

function x_momentum_upwind_step(u, P, dy, Nx, Ny, au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1, stepindx, stepindy, omega=1, niter=2000)
    U = deepcopy(u)
    for k in 1:niter
        # res = 0
        for j in 2:(Ny+1)
            for i in 2:(Nx+1)
                if j <= stepindy && i <= stepindx
                    continue
                end
                test = -(dy * (-P[i, j] + P[i+1, j]) - au_ip1j[i, j] * U[i+1, j]
                         -
                         au_im1j[i, j] * U[i-1, j] - au_ijp1[i, j] * U[i, j+1]
                         -
                         au_ijm1[i, j] * U[i, j-1]) / au_ij[i, j]
                U[i, j] = (1 - omega) * U[i, j] + omega * test
            end
        end
        U[end, :] = U[end-1, :]
        U[:, 1] = -U[:, 2]
        U[:, end] = -U[:, end-1]
        U[1:(stepindx-1), stepindy] = -U[1:(stepindx-1), stepindy+1]
        U[stepindx, 1:(stepindy-1)] .= 0

    end
    return U
end

function y_momentum_upwind_step(v, P, dx, Nx, Ny, av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1, stepindx, stepindy, omega=1, niter=2000)
    V = deepcopy(v)
    for k in 1:niter
        for j in 2:Ny
            for i in 2:(Nx+1)
                if i <= stepindx && j <= stepindy
                    continue
                end
                # if i == Nx+2
                #     test = V[i-1,j] + (V[i])
                # end
                test = -(dx * (-P[i, j] + P[i, j+1]) - av_ip1j[i, j] * V[i+1, j]
                         -
                         av_im1j[i, j] * V[i-1, j] - av_ijp1[i, j] * V[i, j+1]
                         -
                         av_ijm1[i, j] * V[i, j-1]) / av_ij[i, j]
                V[i, j] = (1 - omega) * V[i, j] + omega * test
            end
        end
        V[end, :] = V[end-1, :]
        V[1, :] .= 0
        V[:, end-1] .= 0
        V[:, end] .= 0
        V[1:(stepindx-1), stepindy] .= 0
        V[stepindx, 1:stepindy] = -V[stepindx+1, 1:stepindy]
    end
    return V
end


end
