module solvers
export poisson_solver

using LinearAlgebra


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
        v[:, 1] .= 0
        v[:, end-1] .= 0 # because of forward staggering
        v[:, end] .= 0
    end

end

function y_momentum_GS_returning(v, P, dx, Nx, Ny, av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1, omega=1.5)
    V = zeros(Nx + 2, Ny + 2)
    for k in 1:1000
        for j in 2:Ny
            for i in 2:Nx+1
                test = (dx * (P[i, j] - P[i, j+1]) - av_ip1j[i, j] * v[i+1, j]
                        -
                        av_im1j[i, j] * v[i-1, j] - av_ijp1[i, j] * v[i, j+1]
                        -
                        av_ijm1[i, j] * v[i, j-1]) / av_ij[i, j]
                V[i, j] = (1 - omega) * v[i, j] + omega * test
            end
        end
        V[:, 1] .= 0
        V[:, end-1] .= 0 # because of forward staggering
        V[:, end] .= 0
        v[1, :] = -v[2, :]
    end
    return V
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

            au_ij[i, j] += 0.25 * dy * (u[i, j] + u[i+1, j]) + dy / (dx * Re)
            au_ij[i, j] += (-0.25) * dy * (u[i, j] + u[i-1, j]) + dy / (dx * Re)
            au_ij[i, j] += 0.25 * dx * (v[i, j] + v[i+1, j]) + dx / (dy * Re)
            au_ij[i, j] += (-0.25) * dx * (v[i, j-1] + v[i+1, j-1]) + dx / (dy * Re)

            if isnan(au_ij[i, j])
                throw(ErrorException("what? -joe biden"))
            end
        end
    end
    return au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1
end

function x_momentum_GS!(u, P, dy, Nx, Ny, au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1, Uref)
    """
    Performs Gauss-Seidel iteration for the x-momentum equation and MUTATES u

    """

    # au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1 = u_coefficients(u, v, dx, dy, Re, Nx, Ny)

    for k in 1:1000
        for j in 2:Ny+1
            for i in 2:Nx+1
                u[i, j] = (dy * (P[i, j] - P[i+1, j]) - au_ip1j[i, j] * u[i+1, j]
                           -
                           au_im1j[i, j] * u[i-1, j] - au_ijp1[i, j] * u[i, j+1]
                           -
                           au_ijm1[i, j] * u[i, j-1]) / au_ij[i, j]
                if isnan(u[i, j])
                    println("i=$i, j=$j")
                    throw(ErrorException("nan in u here"))
                end
            end
        end
        # Inlet
        u[1, :] .= Uref
        # Oulet
        u[end, :] = u[end-1, :]
        # No-slip
        u[:, 1] = -u[:, 2]
        u[:, end] = -u[:, end-1]
    end
end

function x_momentum_GS_returning(u, P, dy, Nx, Ny, au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1, u_in, omega=1.5)
    U = zeros(Nx + 2, Ny + 2)
    for k in 1:100
        for j in 2:Ny+1
            for i in 2:Nx+1
                test = (dy * (P[i, j] - P[i+1, j]) - au_ip1j[i, j] * u[i+1, j]
                        -
                        au_im1j[i, j] * u[i-1, j] - au_ijp1[i, j] * u[i, j+1]
                        -
                        au_ijm1[i, j] * u[i, j-1]) / au_ij[i, j]
                U[i, j] = (1 - omega) * u[i, j] + omega * test
                if isnan(test)
                    println("i=$i, j=$j")
                    throw(ErrorException("nan in u here"))
                end
            end
        end
        # Inlet
        U[1, :] = u_in
        # Oulet
        U[end, :] = U[end-1, :]
        # No-slip
        U[:, 1] = -U[:, 2]
        U[:, end] = -U[:, end-1]
    end
    return U
end

function pressure_correction(u, v, P, dx, dy, Nx, Ny, auij, avij, omega=1.5)
    """
    Solves the pressure correction equation and returns the pressure correction field
    """

    p_prime = zeros(Float64, Nx + 2, Ny + 2)
    scale = norm(P)

    for k in 1:1000
        for j in 3:Ny+1
            for i in 3:Nx+1
                a_e = -dy^2 / auij[i, j]
                a_w = -dy^2 / auij[i-1, j]
                a_n = -dx^2 / avij[i, j]
                a_s = -dx^2 / avij[i, j-1]
                a_p = a_e + a_w + a_n + a_s
                b = dy * (u[i-1, j] - u[i, j]) + dx * (v[i, j-1] - v[i, j])
                test = (a_e * p_prime[i+1, j] + a_w * p_prime[i-1, j] + a_n * p_prime[i, j+1] + a_s * p_prime[i, j-1] + b) / a_p
                p_prime[i, j] = (1 - omega) * p_prime[i, j] + omega * test
                if isnan(p_prime[i, j])
                    println("a_e=$a_e, a_w=$a_w, a_n=$a_n, a_s=$a_s, a_p=$a_p")
                    println("i=$i, j=$j")
                    # println("u[i-1, j]=$(u[i-1, j]), u[i, j]=$(u[i, j]), v[i, j-1]=$(v[i, j-1]), v[i, j]=$(v[i, j])")
                    println("b=$b")
                    println("p_prime[i+1, j]=$(p_prime[i+1, j]), p_prime[i-1, j]=$(p_prime[i-1, j]), p_prime[i, j+1]=$(p_prime[i, j+1]), p_prime[i, j-1]=$(p_prime[i, j-1])")
                    throw(ErrorException("NaN in pressure correction"))
                end
            end
        end
        p_prime[1, :] = -p_prime[2, :]
        if norm(p_prime) / scale < 1e-6
            return p_prime
        end
    end
    return p_prime
end
end
