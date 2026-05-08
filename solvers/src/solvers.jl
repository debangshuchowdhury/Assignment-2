module solvers
export poisson_solver

using LinearAlgebra

function apply_BC!(P, dy, M)
    """
    Applying BC to the array P. ( Applying random stuff for now)
    Args:
        P: (Array) 2D array with ij indexing
    Returns:
        P: (Array) 2D array with BC applied
    """

    y = [(j - 0.5) * dy for j in 0:M+1]

    # Top Boundary
    P[:, end] = 2 .- P[:, end-1]
    # Bottom Boundary
    P[:, 1] = -P[:, 2]
    # Right Boundary
    P[end, :] = 2 .* y - P[end-1, :]
    # Left Boundary
    P[1, :] = -P[2, :]
    # return P

end

function poisson_solver(shape, f, dx, dy, tol)
    """
    Solves the poission equation on a 2D grid with N x M cells and source term f.
    Uses a slow algorithm (Gauss Siedel)

    Args:
        shape: (Tuple) dimensions of the grid (N, M)
        f: (Array) source term with ij indexing
        dx: (Float) grid spacing in x direction
        dy: (Float) grid spacing in y direction
        tol: (Float) tolerance for convergence

    Returns:
        P: (Array) solution to the Poisson equation with ij indexing
    """


    N, M = shape
    P = zeros(N + 2, M + 2)
    # P = apply_BC(P, dy, M)
    apply_BC!(P, dy, M)
    # F = zeros(N + 2, M + 2)
    # F[2:N+1, 2:M+1] = f
    # pold = deepcopy(P)

    for k in 1:100000
        # println("Iteration $k")
        res = 0
        for j in 2:M+1
            for i in 2:N+1
                Pnew = (-f[i-1, j-1] + (1 / dx^2) * (P[i+1, j] + P[i-1, j]) + (1 / dy^2) * (P[i, j+1] + P[i, j-1])) / (2 / dx^2 + 2 / dy^2)
                res += (P[i, j] - Pnew)^2
                P[i, j] = Pnew
            end
        end
        # P = apply_BC(P, dy, M)
        apply_BC!(P, dy, M)
        # display(P)
        if res / norm(P) <= tol
            println("Converged in $k iterations")
            break
        end
        # pold = deepcopy(P)
    end

    return P[2:N+1, 2:M+1]

end

end # module solvers
