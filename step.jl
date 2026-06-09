"""
Implementation of 2D backwards facing step using SIMPLE. Using forward staggering for both u and v.
"""


using LinearAlgebra, Plots, Statistics, Pkg
Pkg.activate("solvers/")
using solvers


# Define grid
begin
    Ly = 1
    Lx = 10.0
    Ny = 10
    Nx = convert(Int32, Ny * (Lx / Ly))
    dx = Lx / Nx
    dy = Ly / Ny
    x = collect(round.(LinRange(-dx / 2, Lx + dx / 2, Nx + 2), digits=6))
    y = collect(round.(LinRange(-dy / 2, Ly + dy / 2, Ny + 2), digits=6))
    println("Nx = $Nx, Ny = $Ny")

    # step dimensions
    step_height = floor(Ny / 2) * dy
    step_length = floor(Nx / 3) * dx
    stepindx = convert(Int32, floor(step_length / dx)) + 1
    stepindy = convert(Int32, floor(step_height / dy)) + 1
    mask = zeros(Bool, Nx + 2, Ny + 2)
    for i in 1:Nx+2
        for j in 1:Ny+2
            if x[i] >= 0 && x[i] <= step_length && y[j] >= 0 && y[j] <= step_height
                mask[i, j] = true
            end
        end
    end
    println("stpindx = $stepindx, stepindy = $stepindy")
end


# Plotting computational domain to check
begin
    p = Plots.plot([0, 0], [0, Ly], color=:black, label=:false)
    Plots.plot!([0, Lx], [0, 0], color=:black, label=:false)
    Plots.plot!([Lx, Lx], [0, Ly], color=:black, label=:false)
    Plots.plot!([0, Lx], [Ly, Ly], color=:black, label=:false)
    for yc in y
        for xc in x
            if xc > 0 && yc > 0 && xc < Lx && yc < Ly && !(xc <= step_length && yc <= step_height)
                Plots.scatter!([xc], [yc], color=:green, markersize=:4, label=false)
            elseif xc >= 0 && yc >= 0 && xc <= step_length && yc <= step_height
                Plots.scatter!([xc], [yc], color=:blue, markersize=:4, label=false)
            else
                Plots.scatter!([xc], [yc], color=:red, markersize=:4, label=false)
            end

        end
    end
    # scatter!(x, y, markersize=:1, label=:false, color=:black)
    Plots.scatter!([x[stepindx]], [y[stepindy]], color=:black, markersize=:6, label="step")
    Plots.plot!(aspect_ratio=1)
    Plots.plot!(ylims=[-0.5, Ly + 0.5], xlims=[-0.5, Lx + 0.5])
    Plots.hline!([step_height], color=:black, label=:false)
    Plots.vline!([step_length], color=:black, label=:false)
    display(p)
end


function apply_BC!(u, v, P, P_ref, u_ref, dx, dy, stepindx, stepindy, Ny)
    """
    Applies relevant boundary conditions to the variables
    """

    # Inlet
    # u[1, stepindy:end-1] = LinRange(0, uref, Ny - stepindy + 2)
    u[1, :] .= u_ref
    # P[2, :] = P[3, :]
    P[1, :] = P[2, :]
    # P[1, :] .= 1
    # P[2, :] .= 1
    v[1, :] .= 0


    # outlet
    u[end, :] = u[end-1, :] # .+ dx / dy * (v[end-1, :] - v[end-1, :])
    v[end, 1:end-1] = v[end-1, 1:end-1] #- (u[])
    P[end, :] .= P_ref
    # P[end, :] = P[end-1, :]


    # walls
    P[:, end] = P[:, end-1]
    P[:, 1] = P[:, 2]
    u[:, 1] = -u[:, 2]
    # u[:, end] = -u[:, end-1]
    u[:, end] .= 0 #u_ref
    v[:, 1] .= 0
    v[:, end-1] .= 0
    v[:, end] .= 0
    # v[:, end] = v[:, end-1]


    # step 
    u[1:stepindx, stepindy] = -u[1:stepindx, stepindy+1]
    u[1:stepindx-1, 1:stepindy-1] .= 0
    v[1:stepindx, 1:stepindy] .= 0
    P[1:stepindx, stepindy] = P[1:stepindx, stepindy+1]
    P[stepindx, 1:stepindy] = P[stepindx+1, 1:stepindy]
    P[1:stepindx-1, 1:stepindy-1] .= 1
    u[1:stepindx, 1:stepindy] .= 0
    v[stepindx, 1:stepindy-1] = -v[stepindx+1, 1:stepindy-1]
end


begin
    alpha_p = 1
    alpha_u = 1.5
    alpha_v = 1.5
    Re = 10
    uref = 1
    Pref = 0
    u = ones(Nx + 2, Ny + 2)
    v = ones(Nx + 2, Ny + 2)
    P = zeros(Nx + 2, Ny + 2)
    apply_BC!(u, v, P, Pref, uref, dx, dy, stepindx, stepindy, Ny)

    kkk = Plots.contourf(x, y, P',
        xlabel="x",
        ylabel="y",
        title="P (initial)",
        color=:viridis,
        levels=20,          # number of contour levels; adjust as needed
        linewidth=0,
        aspect_ratio=1,
        # clims=(-clim, clim)
        # xlims=[7.5, 15],
        # ylims=[-dy, 2]
    )
    # Plots.hline!([step_height], color=:black, label=:false)
    # Plots.vline!([step_length], color=:black, label=:false)
    Plots.plot!([0, step_length], [step_height, step_height], color=:black, label=:false)
    Plots.plot!([step_length, step_length], [0, step_height], color=:black, label=:false)
    display(kkk)
    # sakdw


    for k in 1:1000000
        # av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1 = solvers.v_coefficients(u, v, dx, dy, Re, Nx, Ny)
        # au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1 = solvers.u_coefficients(u, v, dx, dy, Re, Nx, Ny)
        av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1 = solvers.v_coefficients_upwind(u, v, dx, dy, Nx, Ny, Re)
        au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1 = solvers.u_coefficients_upwind(u, v, dx, dy, Nx, Ny, Re)
        # vstar = solvers.y_momentum_GS_returning_step(v, P, dx, Nx, Ny, av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1, stepindx, stepindy, alpha_v, 2000)
        # ustar = solvers.x_momentum_GS_returning_step(u, P, dy, Nx, Ny, au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1, stepindx, stepindy, alpha_u, 2000)
        vstar = solvers.y_momentum_upwind_step(v, P, dx, Nx, Ny, av_ij, av_ip1j, av_im1j, av_ijp1, av_ijm1, stepindx, stepindy, alpha_v, 100)
        ustar = solvers.x_momentum_upwind_step(u, P, dy, Nx, Ny, au_ij, au_ip1j, au_im1j, au_ijp1, au_ijm1, stepindx, stepindy, alpha_u, 100)
        # solvers.x_momentum_upwind!(ustar, v, P, dx, dy, Nx, Ny, Re, U_inlet)
        p_prime, b = solvers.pressure_correction_step(ustar, vstar, dx, dy, Nx, Ny, au_ij, av_ij, stepindx, stepindy, alpha_p, 200)

        P = P + alpha_p * p_prime

        uc = zeros(Nx + 2, Ny + 2)
        vc = zeros(Nx + 2, Ny + 2)
        for j in 2:Ny+1
            for i in 2:Nx+1
                if i <= stepindx && j <= stepindy
                    continue
                end
                uc[i, j] = (dy / au_ij[i, j]) * (p_prime[i, j] - p_prime[i+1, j])
                vc[i, j] = (dx / av_ij[i, j]) * (p_prime[i, j] - p_prime[i, j+1])
            end
        end

        u = ustar + uc
        v = vstar + uc

        apply_BC!(u, v, P, Pref, uref, dx, dy, stepindx, stepindy, Ny)


        error = max(norm(uc) / norm(u), norm(vc) / norm(v))

        if k % 100 == 0
            println("Iteration $k: Error=$error, b=$b")
            println(", pprime max = ", maximum(abs.(p_prime)))
            println("ustar max = $(maximum(abs.(ustar)))")
            println("vstar max = $(maximum(abs.(vstar)))")
            println("--------------------------------\n")
            lal = Plots.contourf(x, y, P',
                xlabel="x",
                ylabel="y",
                title="P (Iteration $k)",
                color=:viridis,
                levels=20,          # number of contour levels; adjust as needed
                linewidth=0,
                aspect_ratio=1,
                # clims=(-clim, clim)
                # xlims=[7.5, 15],
                # ylims=[-dy, 2]
            )
            # Plots.hline!([step_height], color=:black, label=:false)
            # Plots.vline!([step_length], color=:black, label=:false)
            Plots.plot!([0, step_length], [step_height, step_height], color=:black, label=:false)
            Plots.plot!([step_length, step_length], [0, step_height], color=:black, label=:false)
            Plots.display(lal)
            # sleep(1)
            # ighg
        end

        if maximum(P) > 1000
            print("iteration $k: max P = $(maximum(P))")
            throw(ErrorException("P is diverging."))
        end
        if b <= 1e-15 #error <= 1e-10 || 
            print("Converged in $k iterations. error = $error")

            break
        end

    end
end

Plots.contourf(x, y, P',
    xlabel="x",
    ylabel="y",
    title="Final P Field",
    color=:viridis,
    levels=1000,
    linewidth=0,
    aspect_ratio=1,
    # colormap=:fire
    # clims=(-clim, clim)
)


Plots.contourf(x, y, sqrt.(u .* u + v .* v)',
    xlabel="x",
    ylabel="y",
    # title="x velocity (Iteration $k)",
    color=:viridis,
    levels=200,
    linewidth=0,
    aspect_ratio=1,
    # clims=(-clim, clim)
)

Plots.plot(u[30, :], y)



# AI GENERATED 

using GLMakie
# using LinearInterpolations  # for bilinear interpolation

function bilinear_interp(xc, yc, f, xi, yi)
    # Find surrounding indices
    i = searchsortedfirst(xc, xi) - 1
    j = searchsortedfirst(yc, yi) - 1

    # Clamp to valid range
    i = clamp(i, 1, length(xc) - 1)
    j = clamp(j, 1, length(yc) - 1)

    # Interpolation weights
    tx = (xi - xc[i]) / (xc[i+1] - xc[i])
    ty = (yi - yc[j]) / (yc[j+1] - yc[j])

    return (1 - tx) * (1 - ty) * f[i, j] + tx * (1 - ty) * f[i+1, j] +
           (1 - tx) * ty * f[i, j+1] + tx * ty * f[i+1, j+1]
end

function plot_streamlines(u, v, x, y, stepindx, stepindy, dx, dy)
    Nx2 = size(u, 1)
    Ny2 = size(u, 2)

    # Interpolate to cell centers (interior only: i=2:Nx+1, j=2:Ny+1)
    u_c = 0.5 .* (u[1:end-1, :] .+ u[2:end, :])   # (Nx+1, Ny+2)
    v_c = 0.5 .* (v[:, 1:end-1] .+ v[:, 2:end])   # (Nx+2, Ny+1)

    u_c = u_c[2:end, 2:end-1]    # (Nx, Ny)
    v_c = v_c[2:end-1, 2:end]    # (Nx, Ny)

    # Interior cell center coordinates
    xc = x[2:end-1]
    yc = y[2:end-1]

    # Mask step
    u_c[1:stepindx-1, 1:stepindy-1] .= 0.0
    v_c[1:stepindx-1, 1:stepindy-1] .= 0.0

    # Velocity function for streamplot
    function velocity_field(xi, yi)
        xi_c = clamp(xi, xc[1], xc[end])
        yi_c = clamp(yi, yc[1], yc[end])
        if xi_c <= xc[stepindx-1] && yi_c <= yc[stepindy-1]
            return Point2f(0.0, 0.0)
        end
        ui = bilinear_interp(xc, yc, u_c, xi_c, yi_c)
        vi = bilinear_interp(xc, yc, v_c, xi_c, yi_c)
        return Point2f(ui, vi)
    end

    fig = Figure(size=(1200, 500))
    ax = Axis(fig[1, 1],
        xlabel="x",
        ylabel="y",
        title="Streamlines — Backwards Facing Step",
        aspect=DataAspect()
    )

    streamplot!(ax, velocity_field,
        xc[1] .. xc[end],
        yc[1] .. yc[end],
        colormap=:viridis,
        gridsize=(80, 40),
        arrow_size=8
    )

    # Draw step as solid rectangle
    poly!(ax,
        Rect(xc[1], yc[1], xc[stepindx-1] - xc[1], yc[stepindy-1] - yc[1]),
        color=:gray30
    )

    display(fig)
    return fig
end


fig = plot_streamlines(u, v, x, y, stepindx, stepindy, dx, dy)

save("streamlines1.png", fig)


function plot_velocity_arrows(u, v, x, y, stepindx, stepindy)

    # Interpolate to cell centers (interior only)
    u_c = 0.5 .* (u[1:end-1, :] .+ u[2:end, :])
    v_c = 0.5 .* (v[:, 1:end-1] .+ v[:, 2:end])

    u_c = u_c[2:end, 2:end-1]    # (Nx, Ny)
    v_c = v_c[2:end-1, 2:end]    # (Nx, Ny)

    # Interior cell center coordinates
    xc = x[2:end-1]
    yc = y[2:end-1]

    # Mask step
    u_c[1:stepindx-1, 1:stepindy-1] .= 0.0
    v_c[1:stepindx-1, 1:stepindy-1] .= 0.0

    # Normalize to unit vectors (direction only)
    mag = sqrt.(u_c .^ 2 .+ v_c .^ 2)
    mag[mag.==0.0] .= 1.0   # avoid division by zero in step region
    u_n = u_c ./ mag
    v_n = v_c ./ mag

    # Build arrow positions and directions
    xs = repeat(xc, 1, length(yc))        # (Nx, Ny)
    ys = repeat(yc', length(xc), 1)       # (Nx, Ny)

    fig = Figure(size=(1200, 500))
    ax = Axis(fig[1, 1],
        xlabel="x",
        ylabel="y",
        title="Velocity Direction — Backwards Facing Step",
        aspect=DataAspect()
    )

    arrows!(ax,
        vec(xs), vec(ys),
        vec(u_n), vec(v_n),
        arrowsize=8,
        lengthscale=0.4 * min(xc[2] - xc[1], yc[2] - yc[1]),
        arrowcolor=:black,
        linecolor=:black
    )

    # Draw step
    poly!(ax,
        Rect(xc[1], yc[1], xc[stepindx-1] - xc[1], yc[stepindy-1] - yc[1]),
        color=:gray30
    )

    display(fig)
    return fig
end

fig2 = plot_velocity_arrows(u, v, x, y, stepindx, stepindy)