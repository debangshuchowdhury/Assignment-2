begin
    using Plots
    x = 0:100
    y = sin.(x)
    plot(x, y; label="sine function")
    plot!(xlab="x", ylab="y")
    # print(x)
end


begin
    x = LinRange(0, 2, 20)
    print("x uncollected = ")
    println(x)
    for j in x
        print(j)
    end
    println()
    println()
    x_collected = collect(x)
    print("x collected = ")
    println(x_collected)
    for i in x_collected
        print(i)
    end
end

# DEFINING THE DOMAIN GEOMETRY
begin
    U_ref = 1
    L = 10
    h_step = 0.25
    H = 1
    Re = 100 # need to parameterise uref and re properly

    num_cells_y = 10
    dy = H / num_cells_y
    dx = dy
    num_cells_x = convert(Int32, L / dx)
    print("dy = ", dy)
    print("numcellsx = ", num_cells_x)
end

# plotting cell nodes

begin
    import Plots.PlotMeasures: px, mm
    # plotting the computational domain
    p = plot([0, 0], [0, H], color=:black)

    for j in -1:num_cells_y+1
        for i in -1:num_cells_x+1
            scatter!([i * dx], [j * dy], color=:black, markersize=:1, label=false)
        end
    end
    plot!(aspect_ratio=:equal)
    plot!(xlim=[-2 * dx, L + 2 * dx], ylim=[-2 * dy, H + 2 * dy])
    plot!(margin=2mm)
    display(p)
end
