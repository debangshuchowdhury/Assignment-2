begin
    using Pkg
    Pkg.activate("solvers/")
    using LinearAlgebra
    using solvers
    using Statistics
    using Plots
    pythonplot()

    N = 200
    dx = 1.0 / N
    dy = 1.0 / N

    # cell-center coordinates
    x = [(i - 0.5) * dx for i in 1:N]
    y = [(j - 0.5) * dy for j in 1:N]

    # manufactured solution and source term on cell centers
    # u_exact = [cos(π * x[i]) * cos(π * y[j]) for i in 1:N, j in 1:N]
    # f = [2π^2 * cos(π * x[i]) * cos(π * y[j]) for i in 1:N, j in 1:N]
    f = zeros(N, N)

    u_num = poisson_solver([N, N], f, dx, dy, 1e-8)
    # display(u_num)

    heatmap(u_num')
    # surface(x, y, u_num)

end


# #  First Print
# print("helo world")


# # Variable assignment
# a = 10
# b = 1 // 123


# # Type() in Julia
# typeof(a)


# # typecasting can be done with convert()
# convert(Float16, a)


# # defining a function
# function addition(a, b)
#     return a + b
# end


# # inline definition of function
# addition_inline(a, b) = a + b


# # anonymous functions. they work like the lambda keyword in python. written like a map in maths
# addition_anon = (a, b) -> a + b


# # integration using anonymous functions
# using Pkg
# Pkg.add("QuadGK")
# begin
#     using QuadGK
#     f(x, y, z) = (x^2 + 2y) * z
#     quadgk(x -> f(x, 42, 4), 3, 4)
# end
# # alternatively can write the integration as
# begin
#     using QuadGK
#     f(x, y, z) = (x^2 + 2y) * z
#     arg(x) = f(x, 42, 4)
#     quadgk(arg, 3, 4)
# end
# # basically avoid the usage of anaonyomus functions as much as possible


# # keyword arguments are separated from positional arguments using a semicolon. 
# # to enforce the usage of keywords when the number or arguments is large
# function add(; a, b)
#     return a + b
# end
# add(a=1, b=1)
# # positional arguments are faster than keyword arguments


# # the if else block is identical
# # and is written as &, or is written as ||


# # string interpolation by preceding variable name with $. also using let deletes from memory after exiting the func
# let
#     name1 = "deb"
#     name2 = "debangshu"
#     println("$name1 and $name2 and \$")
# end


# # splatting: adding a ... after a tuble automatically unpacks integration
# begin
#     println((1, 2, 3))
#     println((1, 2, 3)...)
# end


# # named tuples are like a combination of tuples and dictionaries
# begin
#     a1 = (a=2, b=3)
#     a1[:a]
# end


# # dictionaries are defined using =>
# let
#     dictin = Dict("name" => "debangshu", "huh" => "what")
#     dictin
# end


# # broadcasting means you apply an operation or map element by element for an array or matrix
# let
#     a = [1, 2, 3]
#     sin.(a)
# end


# # defining a module
# module attempt
# export qwerty
# qwerty = 12345
# asdf = 0
# end
# let
#     using .attempt
#     println(qwerty)
#     println(attempt.asdf)
# end


# # 'using' imports all the variables or functions that are exported by the module
# # 'import' imports the module but does not import the funcs and variables. they need to be called with modul_name.variable_name
# begin
#     using SpecialFunctions
#     gamma(3)
# end
# begin
#     import SpecialFunctions
#     SpecialFunctions.gamma(3)
# end
# begin
#     using SpecialFunctions: gamma, sinint
#     gamma(3)
# end