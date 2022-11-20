
module MID_78
  
using CairoMakie, DataFrames, DifferentialEquations

export sis78!, u0_sis78, run_sis78, dataframe_sis78, plot_sis78, plot_sis78!

function u0_sis78(X0, Y0, n) 
    N0 = X0 + Y0
    XY0 = n * Y0 * X0 / N0
    u0 = [ X0, XY0, N0 ]
    return u0
end 

function sis78!(du, u, p, t) 
    # compartments 
    X, XY, N = u

    # parameters 
    gamma, tau, n = p
    
    # ODEs
    # dX
    du[1] = gamma * (N - X) - tau * XY

    # dXY
    du[2] = tau * (n - 1) * (n * X - XY) * XY / (n * X) + gamma * (n * N - n * X - XY) - 
        tau * XY - tau * (n - 1) * XY^2 / (n * X) - gamma * XY 

    # N is constant so dN = 0
    du[3] = 0
end 

function run_sis78(u0, p, duration; saveat = 1)
    @assert minimum(u0) >= 0 "Input u0 = $u0: no compartments can contain negative proportions"
    @assert minimum(p) >= 0 "Input p = $p: cannot run with negative parameters" 
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"

    tspan = ( 0., Float64(duration) )

    prob = ODEProblem(sis78!, u0, tspan, p)
    sol = solve(prob; saveat)

    return sol
end 

function dataframe_sis78(sol)
    result = DataFrame(t = sol.t)
    insertcols!(result, :X => [ sol[i][1] for i ∈ axes(sol, 2) ])
    insertcols!(result, :Y => [ sol[i][3] - sol[i][1] for i ∈ axes(sol, 2) ])
    insertcols!(result, :XY => [ sol[i][2] for i ∈ axes(sol, 2) ])
    return result 
end 

function plot_sis78(result; kwargs...)
    fig = Figure()
    plot_sis78!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_sis78!(any, sol; kwargs...) = plot_sis78!(any, dataframe_sis78(sol); kwargs...)

function plot_sis78!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_sis78!(gl, result; kwargs...)
end 

function plot_sis78!(gl::GridLayout, result::DataFrame; 
        label = "p7.8.jl: Pairwise SIS model")

    ax1 = Axis(gl[1, 1])
    lines!(ax1, result.t, result.X; label = "Susceptible")
    lines!(ax1, result.t, result.Y; label = "Infectious")
    ax1.ylabel = "Individuals"
    leg = Legend(gl[1, 2], ax1)
    hidexdecorations!(ax1; grid = false, ticks = false)

    ax2 = Axis(gl[2, 1])
    lines!(ax2, result.t, result.XY)
    ax2.ylabel = "XY pairs" 
    ax2.xlabel = "Time"
    linkaxes!(ax1, ax2)

    Label(gl[0, :], label)
end 

end # module MID_78
