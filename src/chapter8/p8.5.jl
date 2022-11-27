
module MID_85
  
using CairoMakie, DataFrames, DifferentialEquations

export siqr85!, run_siqr85, dataframe_siqr85, plot_siqr85, plot_siqr85!

function siqr85!(du, u, p, t) 
    # compartments 
    X, Xq, Y, Q, Z = u
    N = +(X, Xq, Y, Q, Z)

    # parameters 
    b, k, di, q, tauq, gamma = p 

    # the ODEs
    du[1] = - (k * b * Y + q * k * (1 - b) * Y) * X / N + tauq * Xq     # dX 
    du[2] = q * k * (1 - b) * Y * X / N - tauq * Xq                     # dXq 
    du[3] = k * b * Y * (1 - q) * X / N - (di + gamma) * Y              # dY 
    du[4] = q * k * b * X * Y / N + di * Y - tauq * Q                   # dQ 
    du[5] = gamma * Y + tauq * Q                                        # dZ
end 

function run_siqr85(u0, p, duration; saveat = 1)
    @assert minimum(u0) >= 0 "Input u0 = $u0: no compartments can contain negative proportions"
    @assert minimum(p) >= 0 "Input p = $p: cannot run with negative parameters"
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"

    tspan = ( 0., Float64(duration) )
    prob = ODEProblem(siqr85!, u0, tspan, p)
    return solve(prob; saveat)
end 

function run_siqr85(; X0, Xq0, Y0, Q0, Z0, b, k, di, q, tauq, gamma, duration, kwargs...)
    u0 = [X0, Xq0, Y0, Q0, Z0]
    p = [b, k, di, q, tauq, gamma]
    return run_siqr85(u0, p, duration; kwargs...)
end

function dataframe_siqr85(sol)
    return DataFrame(
        t = sol.t,
        X = [ sol[i][1] for i ∈ axes(sol, 2) ],
        Xq = [ sol[i][2] for i ∈ axes(sol, 2) ],
        Y = [ sol[i][3] for i ∈ axes(sol, 2) ],
        Q = [ sol[i][4] for i ∈ axes(sol, 2) ],
        Z = [ sol[i][5] for i ∈ axes(sol, 2) ],
    ) 
end 

function plot_siqr85(result; kwargs...)
    fig = Figure()
    plot_siqr85!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_siqr85!(any, sol; kwargs...) = plot_siqr85!(any, dataframe_sir84(sol); kwargs...)

function plot_siqr85!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_siqr85!(gl, result; kwargs...)
end 

function plot_siqr85!(gl::GridLayout, result::DataFrame; 
        label = "p8.5.jl: SIR model with quarantine"
    )
    ax = Axis(gl[1, 1])
    plot_siqr85!(ax, result)
    ax.xlabel = "Time"
    ax.ylabel = "Number"
    leg = Legend(gl[2, 1], ax; orientation = :horizontal)
    Label(fig85[0, 1], label)
end 

function plot_siqr85!(ax::Axis, result::DataFrame)
    lines!(ax, result.t, result.X; label = "Susceptible")
    lines!(ax, result.t, result.Xq; label = "Quarantined contacts")
    lines!(ax, result.t, result.Y; label = "Infected")
    lines!(ax, result.t, result.Q; label = "Quarantined cases")
    lines!(ax, result.t, result.Z; label = "Recovered")
end 

end # module MID_85
