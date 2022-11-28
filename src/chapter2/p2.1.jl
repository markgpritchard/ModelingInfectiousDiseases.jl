
module MID_21

# Simple SIR model (page 19)
  
using CairoMakie, DataFrames, DifferentialEquations

export sir21!, run_sir21, dataframe_sir21, plot_sir21, plot_sir21!

function sir21!(du, u, p, t) 
    # compartments 
    S, I, R = u
    # parameters
    beta, gamma = p 

    # the ODEs
    du[1] = -beta * S * I               # dS
    du[2] = beta * S * I - gamma * I    # dI
    du[3] = gamma * I                   # dR
end 

function run_sir21(; S0, I0, R0 = 1 - (S0 + I0), beta, gamma, duration, kwargs...)
    u0 = [S0, I0, R0]
    p = [beta, gamma]
    return run_sir21(u0, p, duration; kwargs...)
end 

function run_sir21(u0, p, duration; saveat = 1)
    @assert minimum(u0) >= 0 "Input u0 = $u0: Cannot run model with negative starting values in any compartment"
    @assert sum(u0) ≈ 1 "Input u0 = $u0: Compartment values are proportions so must sum to 1"
    @assert minimum(p) >= 0 "Input p = $p: Cannot run with negative values for any parameters"
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"

    tspan = ( 0., Float64(duration) )
    prob = ODEProblem(sir21!, u0, tspan, p)
    sol = solve(prob; saveat)
    return sol
end 

function dataframe_sir21(sol)
    result = DataFrame(t = Float64[], S = Float64[], I = Float64[], R = Float64[])
    for i ∈ eachindex(sol.u)
        push!( result, Dict(
            :t => sol.t[i], 
            :S => sol.u[i][1], 
            :I => sol.u[i][2], 
            :R => sol.u[i][3]
        ) )
    end 
    return result 
end 

function plot_sir21(result; kwargs...)
    fig = Figure()
    plot_sir21!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_sir21!(any, sol; kwargs...) = plot_sir21!(any, dataframe_sir21(sol); kwargs...)

function plot_sir21!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_sir21!(gl, result; kwargs...)
end 

function plot_sir21!(gl::GridLayout, result::DataFrame; 
        label = "p2.1.jl: Susceptible--infectious--resistant model with a constant population", 
        legend = :right
    )
    ax = Axis(gl[1, 1])
    plot_sir21!(ax, result)

    Label(gl[0, :], label)

    if legend == :right
        leg = Legend(gl[1, 2], ax)
    elseif legend == :below 
        leg = Legend(gl[2, 1], ax)
    elseif legend == :none 
        # no legend 
    else 
        @info "Unrecognised legend position given. Recognised options are `:right`, `:below` and `:none`"
    end
end 

function plot_sir21!(ax::Axis, result::DataFrame)
    xs = result.t ./ 7  # to plot time in weeks
    
    lines!(ax, xs, result.S, label = "Susceptible")
    lines!(ax, xs, result.I, label = "Infectious")
    lines!(ax, xs, result.R, label = "Recovered")
    ax.xlabel = "Time, weeks"
    ax.ylabel = "Fraction of population"
end 

end # module MID_21
