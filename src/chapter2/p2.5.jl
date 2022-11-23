
module MID_25
  
using CairoMakie, DataFrames, DifferentialEquations

export sis25!, run_sis25, dataframe_sis25, plot_sis25, plot_sis25!

function sis25!(du, u, p, t)
    # compartments 
    S, I = u 
    # parameters 
    beta, gamma = p 

    # the ODEs
    du[1] = - beta * S * I + gamma * I  # dS
    du[2] = beta * S * I - gamma * I    # dI
end 

function run_sis25(u0, p, duration; saveat = 1)
    @assert minimum(u0) >= 0 "Input u0 = $u0: Cannot run model with negative starting values in any compartment"
    @assert sum(u0) ≈ 1 "Input u0 = $u0: Compartment values are proportions so must sum to 1"
    @assert minimum(p) >= 0 "Input p = $p: Cannot run with negative values for any parameters"
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"

    tspan = ( 0., Float64(duration) )

    prob = ODEProblem(sis25!, u0, tspan, p)
    sol = solve(prob; saveat)

    return sol
end 

function run_sis25(; I0, S0 = 1 - I0, beta, gamma, duration, kwargs...)
    u0 = [S0, I0]
    p = [beta, gamma]
    return run_sis25(u0, p, duration; kwargs...)
end 

function dataframe_sis25(sol)
    result = DataFrame(t = Float64[], S = Float64[], I = Float64[])
    for i ∈ eachindex(sol.u)
        push!( result, Dict(
            :t => sol.t[i], 
            :S => sol.u[i][1], 
            :I => sol.u[i][2], 
        ) )
    end 
    return result 
end 

function plot_sis25(result; kwargs...)
    fig = Figure()
    plot_sis25!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_sis25!(any, sol; kwargs...) = plot_sis25!(any, dataframe_sis25(sol); kwargs...)

function plot_sis25!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_sir22!(gl, result; kwargs...)
end 

function plot_sis25!(gl::GridLayout, result::DataFrame; 
        label = "p2.5.jl: Susceptible--infectious--susceptible model", 
        legend = :right, kwargs...
    )
    ax = Axis(gl[1, 1])
    plot_sis25!(ax, result; kwargs...)

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

function plot_sis25!(ax::Axis, result::DataFrame)
    xs = sol.t ./ 7 # to plot time in weeks
    
    lines!(ax, xs, result.S, label = "Susceptible")
    lines!(ax, xs, result.I, label = "Infectious")
    ax.xlabel = "Time, weeks"
    ax.ylabel = "Fraction of population"
end 

end # module MID_25
