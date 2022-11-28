
module MID_26

# SEIR model (page 41)
  
using CairoMakie, DataFrames, DifferentialEquations

export seir26!, run_seir26, dataframe_seir26, plot_seir26, plot_seir26!

function seir26!(du, u, p, t)
    # compartments 
    S, E, I = u 
    # parameters 
    beta, gamma, mu, nu, sigma = p 

    # the ODEs
    du[1] = nu - beta * S * I - mu * S          # dS
    du[2] = beta * S * I - sigma * E - mu * E   # dE
    du[3] = sigma * E - gamma * I - mu * I      # dI
end 

function run_seir26(; S0, E0, I0, beta, gamma, mu, nu = mu, sigma, duration, kwargs...)
    u0 = [S0, E0, I0]
    p = [beta, gamma, mu, nu, sigma]
    return run_seir26(u0, p, duration; kwargs...)
end 

function run_seir26(u0, p, duration; reltol = 1e-6, saveat = 1)
    # set a low tolerance for the solver, otherwise the oscillations don't converge appropriately
    @assert minimum(u0) >= 0 "Input u0 = $u0: Cannot run model with negative starting values in any compartment"
    @assert sum(u0) <= 1 "Input u0 = $u0: Compartment values are proportions so must sum to 1 or less"
    @assert minimum(p) >= 0 "Input p = $p: Cannot run with negative values for any parameters"
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"

    tspan = ( 0., Float64(duration) )

    prob = ODEProblem(seir26!, u0, tspan, p)
    sol = solve(prob; reltol, saveat)
    return sol
end 

function dataframe_seir26(sol)
    result = DataFrame(t = Float64[], S = Float64[], E = Float64[], I = Float64[])
    for i ∈ eachindex(sol.u)
        push!( result, Dict(
            :t => sol.t[i], 
            :S => sol.u[i][1], 
            :E => sol.u[i][2], 
            :I => sol.u[i][3]
        ) )
    end 
    insertcols!(result, :_tempN => +(result.S, result.E, result.I) )
    insertcols!(result, :R => 1 .- result._tempN )
    select!(result, Not(:_tempN))
    return result 
end 

function plot_seir26(result; kwargs...)
    fig = Figure()
    plot_seir26!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_seir26!(any, sol; kwargs...) = plot_seir26!(any, dataframe_seir26(sol); kwargs...)

function plot_seir26!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_seir26!(gl, result; kwargs...)
end 

function plot_seir26!(gl::GridLayout, result::DataFrame; 
        label = "p2.6.jl: Susceptible--exposed--infectious--resistant model", 
        legend = :right, kwargs...
    )
    ax = Axis(gl[1, 1])
    plot_seir26!(ax, result; kwargs...)

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

function plot_seir26!(ax::Axis, result::DataFrame; plotR = true)
    xs = result.t ./ 365  # to plot time in years
    
    lines!(ax, xs, result.S, label = "Susceptible")
    lines!(ax, xs, result.E, label = "Exposed")
    lines!(ax, xs, result.I, label = "Infectious")
    plotR && lines!(ax, xs, result.R, label = "Recovered")
    ax.xlabel = "Time, years"
    ax.ylabel = "Fraction of population"
end 

end # module MID_26
