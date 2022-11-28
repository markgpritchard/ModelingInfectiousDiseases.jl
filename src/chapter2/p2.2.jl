
module MID_22
  
using CairoMakie, DataFrames, DifferentialEquations

export sir22!, run_sir22, dataframe_sir22, plot_sir22, plot_sir22!

function sir22!(du, u, p, t) 
    # Compartments 
    S, I, R = u
    # Parameters
    beta, gamma, mu, nu = p 

    # The ODEs
    du[1] = nu - beta * S * I - mu * S          # dS
    du[2] = beta * S * I - gamma * I - mu * I   # dI 
    du[3] = gamma * I - mu * R                  # dR
end 

function run_sir22(u0, p, duration; reltol = 1e-12, saveat = 1)
    # set low tolerance for solver, otherwise oscillations don't converge appropriately
    @assert minimum(u0) >= 0 "Input u0 = $u0: Cannot run model with negative starting values in any compartment"
    @assert sum(u0) ≈ 1 "Input u0 = $u0: Compartment values are proportions so must sum to 1"
    @assert minimum(p) >= 0 "Input p = $p: Cannot run with negative values for any parameters"
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"

    tspan = ( 0., Float64(duration) )

    prob = ODEProblem(sir22!, u0, tspan, p)
    sol = solve(prob; reltol, saveat)
    return sol
end 

function run_sir22(; S0, I0, R0 = 1 - (S0 + I0), beta, gamma, mu, nu = mu, duration, kwargs...)
    u0 = [S0, I0, R0]
    p = [beta, gamma, mu, nu]
    return run_sir22(u0, p, duration; kwargs...)
end 

function dataframe_sir22(sol)
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

function plot_sir22(result; kwargs...)
    fig = Figure()
    plot_sir22!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_sir22!(any, sol; kwargs...) = plot_sir22!(any, dataframe_sir22(sol); kwargs...)

function plot_sir22!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_sir22!(gl, result; kwargs...)
end 

function plot_sir22!(gl::GridLayout, result::DataFrame; 
        label = "p2.2.jl: Susceptible--infectious--resistant model with a variable population", 
        legend = :right, kwargs...
    )
    ax = Axis(gl[1, 1])
    plot_sir22!(ax, result; kwargs...)

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

function plot_sir22!(ax::Axis, result::DataFrame; plotR = true)
    xs = result.t ./ 365  # to plot time in years
    
    lines!(ax, xs, result.S, label = "Susceptible")
    lines!(ax, xs, result.I, label = "Infectious")
    plotR && lines!(ax, xs, result.R, label = "Recovered")
    ax.xlabel = "Time, years"
    ax.ylabel = "Fraction of population"
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Help text 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


"""
    run_sir22([; beta, gamma, mu, S0, I0, duration, saveat])

Run the model `sir22`

# Keyword arguments 

All keyword arguments are optional with default values supplied for each. 

* `beta`: The beta parameter in the model (infectiousness of infectives). Default is 520 / 365 (520 per year).
* `gamma`: The gamma parameter in the model (recovery rate). Default is 1 / 7.
* `mu`: The model's birth and mortality rate. Defaults is 1 / (70 * 365) 
    (i.e. mean life duration 70 years and birth rate to match mortality rate).
* `S0`: Proportion of the population susceptible at time = 0. Default is 0.1.
* `I0`: Proportion of the population infectious at time = 0. Default is 1e-4.
* `duration`: How long the model will run (time units are interpretted as days). Default is 60 years.
* `saveat`: How frequently the model should save values. Default is 1 (day).
"""

end # module MID_22
