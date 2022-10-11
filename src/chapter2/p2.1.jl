
module MID_21
  
using CairoMakie, DifferentialEquations

export sir21!, run_sir21, print_sir21, plot_sir21

"""
    sir21!(du, u, p, t) 

A simple compartmental susceptible--infectious--resistant model with a constant 
    population. This is the ordinary differential equations function for programme 2.1 in 
    `Modeling Infectious Diseases in Humans and Animals`

# Example 

    julia> u0 = [.999, .001, .0] # S, I, R
    julia> p = [1., .2] # beta, gamma
    julia> tspan = (0., 100.)
    julia> prob = ODEProblem(sir21!, u0, tspan, p)
    julia> sol = solve(prob)
"""
function sir21!(du, u, p, t) 
    # compartments 
    S, I, R = u
    # parameters
    beta, gamma = p 

    # the ODEs
    du[1] = dS = -beta * S * I 
    du[2] = dI = beta * S * I - gamma * I 
    du[3] = dR = gamma * I 
end 

"""
    run_sir21([; beta, gamma, S0, I0, duration, saveat])

Run the model `sir21`

# Keyword arguments 

All keyword arguments are optional with default values supplied for each. 

* `beta`: The beta parameter in the model (infectiousness of infectives). Default is 520 / 365 (520 per year).
* `gamma`: The gamma parameter in the model (recovery rate). Default is 1 / 7.
* `S0`: Proportion of the population susceptible at time = 0. Default is 1 - 1e-6.
* `I0`: Proportion of the population infectious at time = 0. Default is 1e-6.
* `duration`: How long the model will run (time units are interpretted as days). Default is 70.
* `saveat`: How frequently the model should save values. Default is 1 (day).

# Examples 
    julia> run_sir21()

    julia> run_sir21(beta = .8, gamma = .6, duration = 1000)
"""
function run_sir21(; beta = 520 / 365, gamma = 1 / 7, S0 = 1 - 1e-6, I0 = 1e-6, duration = 70, saveat = 1)
    @assert S0 >= 0 "Input S0 = $S0: cannot run with negative initial proportion susceptible"
    @assert I0 >= 0 "Input I0 = $I0: cannot run with negative initial proportion resistant"
    @assert beta >= 0 "Input beta = $beta: cannot run with negative parameter beta"
    @assert gamma >= 0 "Input gamma = $gamma: cannot run with negative parameter gamma"
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"
    @assert S0 + I0 <= 1 "Initial proportions susceptible ($S0) and infectious ($I0) must not sum to more than 1" 
    if beta < gamma @info "R₀ < 1 ($(beta / gamma))" end

    R_at_time0 = 1 - S0 - I0 # given a long name to ensure distinction between proportion 
        # resistant at time = 0 (R(0)) and the basic reproduction number (R₀)

    u0 = [S0, I0, R_at_time0]
    tspan = ( 0., Float64(duration) )
    p = [beta, gamma]
    prob = ODEProblem(sir21!, u0, tspan, p)
    sol = solve(prob; saveat)

    return sol
end 

"""
    print_sir21([; kwargs...])

Print the saved values after running the model `sir21`. 

Keyword arguments are all optional. See `run_sir21` for details of the arguments and their default values.
"""
function print_sir21(; kwargs...)
    sol = run_sir21(; kwargs...)
    print_sir21(sol)
end 

function print_sir21(sol)
    for i ∈ eachindex(sol.u)
        println("t = $(sol.t[i]): $(sol.u[i])")
    end 
    return nothing 
end 

"""
    plot_sir21([; kwargs...])
    plot_sir21(sol)

Plot the results of running the model `sir21`. 

Can take optional keyword arguments, `run_sir21` (see that function for details 
    of the arguments and their default values), or the solution from an ODE model. 
    The advantage of allowing the plot function to run the ODE is that `saveat` 
    is selected to provide a smooth line on the plot.

# Examples
    julia> plot_sir21()

    julia> plot_sir21(beta = .8, gamma = .6, duration = 100) 
"""
function plot_sir21(; duration = 70, kwargs...)
    saveat = duration / 500 # to give a smoother line for plotting
    sol = run_sir21(; duration, saveat, kwargs...)
    return plot_sir21(sol)
end 

function plot_sir21(sol)
    xs = sol.t ./ 7 # to plot time in weeks
    S = Float64[]; I = Float64[]; R = Float64[]
    for i ∈ eachindex(sol.u)
        push!(S, sol.u[i][1])
        push!(I, sol.u[i][2])
        push!(R, sol.u[i][3])
    end 

    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(ax, xs, S, label = "Susceptible")
    lines!(ax, xs, I, label = "Infectious")
    lines!(ax, xs, R, label = "Recovered")
    Label(
        fig[0, :], 
        "p2.1.jl: Susceptible--infectious--resistant model with a constant population"
    )
    ax.xlabel = "Time, weeks"
    ax.ylabel = "Fraction of population"
    fig[1, 2] = Legend(fig, ax)
    resize_to_layout!(fig)

    return fig
end 

end # module MID_21
