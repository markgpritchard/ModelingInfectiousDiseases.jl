
module MID_25
  
using CairoMakie, DifferentialEquations

export plot_sis25, print_sis25, run_sis25, sis25!

"""
    sis25!(du, u, p, t) 

A compartmental susceptible--infectious--susceptible model.
    
This is the ordinary differential equations function for programme 2.5 in 
    `Modeling Infectious Diseases in Humans and Animals`
"""
function sis25!(du, u, p, t)
    # compartments 
    S, I = u 
    # parameters 
    beta, gamma = p 

    # the ODEs
    du[1] = dS = - beta * S * I + gamma * I
    du[2] = dI = beta * S * I - gamma * I
end 

"""
    run_sis25([; beta, gamma, I0, duration, saveat])

Run the model `sis25`

# Keyword arguments 

All keyword arguments are optional with default values supplied for each. 

* `beta`: The beta parameter in the model (infectiousness of infectives). Default is 520 / 365 (520 per year).
* `gamma`: The gamma parameter in the model (recovery rate). Default is 1 / 7.
* `I0`: Proportion infectious at time = 0. Default is 1e-6.
* `duration`: How long the model will run (time units are interpretted as days). 
    Default is 70 days
* `saveat`: How frequently the model should save values. Default is 1 (day).
"""
function run_sis25(; beta = 520 / 365, gamma = 1 / 7, I0 = 1e-6, duration = 70, saveat = 1)
    @assert I0 >= 0 "Input Y0 = $I0: cannot run with negative initial proportion infectious"
    @assert beta >= 0 "Input beta = $beta: cannot run with negative parameter beta"
    @assert gamma >= 0 "Input gamma = $gamma: cannot run with negative parameter gamma"
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"
    @assert I0 <= 1 "Initial propotion infectious ($I0) cannot be greater than 1" 
    if beta < gamma @info "R₀ < 1 ($(beta / gamma))" end

    S0 = 1 - I0
    u0 = [S0, I0]
    tspan = ( 0., Float64(duration) )
    p = [beta, gamma]
    prob = ODEProblem(sis25!, u0, tspan, p)
    # set a low tolerance for the solver, otherwise the oscillations don't converge appropriately
    sol = solve(prob; saveat)

    return sol
end 

"""
    print_sis25([; kwargs...])

Print the saved values after running the model `sis25`. 

Keyword arguments are all optional. See `run_sis25` for details of the arguments and their default values.
"""
function print_sis25(; kwargs...)
    sol = run_sis25(; kwargs...)
    print_sis25(sol)
end 

function print_sis25(sol)
    for i ∈ eachindex(sol.u)
        println("t = $(sol.t[i]): $(sol.u[i])")
    end 
    return nothing 
end 

"""
    plot_sis25([; kwargs...])
    plot_sis25(sol)

Plot the results of running the model `sis25`. 

Can take optional keyword arguments, `run_sis25` (see that function for details 
    of the arguments and their default values), or the solution from an ODE model. 
    The advantage of allowing the plot function to run the ODE is that `saveat` 
    is selected to provide a smooth line on the plot.
"""
function plot_sis25(; duration = 70, kwargs...)
    saveat = duration / 500 # to give a smoother line for plotting
    sol = run_sis25(; duration, saveat, kwargs...)
    return plot_sis25(sol)
end 

function plot_sis25(sol)
    xs = sol.t ./ 7 # to plot time in weeks
    S = Float64[]; I = Float64[]
    for i ∈ eachindex(sol.u)
        push!(S, sol.u[i][1])
        push!(I, sol.u[i][2])
    end 

    fig = Figure()
    ax = Axis(fig[1, 1])
    lines!(ax, xs, S, label = "Susceptible")
    lines!(ax, xs, I, label = "Infectious")
    Label(fig[0, :], "p2.5.jl: Susceptible--infectious--susceptible model")
    ax.xlabel = "Time, weeks"
    ax.ylabel = "Fraction of population"
    fig[1, 2] = Legend(fig, ax)
    resize_to_layout!(fig)

    return fig
end 

end # module MID_25
