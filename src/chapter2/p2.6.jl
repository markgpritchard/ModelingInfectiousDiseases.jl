
module MID_26
  
using CairoMakie, DifferentialEquations

export plot_seir26, print_seir26, run_seir26, seir26!

"""
    seir26!(du, u, p, t) 

A compartmental susceptible--exposed--infectious--resistant model.
    
This is the ordinary differential equations function for programme 2.6 in 
    `Modeling Infectious Diseases in Humans and Animals`
"""
function seir26!(du, u, p, t)
    # compartments 
    S, E, I = u 
    # parameters 
    beta, gamma, mu, sigma = p 

    # the ODEs
    du[1] = dS = mu - beta * S * I - mu * S
    du[2] = dE = beta * S * I - sigma * E - mu * E
    du[3] = dI = sigma * E - gamma * I - mu * I
end 

"""
    run_seir26([; beta, gamma, I0, duration, saveat])

Run the model `seir26`

# Keyword arguments 

All keyword arguments are optional with default values supplied for each. 

* `beta`: The beta parameter in the model (infectiousness of infectives). Default is 520 / 365 (520 per year).
* `gamma`: The gamma parameter in the model (recovery rate). Default is 1 / 7.
* `mu`: Per capita birth and death rate. Default is 1 / 70 years.
* `sigma`: The rate at which exposed individuals become infectious. Default is 1 / 14.
* `S0`: Proportion susceptible at time = 0. Default is 0.1.
* `E0`: Proportion exposed at time = 0. Default is 1e-4.
* `I0`: Proportion infectious at time = 0. Default is 1e-4.
* `duration`: How long the model will run (time units are interpretted as days). 
    Default is 60 years
* `saveat`: How frequently the model should save values. Default is 1 (day).
"""
function run_seir26(; beta = 520 / 365, gamma = 1 / 7, mu = 1 / (70 * 365), sigma = 1 / 14, 
        S0 = .1, E0 = 1e-4, I0 = 1e-4, duration = 60 * 365, saveat = 1)
    
    @assert S0 >= 0 "Input S0 = $S0: cannot run with negative initial proportion susceptible"
    @assert E0 >= 0 "Input E0 = $E0: cannot run with negative initial proportion exposed"
    @assert I0 >= 0 "Input I0 = $I0: cannot run with negative initial proportion infectious"
    @assert beta >= 0 "Input beta = $beta: cannot run with negative parameter beta"
    @assert gamma >= 0 "Input gamma = $gamma: cannot run with negative parameter gamma"
    @assert mu >= 0 "Input mu = $mu: cannot run with negative parameter mu"
    @assert sigma >= 0 "Input sigma = $sigma: cannot run with negative parameter sigma"
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"
    @assert S0 + E0 + I0 <= 1 "Initial propotions susceptible ($S0), exposed ($E0) and infectious ($I0) cannot sum to greater than 1" 
    if beta * sigma < ( (gamma + mu) * (sigma + mu) ) 
        @info "R₀ < 1 ($( beta * sigma / ( (gamma + mu) * (sigma + mu) ) ))" 
    end

    u0 = [S0, E0, I0]
    tspan = ( 0., Float64(duration) )
    p = [beta, gamma, mu, sigma]
    prob = ODEProblem(seir26!, u0, tspan, p)
    # set a low tolerance for the solver, otherwise the oscillations don't converge appropriately
    sol = solve(prob; saveat, reltol = 1e-6)

    return sol
end 

"""
    print_seir26([; kwargs...])

Print the saved values after running the model `seir26`. 

Keyword arguments are all optional. See `run_seir26` for details of the arguments and their default values.
"""
function print_seir26(; kwargs...)
    sol = run_seir26(; kwargs...)
    print_seir26(sol)
end 

function print_seir26(sol)
    for i ∈ eachindex(sol.u)
        println("t = $(sol.t[i]): $(sol.u[i])")
    end 
    return nothing 
end 

"""
    plot_seir26([; kwargs...])
    plot_seir26(sol)

Plot the results of running the model `seir26`. 

Can take optional keyword arguments, `run_seir26` (see that function for details 
    of the arguments and their default values), or the solution from an ODE model. 
    The advantage of allowing the plot function to run the ODE is that `saveat` 
    is selected to provide a smooth line on the plot.
"""
function plot_seir26(; kwargs...)
    sol = run_seir26(; kwargs...)
    return plot_seir26(sol)
end 

function plot_seir26(sol)
    xs = sol.t ./ 365 # to plot time in years
    S = Float64[]; E = Float64[]; I = Float64[]; R = Float64[]
    for i ∈ eachindex(sol.u)
        push!(S, sol.u[i][1])
        push!(E, sol.u[i][2])
        push!(I, sol.u[i][3])
        push!(R, 1 - (S[i] + E[i] + I[i]))
    end 

    fig = Figure()
    ax1 = Axis(fig[1, 1]); ax2 = Axis(fig[2, 1]); ax3 = Axis(fig[3, 1])
    lines!(ax1, xs, S)
    lines!(ax2, xs, E, label = "Exposed")
    lines!(ax2, xs, I, label = "Infectious")
    lines!(ax3, xs, R)
    Label(fig[0, :], "p2.6.jl: Susceptible--exposed--infectious--resistant model")
    ax3.xlabel = "Time, years"
    ax1.ylabel = "Proportion susceptible"
    ax2.ylabel = "Fraction of population"
    ax3.ylabel = "Proportion resistant"
    fig[2, 2] = Legend(fig, ax2)
    linkxaxes!(ax1, ax2, ax3)
    hidexdecorations!(ax1; grid = false, ticks = false)
    hidexdecorations!(ax2; grid = false, ticks = false)
    resize_to_layout!(fig)

    return fig
end 

end # module MID_26
