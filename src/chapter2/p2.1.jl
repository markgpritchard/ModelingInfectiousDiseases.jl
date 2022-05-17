
module MID_2_1
  
using CairoMakie, DifferentialEquations, StaticArrays

export sir_21, run_sir_21, print_sir_21, plot_sir_21

"""
    sir_21(u, p, t) 

A simple compartmental susceptible--infectious--resistant model with a constant 
    population. This is the ordinary differential equations function for programme 2.1 in 
    `Modeling Infectious Diseases in Humans and Animals`
"""
function sir_21(u, p, t)
    # compartments 
    S, I, R = u
    # parameters
    beta, gamma = p 

    dS = -beta * S * I 
    dI = beta * S * I - gamma * I 
    dR = gamma * I 
    return @SVector [dS, dI, dR]
end 

"""
    run_sir_21([; beta, gamma, S0, I0, duration, saveat])

Run the model `sir_21`

# Keyword arguments 

All keyword arguments are optional with default values supplied for each. 

* `beta`: The beta parameter in the model (infectiousness of infectives). Default is 520 / 365 (520 per year).
* `gamma`: The gamma parameter in the model (recovery rate). Default is 1 / 7.
* `S0`: Proportion of the population susceptible at time = 0. Default is 1 - 1e-6.
* `I0`: Proportion of the population infectious at time = 0. Default is 1e-6.
* `R_at_time0`: Proportion of the population resistant at time = 0. Has a long 
    name to avoid confusion with the basic reproduction number R₀. Default value is 0.
* `duration`: How long the model will run (time units are interpretted as days). Default is 70.
* `saveat`: How frequently the model should save values. Default is 1 (day).
"""
function run_sir_21(; beta = 520 / 365, gamma = 1 / 7, S0 = 1 - 1e-6, I0 = 1e-6, 
        R_at_time0 = 0, duration = 70, saveat = 1)
  
    @assert S0 >= 0 "Input S0 = $S0: cannot run with negative initial proportion susceptible"
    @assert I0 >= 0 "Input I0 = $I0: cannot run with negative initial proportion resistant"
    @assert beta >= 0 "Input beta = $beta: cannot run with negative parameter beta"
    @assert gamma >= 0 "Input gamma = $gamma: cannot run with negative parameter gamma"
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"
    @assert S0 + I0 <= 1 "Initial proportions susceptible ($S0) and infectious ($I0) must not sum to more than 1" 
    if beta < gamma @info "R₀ < 1 ($(beta / gamma))" end

    R_at_time0 = 1 - S0 - I0 # given a long name to ensure distinction between proportion 
        # resistant at time = 0 (R(0)) and the basic reproduction number (R₀)

    u0::SVector{3, Float64} = @SVector [S0, I0, R_at_time0]
    tspan::Tuple{Float64, Float64} = ( 0., Float64(duration) )
    p::Vector{Float64} = [beta, gamma]
    prob = ODEProblem(sir_21, u0, tspan, p)
    sol = solve(prob; saveat)

    return sol
end 

"""
    print_sir_21([; kwargs...])

Print the saved values after running the model `sir_21`. 

Keyword arguments are all optional. See `run_sir_21` for details of the arguments and their default values.
"""
function print_sir_21(; kwargs...)
    sol = run_sir_21(; kwargs...)
    for i in eachindex(sol.u)
        println(sol.u[i])
    end 
    return sol 
end 

"""
    plot_sir_21([; kwargs...])
    plot_sir_21(sol)

Plot the results of running the model `sir_21`. 

Can take optional keyword arguments, `run_sir_21` (see that function for details 
    of the arguments and their default values), or the solution from an ODE model. 
    The advantage of allowing the plot function to run the ODE is that `saveat` 
    is selected to provide a smooth line on the plot.
"""
function plot_sir_21(; duration = 70, kwargs...)
    saveat = duration / 500 # to give a smoother line for plotting
    sol = run_sir_21(; duration, saveat, kwargs...)
    return plot_sir_21(sol)
end 

function plot_sir_21(sol)
  xs = sol.t ./ 7 # to plot time in weeks
  S = Float64[]; I = Float64[]; R = Float64[]
  for i in eachindex(sol.u)
      push!(S, sol.u[i][1])
      push!(I, sol.u[i][2])
      push!(R, sol.u[i][3])
  end 

  fig = Figure()
  ax = Axis(fig[1, 1])
  lines!(ax, xs, S, label = "Susceptible")
  lines!(ax, xs, I, label = "Infectious")
  lines!(ax, xs, R, label = "Recovered")
  ax.xlabel = "Time, weeks"
  ax.ylabel = "Fraction of population"
  fig[1, 2] = Legend(fig, ax)

  return fig
end 

end # module MID_2_1
