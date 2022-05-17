
module MID_2_2
  
using CairoMakie, DifferentialEquations, StaticArrays

export sir_22, run_sir_22, print_sir_22, plot_sir_22

"""
    sir_22(u, p, t) 

A compartmental susceptible--infectious--resistant model with a variable 
    population. This is the ordinary differential equations function for programme 2.2 in 
    `Modeling Infectious Diseases in Humans and Animals`
"""
function sir_22(u, p, t)
    S, I, R = u
    beta, gamma, mu = p 

    dS = mu - beta * S * I - mu * S 
    dI = beta * S * I - gamma * I - mu * I
    dR = gamma * I - mu * R
    return @SVector [dS, dI, dR]
end 

"""
    run_sir_22([; beta, gamma, mu, S0, I0, duration, saveat])

Run the model `sir_22`

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
function run_sir_22(; beta = 520 / 365, gamma = 1 / 7, mu = 1 / (70 * 365), S0 = .1, I0 = 1e-4, 
        duration = 60 * 365, saveat = 1)
  
    @assert S0 >= 0 "Input S0 = $S0: cannot run with negative initial proportion susceptible"
    @assert I0 >= 0 "Input I0 = $I0: cannot run with negative initial proportion resistant"
    @assert beta >= 0 "Input beta = $beta: cannot run with negative parameter beta"
    @assert gamma >= 0 "Input gamma = $gamma: cannot run with negative parameter gamma"
    @assert mu >= 0 "Input mu = $mu: cannot run with negative parameter mu"
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"
    @assert S0 + I0 <= 1 "Initial proportions susceptible ($S0) and infectious ($I0) must not sum to more than 1" 
    if beta < gamma @info "R₀ < 1 ($(beta / gamma))" end

    R_at_time0 = 1 - S0 - I0 # given a long name to ensure distinction between proportion 
        # resistant at time = 0 (R(0)) and the basic reproduction number (R₀)

    u0 = @SVector [S0, I0, R_at_time0]
    tspan = ( 0., Float64(duration) )
    p = [beta, gamma, mu]
    prob = ODEProblem(sir_22, u0, tspan, p)
    # set a low tolerance for the solver, otherwise the oscillations don't converge appropriately
    sol = solve(prob; saveat, abstol = 1e-12, reltol = 1e-12)

    return sol
end 

"""
    print_sir_22([; kwargs...])

Print the saved values after running the model `sir_22`. 

Keyword arguments are all optional. See `run_sir_22` for details of the arguments and their default values.
"""
function print_sir_22(; kwargs...)
    sol = run_sir_22(; kwargs...)
    for i in eachindex(sol.u)
        println(sol.u[i])
    end 
    return sol 
end 

"""
    plot_sir_22([; kwargs...])
    plot_sir_22(sol)

Plot the results of running the model `sir_22`. 

Can take optional keyword arguments, `run_sir_22` (see that function for details 
    of the arguments and their default values), or the solution from an ODE model. 
    The advantage of allowing the plot function to run the ODE is that `saveat` 
    is selected to provide a smooth line on the plot.
"""
function plot_sir_22(; kwargs...)
    sol = run_sir_22(; kwargs...)
    return plot_sir_22(sol)
end 

function plot_sir_22(sol)
  xs = sol.t ./ 365 # to plot time in years
  S = Float64[]; I = Float64[]; R = Float64[]
  for i in eachindex(sol.u)
      push!(S, sol.u[i][1])
      push!(I, sol.u[i][2])
      push!(R, sol.u[i][3])
  end 

  fig = Figure()
  ax1 = Axis(fig[1, 1]); ax2 = Axis(fig[2, 1]); ax3 = Axis(fig[3, 1])
  lines!(ax1, xs, S)
  lines!(ax2, xs, I)
  lines!(ax3, xs, R)
  ax3.xlabel = "Time, years"
  ax1.ylabel = "Fraction susceptible"
  ax2.ylabel = "Fraction infectious"
  ax3.ylabel = "Fraction resistant"
  hidexdecorations!(ax1; grid = false); hidexdecorations!(ax2; grid = false)

  return fig
end 

end # module MID_2_2
