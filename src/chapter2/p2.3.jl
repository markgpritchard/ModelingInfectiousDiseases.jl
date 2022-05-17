
module MID_2_3
  
using CairoMakie, DifferentialEquations, StaticArrays

export sir_23, run_sir_23, print_sir_23, plot_sir_23

"""
    sir_23(u, p, t) 

A compartmental susceptible--infectious--resistant model with infection-induced
    mortality. This is the ordinary differential equations function for programme 2.3 in 
    `Modeling Infectious Diseases in Humans and Animals`
"""
function sir_23(u, p, t)
    X, Y, Z = u # following convention in book, {S, I, R} refer to proportions and {X, Y, Z} refer to numbers
    beta, gamma, mu, nu, rho = p 

    dX = nu - beta * X * Y - mu * X 
    dY = beta * X * Y - (gamma + mu) * Y / (1 - rho)
    dZ = gamma * Y - mu * Z
    return @SVector [dX, dY, dZ]
end 

"""
    run_sir_23([; beta, gamma, mu, nu, rho, X0, Y0, N0, duration, saveat])

Run the model `sir_23`

# Keyword arguments 

All keyword arguments are optional with default values supplied for each. 

* `beta`: The beta parameter in the model (infectiousness of infectives). Default is 520 / 365 (520 per year).
* `gamma`: The gamma parameter in the model (recovery rate). Default is 1 / 7.
* `mu`: The model's mortality rate not due to the pathogen. Default is 1 / (70 * 365) 
    (i.e. mean life duration 70 years).
* `nu`: The model's birth rate. Default is 1 / (70 * 365) (i.e. match the mortality in the absence of the pathogen).
* `rho`: The mortality probability for infecteds. Default is 0.5.
* `X0`: Number susceptible at time = 0. Default is 0.2. Note, we follow the convention 
    used in the book that {S, I, R} represent proportions and {X, Y, Z} represent 
    numbers susceptible, infectious and resistant.
* `Y0`: Number infectious at time = 0. Default is 1e-6.
* `N0`: The initial total population size. Default is 1.
* `duration`: How long the model will run (time units are interpretted as days). 
    Default is 100 years. (The default differs between the programmes in different 
    languages -- 100 years appears to give a reasonable illustration of the results.)
* `saveat`: How frequently the model should save values. Default is 1 (day).
"""
function run_sir_23(; beta = 520 / 365, gamma = 1 / 7, mu = 1 / (70 * 365), nu = 1 / (70 * 365), 
        rho = .5, X0 = .2, Y0 = 1e-6, N0 = 1, duration = 100 * 365, saveat = 1)
  
    @assert X0 >= 0 "Input X0 = $X0: cannot run with negative initial number susceptible"
    @assert Y0 >= 0 "Input Y0 = $Y0: cannot run with negative initial number resistant"
    @assert beta >= 0 "Input beta = $beta: cannot run with negative parameter beta"
    @assert gamma >= 0 "Input gamma = $gamma: cannot run with negative parameter gamma"
    @assert mu >= 0 "Input mu = $mu: cannot run with negative parameter mu"
    @assert nu >= 0 "Input nu = $nu: cannot run with negative parameter nu"
    @assert rho >= 0 && rho <= 1 "Input rho = $rho: rho is a probability and cannot be <0 or >1"
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"
    @assert X0 + Y0 <= N0 "Initial numbers susceptible ($X0) and infectious ($Y0) must not sum to more than the total population size ($N0)" 
    if beta * (1 - rho) * nu < (gamma + nu) * mu @info "R₀ < 1 ($(beta * (1 - rho) * nu / (gamma + nu) * mu))" end

    Z0 = N0 - X0 - Y0 # given a long name to ensure distinction between proportion 
        # resistant at time = 0 (R(0)) and the basic reproduction number (R₀)

    u0 = @SVector [X0, Y0, Z0]
    tspan = ( 0., Float64(duration) )
    p = [beta, gamma, mu, nu, rho]
    prob = ODEProblem(sir_23, u0, tspan, p)
    # set a low tolerance for the solver, otherwise the oscillations don't converge appropriately
    sol = solve(prob; saveat, abstol = 1e-12, reltol = 1e-12)

    return sol
end 

"""
    print_sir_23([; kwargs...])

Print the saved values after running the model `sir_23`. 

Keyword arguments are all optional. See `run_sir_23` for details of the arguments and their default values.
"""
function print_sir_23(; kwargs...)
    sol = run_sir_23(; kwargs...)
    for i in eachindex(sol.u)
        println(sol.u[i])
    end 
    return sol 
end 

"""
    plot_sir_23([; kwargs...])
    plot_sir_23(sol)

Plot the results of running the model `sir_23`. 

Can take optional keyword arguments, `run_sir_23` (see that function for details 
    of the arguments and their default values), or the solution from an ODE model. 
    The advantage of allowing the plot function to run the ODE is that `saveat` 
    is selected to provide a smooth line on the plot.
"""
function plot_sir_23(; kwargs...)
    sol = run_sir_23(; kwargs...)
    return plot_sir_23(sol)
end 

function plot_sir_23(sol)
  xs = sol.t ./ 365 # to plot time in years
  X = Float64[]; Y = Float64[]; Z = Float64[]; N = Float64[]
  for i in eachindex(sol.u)
      push!(X, sol.u[i][1])
      push!(Y, sol.u[i][2])
      push!(Z, sol.u[i][3])
      push!(N, (X[i] + Y[i] + Z[i]))
  end 

  fig = Figure()
  ax1 = Axis(fig[1, 1]); ax2 = Axis(fig[2, 1]); ax3 = Axis(fig[3, 1])
  lines!(ax1, xs, X)
  lines!(ax2, xs, Y)
  lines!(ax3, xs, Z, label = "Recovered")
  lines!(ax3, xs, N, label = "Total population")
  ax3.xlabel = "Time, years"
  ax1.ylabel = "Number susceptible"
  ax2.ylabel = "Number infectious"
  ax3.ylabel = "Numbers"
  fig[3, 2] = Legend(fig, ax3)
  linkxaxes!(ax1, ax2, ax3)
  hidexdecorations!(ax1; grid = false); hidexdecorations!(ax2; grid = false)

  return fig
end 

end # module MID_2_3
