
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
    S, I, R = u
    beta, gamma = p 

    dS = -beta * S * I 
    dI = beta * S * I - gamma * I 
    dR = gamma * I 
    return @SVector [dS, dI, dR]
end 

"""
    run_sir_21([; beta, gamma, S0, I0, R_at_time0, duration, saveat])

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
  
    @assert S0 >= 0 "Input S0 = $S0: cannot run with negative value of S0"
    @assert I0 >= 0 "Input I0 = $I0: cannot run with negative value of I0"
    @assert R_at_time0 >= 0 "Input R_at_time0 = $R_at_time0: cannot run with negative value of R_at_time0"
    @assert beta >= 0 "Input beta = $beta: cannot run with negative value of beta"
    @assert gamma >= 0 "Input gamma = $gamma: cannot run with negative value of gamma"
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero value of duration"
    if S0 + I0 > 1 @warn "Values of S0 ($S0) and I0 ($I0) sum to more than 1" end 
    if beta < gamma @info "R₀ < 1 ($(beta / gamma))" end

    u0::SVector{3, Float64} = @SVector [S0, I0, R_at_time0]
    tspan::Tuple{Float64, Float64} = ( 0., Float64(duration) )
    p::Vector{Float64} = [beta, gamma]
    prob = ODEProblem(sir_21, u0, tspan, p)
    sol = solve(prob; saveat)

    return sol
end 

"""
    print_sir_21([; kwargs...])

Print the saved values after running the model `sir_21`. Keyword arguments are those of `run_sir_21`
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

Plot the results of running the model `sir_21`. Keyword arguments are those of 
    `run_sir_21`. It is recommended not to set `saveat`, as the plotting function 
    automatically uses a value that plots a smooth line.
"""
function plot_sir_21(; duration = 70, kwargs...)
    saveat = duration / 500 # to give a smoother line for plotting
    sol = run_sir_21(; duration, saveat, kwargs...)
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
