
module MID_27

using CairoMakie, DifferentialEquations

export plot_sir27, print_sir27, run_sir27, sir27!

"""
    sir27!(du, u, p, t) 

A compartmental susceptible--infectious--resistant model with a carrier state.
    
This is the ordinary differential equations function for programme 2.7 in 
    `Modeling Infectious Diseases in Humans and Animals`
"""
function sir27!(du, u, p, t)
    # compartments 
    S, I, C = u 
    # parameters 
    beta, gamma, Gamma, epsilon, mu, q = p # note that gamma and Gamma are different parameters

    # the ODEs
    du[1] = dS = mu - beta * S * (I + epsilon * C) - mu * S
    du[2] = dI = beta * S * (I + epsilon * C) - gamma * I - mu * I
    du[3] = dC = gamma * q * I - Gamma * C - mu * C
end 

"""
    run_sir27([; beta, gamma, I0, duration, saveat])

Run the model `sir27`

# Keyword arguments 

All keyword arguments are optional with default values supplied for each. 

* `beta`: The beta parameter in the model (infectiousness of infectives). Default is 0.2.
* `gamma`: The gamma (γ) parameter in the model (recovery rate). Default is 1 / 100.
* `Gamma`: The capital Gamma (Γ) parameter in the model (recovery rate of carriers). Default is 1 / 1000.
* `epsilon`: The epsilon parameter in the model (proportion reduction in transmission 
    from carriers compared with standard infectious compartment). Default is 0.1.
* `mu`: Per capita birth and death rate. Default is 1 / 50 years.
* `q`: The proportion of infected individuals who become carriers. Default is 0.4.
* `S0`: Proportion susceptible at time = 0. Default is 0.1.
* `I0`: Proportion infectious at time = 0. Default is 1e-4.
* `C0`: Proportion carriers at time = 0. Default is 1e-3.
* `duration`: How long the model will run (time units are interpretted as days). 
    Default is 60 years
* `saveat`: How frequently the model should save values. Default is 1 (day).
"""
function run_sir27(; beta = 0.2, gamma = 1 / 100, Gamma = 1 / 1000, epsilon = 0.1, mu = 1 / (50 * 365), q = .4, 
        S0 = .1, I0 = 1e-4, C0 = 1e-3, duration = 60 * 365, saveat = 1)
    
    @assert S0 >= 0 "Input S0 = $S0: cannot run with negative initial proportion susceptible"
    @assert I0 >= 0 "Input I0 = $I0: cannot run with negative initial proportion infectious"
    @assert C0 >= 0 "Input C0 = $C0: cannot run with negative initial proportion carriers"
    @assert beta >= 0 "Input beta = $beta: cannot run with negative parameter beta"
    @assert gamma >= 0 "Input gamma = $gamma: cannot run with negative parameter gamma (γ)"
    @assert Gamma >= 0 "Input Gamma = $Gamma: cannot run with negative parameter Gamma (Γ)"
    @assert epsilon >= 0 "Input epsilon = $epsilon: cannot run with negative parameter epsilon"
    @assert mu >= 0 "Input mu = $mu: cannot run with negative parameter mu"
    @assert q >= 0 "Input q = $q: cannot run with negative parameter sigqma"
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"
    @assert S0 + I0 + C0 <= 1 "Initial propotions susceptible ($S0), infectious ($I0) and carriers ($C0) cannot sum to greater than 1" 
    if beta / (gamma + mu) + beta * q * epsilon / (Gamma + mu) < 1 
        @info "R₀ < 1 ($( beta / (gamma + mu) + beta * q * epsilon / (Gamma + mu) ))" 
    end

    u0 = [S0, I0, C0]
    tspan = ( 0., Float64(duration) )
    p = [beta, gamma, Gamma, epsilon, mu, q]
    prob = ODEProblem(sir27!, u0, tspan, p)
    # set a low tolerance for the solver, otherwise the oscillations don't converge appropriately
    sol = solve(prob; saveat)

    return sol
end 

"""
    print_sir27([; kwargs...])

Print the saved values after running the model `sir27`. 

Keyword arguments are all optional. See `run_sir27` for details of the arguments and their default values.
"""
function print_sir27(; kwargs...)
    sol = run_sir27(; kwargs...)
    print_sir27(sol)
end 

function print_sir27(sol)
    for i ∈ eachindex(sol.u)
        println("t = $(sol.t[i]): $(sol.u[i])")
    end 
    return nothing 
end 

"""
    plot_sir27([; kwargs...])
    plot_sir27(sol)

Plot the results of running the model `sir27`. 

Can take optional keyword arguments, `run_sir27` (see that function for details 
    of the arguments and their default values), or the solution from an ODE model. 
    The advantage of allowing the plot function to run the ODE is that `saveat` 
    is selected to provide a smooth line on the plot.
"""
function plot_sir27(; kwargs...)
    sol = run_sir27(; kwargs...)
    return plot_sir27(sol)
end 

function plot_sir27(sol)
    xs = sol.t ./ 365 # to plot time in years
    S = Float64[]; I = Float64[]; C = Float64[]; R = Float64[]
    for i ∈ eachindex(sol.u)
        push!(S, sol.u[i][1])
        push!(I, sol.u[i][2])
        push!(C, sol.u[i][3])
        push!(R, 1 - (S[i] + I[i] + C[i]))
    end 

    fig = Figure()
    ax1 = Axis(fig[1, 1]); ax2 = Axis(fig[2, 1]); ax3 = Axis(fig[3, 1])
    lines!(ax1, xs, S)
    lines!(ax2, xs, I)
    lines!(ax3, xs, C)
    Label(
        fig[0, :], 
        "p2.7.jl: Susceptible--infectious--resistant model with a carrier state"
    )
    ax3.xlabel = "Time, years"
    ax1.ylabel = "Proportion susceptible"
    ax2.ylabel = "Proportion infectious"
    ax3.ylabel = "Proportion carriers"
    linkxaxes!(ax1, ax2, ax3)
    hidexdecorations!(ax1; grid = false, ticks = false)
    hidexdecorations!(ax2; grid = false, ticks = false)
    resize_to_layout!(fig)
    
    return fig
end 

end # module MID_27
