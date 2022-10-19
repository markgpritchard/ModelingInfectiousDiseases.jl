
module MID_31
  
using CairoMakie, DifferentialEquations

export Parameters31, sir31!, run_sir31, print_sir31, plot_sir31

struct Parameters31 # beta will be a matrix but gamma is a float, so group them in a structure
    beta    :: Matrix{<:Real} 
    gamma   :: Real
end

"""
    sir31!(du, u, p, t) 

A compartmental susceptible--infectious--resistant model with a high-risk and low-risk
    populations. This is the ordinary differential equations function for programme 3.1 in 
    `Modeling Infectious Diseases in Humans and Animals`


"""
function sir31!(du, u, p, t) 
    # compartments 
    Sh, Ih, Sℓ, Iℓ = u

    # parameters 
    (βhh, βhl, βlh, βll) = p.beta 
    γ = p.gamma
    
    # the ODEs
    du[1] = -(βhh * Ih + βhl * Iℓ) * Sh + γ * Ih    # dSh
    du[2] = (βhh * Ih + βhl * Iℓ) * Sh - γ * Ih     # dIh
    du[3] = -(βlh * Ih + βll * Iℓ) * Sℓ + γ * Iℓ    # dSℓ
    du[4] = (βlh * Ih + βll * Iℓ) * Sℓ - γ * Iℓ     # dIℓ
end 

"""
    run_sir31
"""
function run_sir31(; beta = [10 .1; .1 1], gamma = 1, nh = .2, Ih = 1e-5, Il = 1e-3, duration = 15, saveat = 1)
    @assert minimum(beta) >= 0 "Input beta cannot include any negative transmission rates"
    @assert Ih >= 0 "Input Ih = $Ih: cannot run with negative initial proportion high-risk infectious"
    @assert Il >= 0 "Input Il = $Il: cannot run with negative initial proportion low-risk infectious"
    @assert nh >= 0 "Input nh = $nh: cannot run with negative initial proportion in high-risk group"
    @assert nh <= 1 "Input nh = $nh: cannot run with >1 proportion in high-risk group"
    @assert Ih <= nh "Input Ih = $Ih: cannot run with greater high-risk infectious proportion than total size of high-risk group" 
    @assert Il <= 1 - nh "Input Il = $Il: cannot run with greater low-risk infectious proportion than total size of low-risk group" 
    @assert gamma >= 0 "Input gamma = $gamma: cannot run with negative parameter gamma" 
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"
    # doesn't seem to be included in the other version of the model but seems important 
    @assert minimum(beta) >= 0 "Smallest input beta = $(minimum(beta)): cannot run with negative parameter beta" 

    nl = 1 - nh
    Sh = nh - Ih
    Sl = nl - Il

    u0 = [Sh, Ih, Sl, Il]
    tspan = ( 0., Float64(duration) )
    p = Parameters31(beta, gamma)
    prob = ODEProblem(sir31!, u0, tspan, p)
    sol = solve(prob; saveat)

    return sol
end 

"""
    print_sir31([; kwargs...])

Print the saved values after running the model `sir31`. 

Keyword arguments are all optional. See `run_sir31` for details of the arguments and their default values.
"""
function print_sir31(; kwargs...)
    sol = run_sir31(; kwargs...)
    print_sir31(sol)
end 

function print_sir31(sol)
    for i ∈ eachindex(sol.u)
        println("t = $(sol.t[i]): $(sol.u[i])")
    end 
    return nothing 
end 

"""
    plot_sir31([; kwargs...])
    plot_sir31(sol)

Plot the results of running the model `sir31`. 

Can take optional keyword arguments, `run_sir31` (see that function for details 
    of the arguments and their default values), or the solution from an ODE model. 
    The advantage of allowing the plot function to run the ODE is that `saveat` 
    is selected to provide a smooth line on the plot.

# Examples
    julia> plot_sir31()

    julia> plot_sir31(beta = .8, gamma = .6, duration = 100) 
"""
function plot_sir31(; duration = 15, kwargs...)
    saveat = duration / 500 # to give a smoother line for plotting
    sol = run_sir31(; duration, saveat, kwargs...)
    return plot_sir31(sol)
end 

function plot_sir31(sol)
    xs = sol.t 
    Ih = Float64[]; Iℓ = Float64[];
    for i ∈ eachindex(sol.u)
        push!(Ih, sol.u[i][2])
        push!(Iℓ, sol.u[i][4])
    end 

    fig = Figure()
    ax1 = Axis(fig[1, 1])
    lines!(ax1, xs, Ih, label = "High risk")
    lines!(ax1, xs, Iℓ, label = "Low risk")
    Label(
        fig[0, :], 
        "p3.1.jl: Susceptible--infectious--resistant model with a constant population"
    )
    ax1.ylabel = "Fraction infectious\n(linear scale)"
    ax2 = Axis(fig[2, 1], yscale = log10)
    lines!(ax2, xs, Ih, label = "High risk")
    lines!(ax2, xs, Iℓ, label = "Low risk")
    ax2.xlabel = "Time, days"
    ax2.ylabel = "Fraction infectious\n(log scale)"
    fig[1:2, 2] = Legend(fig, ax1)
    resize_to_layout!(fig)

    return fig
end 

end # module MID_31
