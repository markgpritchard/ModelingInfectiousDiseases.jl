
module MID_61

# SIR model with Constant additive noise (page 194)
  
using CairoMakie, DataFrames, DifferentialEquations, Distributions, Random
import Base: minimum

export Parameters61, sir61!, run_sir61, plot_sir61, plot_sir61!


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Structures that pass the model parameters to the relevant functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

struct Parameters61 
    beta    :: Float64
    gamma   :: Float64
    mu      :: Float64
    nu      :: Float64
    xi      :: Float64
end

# The model runs the ODE solver for a series of small time intervals. For each time 
# interval a new `Noise` parameter is calculated from Parameters61.xi . This structure, 
# with the calculated `Noise` parameter is then passed to the ODE solver. 
mutable struct InputParameters61
    beta    :: Float64
    gamma   :: Float64
    mu      :: Float64
    nu      :: Float64
    Noise   :: Float64
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The function that is passed to the DifferentialEquations solver 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function sir61!(du, u, p, t) 
    # compartments 
    X, Y, Z = u
    N = X + Y + Z
    
    # ODEs
    du[1] = p.nu * N - ( p.beta * X * Y / N + p.Noise ) - p.mu * X  # dX
    du[2] = p.beta * X * Y / N + p.Noise - ( p.gamma + p.mu ) * Y   # dY
    du[3] = p.gamma * Y - p.mu * Z                                  # dZ
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to run the whole model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

run_sir61(u0, p, duration; δt = 1, seed = nothing) = _run_sir61(u0, p, duration, seed; δt)

function run_sir61(; N0 = 0, X0, Y0, Z0 = N0 - (X0 + Y0), beta, gamma, mu, nu = mu, xi, duration, kwargs...)
    u0 = [X0, Y0, Z0]
    p = Parameters61(beta, gamma, mu, nu, xi)
    return run_sir61(u0, p, duration; kwargs...)
end 

function _run_sir61(u0, p, duration, seed::Real; δt)
    Random.seed!(seed)
    return _run_sir61(u0, p, duration, nothing; δt)
end 

_run_sir61(u0, p, duration, seed::Nothing; δt) = __run_sir61(u0, p, duration, δt)

# It is much more convenient to ensure at this stage that δt will always be a Float64
__run_sir61(u0, p, duration, δt) = __run_sir61(u0, p, duration, Float64(δt))

function __run_sir61(u0, p, duration, δt::Float64)
    @assert minimum(u0) >= 0 "Input u0 = $u0: Cannot run model with negative starting values in any compartment"
    @assert minimum(p) >= 0 "Input p = $p: Cannot run with negative values for any parameters"
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"
    @assert δt > 0 "Input δt = $δt: model needs δt > 0"
    @assert δt <= duration "Input δt = $δt, duration = $duration: Model needs δt ≤ duration"

    # The model runs deterministically for a period δt, then a new `Noise` is calculated 
    # and the model runs for a further δt until duration is reached 
    u = u0
    parameters = InputParameters61(
        p.beta,
        p.gamma,
        p.mu,
        p.nu,
        .0      # the Noise parameter is added in the loop below before the ODE is run
    )
    τ0 = .0
    τ1 = δt 
    results = DataFrame(t = Float64[], X = Float64[], Y = Float64[], Z = Float64[])
    push!( results, Dict(:t => τ0, :X => u[1], :Y => u[2], :Z => u[3]) )
    while τ1 <= duration 
        tspan = ( τ0, τ1 )
        noise!(parameters, p, δt)   # add the noise parameter
        prob = ODEProblem(sir61!, u, tspan, parameters)
        sol = solve(prob)
        u = last(sol)
        τ0 = τ1
        τ1 += δt 
        push!( results, Dict(:t => τ0, :X => u[1], :Y => u[2], :Z => u[3]) )
    end

    return results
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Accessory functions for sir61! and run_sir61 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the Base function minimum for Parameters61

minimum(p::Parameters61) = min(p.beta, p.gamma, p.mu, p.nu, p.xi)

## Functions that calculate the appropriate `noise` value

noise(p, δt) = p.xi * rand(Normal(0, 1)) / sqrt(δt)
noise!(parameters, p, δt) = parameters.Noise = noise(p, δt)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plotting functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function plot_sir61(results, noise::Real)
    return plot_sir61(
        results, 
        "p6.1.jl: SIR model with random noise added to the transmission term\nNoise magnitude = $noise"
    )
end

plot_sir61(results, p::Parameters61) = plot_sir61(results, p.xi) 

function plot_sir61(results, label = "p6.1.jl: SIR model with random noise added to the transmission term")
    fig = Figure()
    plot_sir61!(fig, results, label)
    resize_to_layout!(fig)
    return fig 
end 

function plot_sir61!(fig::Figure, results, label)
    gl = GridLayout(fig[1, 1])
    plot_sir61!(gl, results, label)
end 

function plot_sir61!(gl::GridLayout, results, label)
    axs = [ Axis(gl[i, 1]) for i ∈ 1:3 ]
    for i ∈ 1:3
        lines!(axs[i], results.t ./ 365, results[:, i+1])
        i <= 2 && hidexdecorations!(axs[i]; grid = false, ticks = false)
    end 
    linkxaxes!(axs...)
    axs[3].xlabel = "Time, years"
    axs[1].ylabel = "Susceptible"
    axs[2].ylabel = "Infected"
    axs[3].ylabel = "Recovered"
    Label(gl[0, :], label; justification = :left)
end 

end # module MID_6_1
