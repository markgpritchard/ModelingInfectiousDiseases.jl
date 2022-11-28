
module MID_62

# SIR model with scaled additive noise (page 197)
  
using CairoMakie, DataFrames, DifferentialEquations, Distributions, Random
import Base: minimum

export Parameters62, sir62!, run_sir62, plot_sir62


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Structures that pass the model parameters to the relevant functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

struct Parameters62
    beta    :: Float64
    gamma   :: Float64
    mu      :: Float64
    nu      :: Float64
    xi      :: Float64
end

# The model runs the ODE solver for a series of small time intervals. For each time 
# interval a new vector of `Noise` parameters is calculated from Parameters62.xi . 
# This structure, with the calculated `Noise` parameter is then passed to the ODE 
# solver. 
mutable struct InputParameters62
    beta    :: Float64
    gamma   :: Float64
    mu      :: Float64
    nu      :: Float64
    p       :: Vector{<:Float64} # name of the Noise vector in this model 
end


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The function that is passed to the DifferentialEquations solver 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function sir62!(du, u, p, t) 
    # compartments 
    X, Y, Z = u
    N = X + Y + Z
    
    # ODEs - these equations are broken down into further functions to make the code 
    # easier to follow
    du[1] = birth(p.nu, N, p.p[1]) - infection(p.beta, X, Y, N, p.p[2]) - 
        death(p.mu, X, p.p[3])                                                  # dX
    du[2] = infection(p.beta, X, Y, N, p.p[2]) - recovery(p.gamma, Y, p.p[4]) - 
        death(p.mu, Y, p.p[5])                                                  # dY
    du[3] = recovery(p.gamma, Y, p.p[4]) - death(p.mu, Z, p.p[6])               # dZ
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function to run the whole model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function run_sir62(; N0 = 0, X0, Y0, Z0 = N0 - (X0 + Y0), beta, gamma, mu, nu = mu, xi, duration, kwargs...)
    u0 = [X0, Y0, Z0]
    p = Parameters62(beta, gamma, mu, nu, xi)
    return run_sir62(u0, p, duration; kwargs...)
end 

run_sir62(u0, p, duration; seed = nothing, kwargs...) = _run_sir62(u0, p, duration, seed; kwargs...)

function _run_sir62(u0, p, duration, seed::Real; kwargs...)
    Random.seed!(seed)
    return _run_sir62(u0, p, duration, nothing; kwargs...)
end 

_run_sir62(u0, p, duration, seed::Nothing; δt = 1, kwargs...) = __run_sir62(u0, p, duration, δt; kwargs...)

# It is much more convenient to ensure at this stage that δt will always be a Float64
__run_sir62(u0, p, duration, δt; kwargs...) = __run_sir62(u0, p, duration, Float64(δt); kwargs...)

function __run_sir62(u0, p, duration, δt::Float64; kwargs...)
    @assert minimum(u0) >= 0 "Model cannot run with negative starting values in `u0`. Model supplied u0 = $u0."
    @assert minimum(p) >= 0 "Model cannot run with negative parameters. Model supplied with p = $p." 
    @assert duration > 0 "Model needs duration > 0. Model supplied duration = $duration."
    @assert δt > 0 "Model needs δt > 0. Model supplied δt = $δt."
    @assert δt <= duration "Model needs δt <= duration. Model supplied δt = $δt and duration = $duration."

    # The model runs deterministically for a period δt, then a new `Noise` is calculated 
    # and the model runs for a further δt until duration is reached 
    u = u0
    parameters = InputParameters62(
            p.beta,
            p.gamma,
            p.mu,
            p.nu,
            zeros(6)    # the Noise parameter is added in the loop below before the ODE is run
        )
    τ0 = .0
    τ1 = δt 
    results = DataFrame(t = Float64[], X = Float64[], Y = Float64[], Z = Float64[])
    push!( results, Dict(:t => τ0, :X => u[1], :Y => u[2], :Z => u[3]) )
    while τ1 <= duration 
        tspan = ( τ0, τ1 )
        noise!(parameters, p, δt)   # add the noise parameter
        prob = ODEProblem(sir62!, u, tspan, parameters)
        sol = solve(prob; kwargs...)
        us = last(sol)
        u = [ max(1e-12, x) for x ∈ us ]
        τ0 = τ1
        τ1 += δt 
        push!( results, Dict(:t => τ0, :X => u[1], :Y => u[2], :Z => u[3]) )
    end

    return results
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Accessory functions for sir62! and run_sir62 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Define the Base function minimum for Parameters62

minimum(p::Parameters62) = min(p.beta, p.gamma, p.mu, p.nu, p.xi)


## Functions that calculate the appropriate `noise` values

noise(p, δt) = p.xi * rand(Normal(0, 1)) / sqrt(δt)
function noise!(parameters, p, δt) 
    for i ∈ 1:6
        parameters.p[i] = noise(p, δt)
    end 
end 

function addnoise(a, p) 
    # With large magnitudes of noise, negative values can be passed to the functions, 
    # such that a negative number of people might be born or might be infected. 
    # In turn, this can lead to negative numbers of individuals in compartments, 
    # even when tolerance is set very low. 

    # This function is designed 
    #   1.  to ensure that a negative value is never passed to the function 
    #   2.  to ensure that if a compartment ever has a negative value, this does 
    #       not lead to an attempt to square-root a negative value.
    if a < 0 
        return 0 
    elseif a + p * sqrt(a) < 0 
        return 0 
    else 
        return a + p * sqrt(a) 
    end 
end 

## Functions for events that occur in sir62!

birth(nu, N, p) = addnoise(nu * N, p)
infection(beta, X, Y, N, p) = addnoise(beta * X * Y / N, p)
recovery(gamma, Y, p) = addnoise(gamma * Y, p)
death(mu, A, p) = addnoise(mu * A, p)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plotting functions
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function plot_sir62(results)
    return plot_sir62(
        results, 
        "p6.2.jl: SIR model with random noise added to each parameter"
    )
end

function plot_sir62(results, noise::Real)
    return plot_sir62(
        results, 
        "p6.2.jl: SIR model with random noise added to each parameter\nNoise magnitude = $noise"
    )
end

plot_sir62(results, p::Parameters62) = plot_sir62(results, p.xi) 

function plot_sir62(results, label::String)
    fig = Figure()
    axs = [ Axis(fig[i, 1]) for i ∈ 1:3 ]
    for i ∈ 1:3
        lines!(axs[i], results.t ./ 365, results[:, i+1])
        i <= 2 && hidexdecorations!(axs[i]; grid = false, ticks = false)
    end 
    linkxaxes!(axs...)
    axs[3].xlabel = "Time, years"
    axs[1].ylabel = "Susceptible"
    axs[2].ylabel = "Infected"
    axs[3].ylabel = "Recovered"
    Label(fig[0, :], label; justification = :left)
    
    return fig
end 

end # module MID_62
