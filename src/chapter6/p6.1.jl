
module MID_61
  
using CairoMakie, DataFrames, DifferentialEquations, Distributions, Random
import Base: minimum

export Parameters61, sir61!, run_sir61, plot_sir61, plot_sir61!


## Structures that pass the model parameters to the relevant functions

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

# Define the Base function minimum for Parameters61
minimum(p::Parameters61) = min(p.beta, p.gamma, p.mu, p.nu, p.xi)


## The function that is passed to the DifferentialEquations solver 

function sir61!(du, u, p, t) 
    # compartments 
    X, Y, Z = u
    N = X + Y + Z
    
    # ODEs
    du[1] = p.nu * N - ( p.beta * X * Y / N + p.Noise ) - p.mu * X  # dX
    du[2] = p.beta * X * Y / N + p.Noise - ( p.gamma + p.mu ) * Y   # dY
    du[3] = p.gamma * Y - p.mu * Z                                  # dZ
end 


## Functions that calculate the appropriate `noise` value

noise(p, δt) = p.xi * rand(Normal(0, 1)) / sqrt(δt)
noise!(parameters, p, δt) = parameters.Noise = noise(p, δt)


## Function to run the whole model

run_sir61(u0, p, duration; δt = 1, seed = nothing) = _run_sir61(u0, p, duration, seed; δt)

function _run_sir61(u0, p, duration, seed::Real; δt)
    Random.seed!(seed)
    return _run_sir61(u0, p, duration, nothing; δt)
end 

_run_sir61(u0, p, duration, seed::Nothing; δt) = __run_sir61(u0, p, duration, δt)

# It is much more convenient to ensure at this stage that δt will always be a Float64
__run_sir61(u0, p, duration, δt) = __run_sir61(u0, p, duration, Float64(δt))

function __run_sir61(u0, p, duration, δt::Float64)
    @assert minimum(u0) >= 0 "Model cannot run with negative starting values in `u0`. Model supplied u0 = $u0."
    if minimum(p) < 0 @warn "Model may be unreliable with negative parameters. Running with p = $p." end
    @assert duration > 0 "Model needs duration > 0. Model supplied duration = $duration."
    @assert δt > 0 "Model needs δt > 0. Model supplied δt = $δt."
    @assert δt <= duration "Model needs δt <= duration. Model supplied δt = $δt and duration = $duration."

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


## Plotting functions

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
        if i <= 2 hidexdecorations!(axs[i]; ticks = false) end
    end 
    linkxaxes!(axs...)
    axs[3].xlabel = "Time, years"
    axs[1].ylabel = "Susceptible"
    axs[2].ylabel = "Infected"
    axs[3].ylabel = "Recovered"
    Label(gl[0, :], label; justification = :left)
end 


## Help text 

"""
    run_sir61(u0, p, duration[; δt, seed])

Run the model `sir61!`.

`sir61!` is an ordinary differential equations (ODE) model, but this function is 
intended to introduce stochastic noise. It does this by introducing a stochastic 
parameter, which is inversely proportional to the square root of `δt`. The model 
runs for a duration `δt` before calculating a new, independent, noise parameter. 
This continues until `duration` has been reached.

## Model parameters 
Parameters are expected in a `Parameters61` type
* `beta`: infection parameter
* `gamma`: recovery rate
* `mu`: birth rate 
* `nu`: death rate 
* `xi`: magnitude of noise to be added

## Function arguments
* `u0`: The starting conditions for the model, a vector of 3 values (`X0`, `Y0`, `Z0`)
* `p`: Parameters for the model, expected in a `Parameters61` type
* `duration`: The time that the model should run for
### Optional keyword arguments
* `δt`: How often the random noise parameter should update. Default value is 1.
* `seed`: Seed for the random number generator. Default is not to supply a seed.

## Example 
```
julia> u0 = [1e5, 500, 1e6]
3-element Vector{Float64}:
 100000.0
    500.0
      1.0e6

julia> p = Parameters61(1., .1, 1 / (50 * 365), 1 / (50 * 365), 10.)
Parameters61(1.0, 0.1, 5.479452054794521e-5, 5.479452054794521e-5, 10.0)

julia> run_sir61(u0, p, 5; seed = 61)
6×4 DataFrame
 Row │ t        X               Y        Z
     │ Float64  Float64         Float64  Float64
─────┼──────────────────────────────────────────────────
   1 │     0.0  100000.0        500.0         1.0e6
   2 │     1.0       1.00008e5  497.192       9.99995e5
   3 │     2.0       1.00006e5  503.268       9.9999e5
   4 │     3.0       1.00024e5  491.14        9.99985e5
   5 │     4.0       1.00041e5  480.24   999979.0
   6 │     5.0       1.00056e5  471.756       9.99972e5
```
"""
function run_sir61() end

"""
    plot_sir61(results[, noise])
    plot_sir61(results, label::String)

Plot the `results` DataFrame output from the function `run_sir61`.
    
A `label` term can be added which will be printed at the top of the figure. If a 
`noise` term is included, the magnitude of the noise is printed on the plot. `noise` 
can be a value or a `Parameters61` structure.
"""
function plot_sir61() end

end # module MID_6_1
