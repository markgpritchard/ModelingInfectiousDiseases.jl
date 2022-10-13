
module MID_65
  
using CairoMakie, DataFrames, Distributions, Random

export sir65, u0_sir65, run_sir65, plot_sir65

const changematrix =    # how contents of each compartment change for each possible event     
    [   -1  1   0       # infection 
        0   -1  1       # recovery
        1   0   0       # birth 
        -1  0   0       # susceptible death 
        0   -1  0       # infectious death 
        0   0   -1  ]   # recovered death

"""
    u0_sir65(N0, p) 

Produce a set of values for `u0` based on the model parameters and the initial population 
size.

## Example
```
julia> p = [1., .1, 5e-4]
3-element Vector{Float64}:
 1.0
 0.1
 0.0005

julia> u0 = u0_sir65(5000, p)
3-element Vector{Int64}:
  500
   25
 4475
```
"""
function u0_sir65(N0::Int, p) 
    β, γ, μ = p
    X0 = round(Int, γ * N0 / β, RoundDown)
    Y0 = round(Int, μ * N0 / γ, RoundUp) 
    Z0 = N0 - ( X0 + Y0 ) 
    return [X0, Y0, Z0]
end 

function sir65(u, p, t, δt, N0) 
    # The model does not enforce a constant population. However, an approximately 
    # constant population is desired. If we calculated N = X + Y + Z then the birth 
    # rate would be proportional to the current population and could lead to a positive 
    # feedback. This function therefore takes an input of N0 to use in calculating 
    # birth rate.
    X, Y, Z = u 
    β, γ, μ = p 
    N = X + Y + Z

    rates = [           # rates at with each possible event occur
        β * X * Y / N,  # infection
        γ * Y,          # recovery 
        μ * N0,         # birth 
        μ * X,          # susceptible death 
        μ * Y,          # infectious death 
        μ * Z           # recovered death
    ]

    events = [ rand(Poisson(rate * δt)) for (i, rate) ∈ enumerate(rates) ]

    for (i, event) ∈ enumerate(events) 
        u += event * changematrix[i, :] 
        if minimum(u) < 0 # then this cycle has made more events than it could have 
            reduction = -minimum(u)                         # how many events were too many 
            u += -event * changematrix[i, :]                # undo what you just did 
            u += ( event - reduction ) * changematrix[i, :] # re-do it with the reduced number of events
        end 
    end 

    t += δt
    return t, u
end 

"""
    run_sir65(N0::Int, p, duration[; δt, seed])
    run_sir65(u0::Vector{<:Int}, p, duration[; δt, seed])

Run the model `sir65`.

This is a susceptible--infectious--recovered model. It calculates rates of infection, 
recovery, births and deaths. It has defined step intervals and estimates the number 
of each event that will occur during each step interval from a Poisson distribution.

## Model parameters 
Parameters can be entered as a vector in this order
* `β`: infection parameter
* `γ`: recovery rate
* `μ`: birth and death rate 

## Function arguments
* `u0`: The starting conditions for the model, a vector of 3 integers (`X0`, `Y0`, `Z0`), 
    or a single integer that will be interpretted as `N0`.
* `p`: Parameters for the model, expected as a vector.
* `duration`: The time that the model should run for
### Optional keyword arguments
* `δt`: Time interval between calculations. Default is 1.
* `seed`: Seed for the random number generator. Default is not to supply a seed.

## Example 
```
julia> u0 = [50, 50, 400]
3-element Vector{Int64}:
  50
  50
 400

julia> p = [1., .1, 5e-4]
3-element Vector{Float64}:
 1.0
 0.1
 0.0005

julia> run_sir65(u0, p, 5; seed = 65)
6×4 DataFrame
 Row │ t      X      Y      Z     
     │ Int64  Int64  Int64  Int64 
─────┼────────────────────────────
   1 │     0     50     50    400
   2 │     1     46     52    401
   3 │     2     42     47    410
   4 │     3     37     48    415
   5 │     4     36     43    423
   6 │     5     35     43    424
```
"""
run_sir65(u0, p, duration; δt = 1, seed = nothing) = _run_sir65(u0, p, duration, seed; δt)

function _run_sir65(u0, p, duration, seed::Int; δt)
    Random.seed!(seed)
    return _run_sir65(u0, p, duration, nothing; δt)
end 

_run_sir65(N0::Int, p, duration, seed::Nothing; δt) = 
    _run_sir65(u0_sir65(N0, p), p, duration, seed; δt)

function _run_sir65(u0::Vector{<:Int}, p, duration, seed::Nothing; δt)
    @assert minimum(u0) >= 0 "Model cannot run with negative starting values in `u0`. Model supplied u0 = $u0."
    @assert minimum(p) >= 0 "Model cannot run with negative parameters. Running with p = $p."
    @assert duration > 0 "Model needs duration > 0. Model supplied duration = $duration."
    @assert δt <= duration "Cannot run model with δt > duration"

    N0 = sum(u0)
    t = zero(typeof(δt)) 
    u = u0
    results = DataFrame(
        t = typeof(t)[], X = typeof(u0[1])[], Y = typeof(u0[2])[], Z = typeof(u0[3])[]
    )
    push!( results, Dict(:t => t, :X => u[1], :Y => u[2], :Z => u[3]) )
    while t < duration 
        t, u = sir65(u, p, t, δt, N0) 
        push!( results, Dict(:t => t, :X => u[1], :Y => u[2], :Z => u[3]) )
    end

    return results
end 


"""
    plot_sir65(results)

Plot the `results` DataFrame output from the function `run_sir65` 
"""
function plot_sir65(results)
    fig = Figure()
    axs = [ Axis(fig[i, 1]) for i ∈ 1:3 ]
    for i ∈ 1:3
        lines!(axs[i], results.t ./ 365, results[:, i+1]) # lines! instead of stairs! 
            # because we're no longer trying to simulate the instance events occurred
        if i <= 2 hidexdecorations!(axs[i]; ticks = false) end
    end 
    linkxaxes!(axs...)
    axs[3].xlabel = "Time, years"
    axs[1].ylabel = "Susceptible"
    axs[2].ylabel = "Infected"
    axs[3].ylabel = "Recovered"
    
    return fig
end 

end # module MID_65
