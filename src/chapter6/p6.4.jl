
module MID_64
  
using CairoMakie, DataFrames, Random

export sir64, u0_sir64, run_sir64, plot_sir64

const changematrix =    # how contents of each compartment change for each possible event     
    [   -1  1   0       # infection 
        0   -1  1       # recovery
        1   0   0       # birth 
        -1  0   0       # susceptible death 
        0   -1  0       # infectious death 
        0   0   -1  ]   # recovered death

"""
    u0_sir64(N0, p) 

Produce a set of values for u0 based on the model parameters and the initial population 
size.

## Example
```
julia> p = [1., .1, 5e-4]
3-element Vector{Float64}:
 1.0
 0.1
 0.0005

julia> u0 = u0_sir64(5000, p)
3-element Vector{Int64}:
  500
   25
 4475
```
"""
function u0_sir64(N0::Int, p) 
    β, γ, μ = p
    X0 = round(Int, γ * N0 / β, RoundDown)
    Y0 = round(Int, μ * N0 / γ, RoundUp) 
    Z0 = N0 - ( X0 + Y0 ) 
    return [X0, Y0, Z0]
end 

function sir64(u, p, t, N0) 
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

    r1 = rand(); r2 = rand() 

    # what to do? 
    # to choose which event will occur we need a cumulative sum of their rates divided 
    # by the total sum of rates
    cumsumratesratio = cumsum(rates) / sum(rates)
    # then we find the first value >= r1 
    i = searchsortedfirst(cumsumratesratio, r1)
    u += changematrix[i, :]

    # when to do it?
    timestep = -log(r2) / sum(rates)
    t += timestep

    return t, u
end 

"""
    run_sir64(N0::Int, p, duration[; seed, pop])
    run_sir64(u0::Vector{<:Int}, p, duration[; seed, pop])

Run the model `sis64`.

This is a susceptible--infectious--recovered model. It calculates rates of infection, 
recovery, births and deaths. It has random step intervals, proportional to these 
rates. The model keeps running until one time step has been taken after `duration` 
has been reached.

The number of infectious and susceptible individuals will always be integer, and 
the model has the potential for all infectious individuals to recover and the epidemic 
ends.

## Model parameters 
Parameters can be entered as a vector in this order
* `β`: infection parameter
* `γ`: recovery rate
* `μ`: birth and death rate 

## Function arguments
* `u0`: The starting conditions for the model, a vector of 3 integers (`X0`, `Y0`, `Z0`),
    or a single integer that will be interpretted as `N0`.
* `p`: Parameters for the model, expected as a vector
* `duration`: The time that the model should run for
### Optional keyword arguments
* `seed`: Seed for the random number generator. Default is not to supply a seed.
* `pop`: Whether the final value of the DataFrame should be removed. Default is `true`

## Example 
```
julia> p = [1., .1, 5e-4]
3-element Vector{Float64}:
 1.0
 0.1
 0.0005

julia> u0 = u0_sir64(50, p)
3-element Vector{Int64}:
  5
  1
 44

julia> run_sir64(u0, p, 10; seed = 64)
4×4 DataFrame
 Row │ t        X      Y      Z     
     │ Float64  Int64  Int64  Int64 
─────┼──────────────────────────────
   1 │ 0.0          5      1     44
   2 │ 4.3724       4      2     44
   3 │ 6.44244      4      1     45
   4 │ 7.86218      3      2     45
```
"""
run_sir64(u0, p, duration; seed = nothing, kwargs...) = _run_sir64(u0, p, duration, seed; kwargs...)

function _run_sir64(u0, p, duration, seed::Int; kwargs...)
    Random.seed!(seed)
    return _run_sir64(u0, p, duration, nothing; kwargs...)
end 

_run_sir64(N0::Int, p, duration, seed::Nothing; kwargs...) = 
    _run_sir64(u0_sir64(N0, p), p, duration, seed; kwargs...)

function _run_sir64(u0::Vector{<:Int}, p, duration, seed::Nothing; pop = true)
    @assert minimum(u0) >= 0 "Model cannot run with negative starting values in `u0`. Model supplied u0 = $u0."
    @assert minimum(p) >= 0 "Model cannot run with negative parameters. Running with p = $p."
    @assert duration > 0 "Model needs duration > 0. Model supplied duration = $duration."

    N0 = sum(u0)
    t = .0 
    u = u0
    results = DataFrame(
        t = Float64[], X = typeof(u0[1])[], Y = typeof(u0[2])[], Z = typeof(u0[3])[]
    )
    push!( results, Dict(:t => t, :X => u[1], :Y => u[2], :Z => u[3]) )
    while t <= duration 
        t, u = sir64(u, p, t, N0) 
        push!( results, Dict(:t => t, :X => u[1], :Y => u[2], :Z => u[3]) )
    end

    if pop pop!(results) end
    return results
end 

"""
    plot_sir64(results[, population])
    plot_sir64(results, label::String)

Plot the `results` DataFrame output from the function `run_sir64` 

A `label` term can be added which will be printed at the top of the figure. If a 
`population` term is included, the size of the population is printed on the plot.
"""
plot_sir64(results) = plot_sir64(results, "p6.4.jl: SIR model with demographic stochasticity")

function plot_sir64(results, population::Real)
    return plot_sir64(
        results, 
        "p6.4.jl: SIR model with demographic stochasticity\nPopulation = $population"
    )
end 

function plot_sir64(results, label::String)
    fig = Figure()
    axs = [ Axis(fig[i, 1]) for i ∈ 1:3 ]
    for i ∈ 1:3
        stairs!(axs[i], results.t ./ 365, results[:, i+1]; step = :post)
        if i <= 2 hidexdecorations!(axs[i]; ticks = false) end
    end 
    linkxaxes!(axs...)
    axs[3].xlabel = "Time, years"
    axs[1].ylabel = "Susceptible"
    axs[2].ylabel = "Infected"
    axs[3].ylabel = "Recovered"
    Label(fig[0, :], label)
    
    return fig
end 

end # module MID_64
