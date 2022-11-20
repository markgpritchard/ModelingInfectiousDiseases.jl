
module MID_63
  
using CairoMakie, DataFrames, Random

export sis63, run_sis63, plot_sis63

function sis63(u, p, t) 
    # This function does not use DifferentialEquations.jl, but I have kept the format 
    # similar for consistency 
    X, Y = u 
    β, γ = p 
    N = X + Y

    infectionrate = β * X * Y / N
    recoveryrate = γ * Y

    r1 = rand(); r2 = rand() 

    # what to do? 
    if r1 < infectionrate / (infectionrate + recoveryrate) 
        X += -1; Y += 1     # infection 
    else 
        X += 1; Y += -1     # recovery 
    end

    # when to do it?
    timestep = -log(r2) / (infectionrate + recoveryrate)
    t += timestep

    return t, [X, Y]
end 

"""
    run_sis63(u0::Vector{<:Int}, p, duration[; seed, pop])

Run the model `sis63`.

This is a susceptible--infectious--susceptible model. It calculates rates of infection 
and recovery. The model has random step intervals, proportional to these rates, 
at which an event takes place. The model keeps running until one time step has been 
taken after `duration` has been reached.

The number of infectious and susceptible individuals will always be integer, and 
the model has the potential for all infectious individuals to recover and the epidemic 
ends. In this case, the final time step value will be `Inf` and the final value 
of `Y` will be `-1`. Therefore, by default the final value from the DataFrame is 
removed. To avoid this, set `pop = false`

## Model parameters 
Parameters can be entered as a vector in this order
* `β`: infection parameter
* `γ`: recovery rate

## Function arguments
* `u0`: The starting conditions for the model, a vector of 2 integers  (`X0`, `Y0`)
* `p`: Parameters for the model, expected as a vector
* `duration`: The time that the model should run for
### Optional keyword arguments
* `seed`: Seed for the random number generator. Default is not to supply a seed.
* `pop`: Whether the final value of the DataFrame should be removed. Default is `true`

## Example 
```
julia> u0 = [30, 70]
2-element Vector{Int64}:
 30
 70

julia> p = [.03, .01]
2-element Vector{Float64}:
 0.03
 0.01

julia> run_sis63(u0, p, 3; seed = 63)
5×3 DataFrame
 Row │ t         X      Y     
     │ Float64   Int64  Int64 
─────┼────────────────────────
   1 │ 0.0          30     70
   2 │ 0.841074     29     71
   3 │ 1.03188      30     70
   4 │ 2.31459      29     71
   5 │ 2.8178       28     72
```
"""
run_sis63(u0::Vector{<:Int}, p, duration; seed = nothing, kwargs...) =
    _run_sis63(u0, p, duration, seed; kwargs...)

function _run_sis63(u0, p, duration, seed::Int; kwargs...)
    Random.seed!(seed)
    return _run_sis63(u0, p, duration, nothing; kwargs...)
end 

function _run_sis63(u0, p, duration, seed::Nothing; pop = true)
    @assert minimum(u0) >= 0 "Model cannot run with negative starting values in `u0`. Model supplied u0 = $u0."
    @assert minimum(p) >= 0 "Model cannot run with negative parameters. Running with p = $p."
    @assert duration > 0 "Model needs duration > 0. Model supplied duration = $duration."

    t = .0 
    u = u0
    results = DataFrame(t = Float64[], X = typeof(u0[1])[], Y = typeof(u0[2])[])
    push!( results, Dict(:t => t, :X => u[1], :Y => u[2]) )
    while t <= duration 
        t, u = sis63(u, p, t) 
        push!( results, Dict(:t => t, :X => u[1], :Y => u[2]) )
    end

    if pop pop!(results) end
    return results
end 

"""
    plot_sis63(results[, label])

Plot the `results` DataFrame output from the function `run_sis63` 

A `label` term can be added which will be printed at the top of the figure.
"""
plot_sis63(results) = plot_sis63(results, "p6.3.jl: SIS model with demographic stochasticity")

function plot_sis63(results, label::String)
    fig = Figure()
    axs = [ Axis(fig[i, 1]) for i ∈ 1:2 ]
    for i ∈ 1:2
        stairs!(axs[i], results.t ./ 365, results[:, i+1]; step = :post)
    end 
    linkxaxes!(axs...)
    hidexdecorations!(axs[1]; grid = false, ticks = false)
    axs[2].xlabel = "Time, years"
    axs[1].ylabel = "Susceptible"
    axs[2].ylabel = "Infected"
    Label(fig[0, :], label; justification = :left)
    
    return fig
end 

end # module MID_63
