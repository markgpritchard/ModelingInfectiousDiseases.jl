
module MID_63
  
using CairoMakie, DataFrames, Random

export sis63!, run_sis63, plot_sis63

function sis63!(u, p, t) 
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
    run_sir61(u0, p, duration[; δt, seed])

Run the model `sir61!`.

`sir61!` is an ordinary differential equations (ODE) model, but this function is 
intended to introduce stochastic noise. It does this by introducing a stochastic 
parameter, which is inversely proportional to the square root of `δt`. The model 
runs for a duration `δt` before calculating a new, independent, noise parameter. 
    This continues until `duration` has been reached.

## Parameters 
* `u0`: The starting conditions for the model, a vector of 3 values.
* `p`: Parameters for the model, including a term `xi` that represents the magnitude 
    of the random noise. Should be supplied as a `Parameters61` structure.
* `duration`: The time that the model should run for
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
run_sis63(u0, p, duration; seed = nothing) = _run_sis63(u0, p, duration, seed)

function _run_sis63(u0, p, duration, seed::Int)
    Random.seed!(seed)
    return _run_sis63(u0, p, duration, nothing)
end 

function _run_sis63(u0, p, duration, seed::Nothing)
    @assert minimum(u0) >= 0 "Model cannot run with negative starting values in `u0`. Model supplied u0 = $u0."
    if minimum(p) < 0 @warn "Model may be unreliable with negative parameters. Running with p = $p." end
    @assert duration > 0 "Model needs duration > 0. Model supplied duration = $duration."

    t = .0 
    u = u0
    results = DataFrame(t = Float64[], X = typeof(u0[1])[], Y = typeof(u0[2])[])
    push!( results, Dict(:t => t, :X => u[1], :Y => u[2]) )
    while t <= duration 
        t, u = sis63!(u, p, t) 
        push!( results, Dict(:t => t, :X => u[1], :Y => u[2]) )
    end

    return results
end 


"""
    plot_sir61(results)

Plot the `results` DataFrame output from the function `run_sir61` 
"""
function plot_sis63(results)
    fig = Figure()
    axs = [ Axis(fig[i, 1]) for i ∈ 1:2 ]
    for i ∈ 1:2
        stairs!(axs[i], results.t ./ 365, results[:, i+1]; step = :post)
        if i == 1 hidexdecorations!(axs[i]; ticks = false) end
    end 
    linkxaxes!(axs...)
    axs[2].xlabel = "Time, years"
    axs[1].ylabel = "Susceptible"
    axs[2].ylabel = "Infected"
    
    return fig
end 

end # module MID_63
