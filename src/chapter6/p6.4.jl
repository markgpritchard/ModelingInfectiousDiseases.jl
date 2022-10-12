
module MID_64
  
using CairoMakie, DataFrames, Random

export sir64!, u0_sir64, run_sir64, plot_sir64

const changematrix =    # how contents of each compartment change for each possible event     
    [   -1  1   0       # infection 
        0   -1  1       # recovery
        1   0   0       # birth 
        -1  0   0       # susceptible death 
        0   -1  0       # infectious death 
        0   0   -1  ]   # recovered death

function u0_sir64(N0::Int, p) 
    β, γ, μ = p
    X0 = round(Int, γ * N0 / β, RoundDown)
    Y0 = round(Int, μ * N0 / γ, RoundUp) 
    Z0 = N0 - ( X0 + Y0 ) 
    return [X0, Y0, Z0]
end 

function sir64!(u, p, t, N0) 
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
run_sir64(u0, p, duration; seed = nothing) = _run_sir64(u0, p, duration, seed)

function _run_sir64(u0, p, duration, seed::Int)
    Random.seed!(seed)
    return _run_sir64(u0, p, duration, nothing)
end 

_run_sir64(N0::Int, p, duration, seed::Nothing) = _run_sir64(u0_sir64(N0, p), p, duration, seed)

function _run_sir64(u0::Vector{<:Int}, p, duration, seed::Nothing)
    @assert minimum(u0) >= 0 "Model cannot run with negative starting values in `u0`. Model supplied u0 = $u0."
    @assert minimum(p) < 0 "Model cannot run with negative parameters. Running with p = $p."
    @assert duration > 0 "Model needs duration > 0. Model supplied duration = $duration."

    N0 = sum(u0)
    t = .0 
    u = u0
    results = DataFrame(
        t = Float64[], X = typeof(u0[1])[], Y = typeof(u0[2])[], Z = typeof(u0[3])[]
    )
    push!( results, Dict(:t => t, :X => u[1], :Y => u[2], :Z => u[3]) )
    while t <= duration 
        t, u = sir64!(u, p, t, N0) 
        push!( results, Dict(:t => t, :X => u[1], :Y => u[2], :Z => u[3]) )
    end

    return results
end 


"""
    plot_sir61(results)

Plot the `results` DataFrame output from the function `run_sir61` 
"""
function plot_sir64(results)
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
    
    return fig
end 

end # module MID_64
