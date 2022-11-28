
module MID_66

# SIR model with two types of imports (page 210)
  
using CairoMakie, DataFrames, Distributions, Random

export sir66, u0_sir66, run_sir66, plot_sir66

const changematrix = # how contents of each compartment change for each possible event     
    [ -1   1   0    # infection 
       0  -1   1    # recovery
       1   0   0    # birth 
      -1   0   0    # susceptible death 
       0  -1   0    # infectious death 
       0   0  -1    # recovered death
       0   1   0    # infectious immigration 
      -1   1   0 ]  # external infection

function sir66(u, p, t, δt, N0) 
    # The model does not enforce a constant population. However, an approximately 
    # constant population is desired. If we calculated N = X + Y + Z then the birth 
    # rate would be proportional to the current population and could lead to a positive 
    # feedback. This function therefore takes an input of N0 to use in calculating 
    # birth rate.
    X, Y, Z = u 
    β, γ, δ, ε, μ = p 
    N = X + Y + Z

    rates = [           # rates at with each possible event occur
        β * X * Y / N,  # infection
        γ * Y,          # recovery 
        μ * N0,         # birth 
        μ * X,          # susceptible death 
        μ * Y,          # infectious death 
        μ * Z,          # recovered death
        δ,              # infectious immigration 
        ε * X           # external infection
    ]

    events = [ rand(Poisson(rate * δt)) for rate ∈ rates ]

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

function run_sir66(; N0, beta, gamma, delta, epsilon, mu, duration, kwargs...)
    p = [beta, gamma, delta, epsilon, mu]
    u0 = u0_sir66(N0, p)
    return run_sir66(u0, p, duration; kwargs...)
end

run_sir66(u0, p, duration; δt = 1, seed = nothing) = _run_sir66(u0, p, duration, seed; δt)

function _run_sir66(u0, p, duration, seed::Int; δt)
    Random.seed!(seed)
    return _run_sir66(u0, p, duration, nothing; δt)
end 

_run_sir66(N0::Int, p, duration, seed::Nothing; δt) = 
    _run_sir65(u0_sir66(N0, p), p, duration, seed; δt)

function _run_sir66(u0::Vector{<:Int}, p, duration, seed::Nothing; δt)
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
        t, u = sir66(u, p, t, δt, N0) 
        push!( results, Dict(:t => t, :X => u[1], :Y => u[2], :Z => u[3]) )
    end

    return results
end 

function u0_sir66(N0::Int, p) 
    β, γ, δ, ε, μ = p
    X0 = round(Int, γ * N0 / β, RoundDown)
    Y0 = round(Int, μ * N0 / γ, RoundUp) 
    Z0 = N0 - ( X0 + Y0 ) 
    return [X0, Y0, Z0]
end 

function plot_sir66(results)
    return plot_sir66(
        results, 
        "p6.6.jl: SIR model with τ-leap stochasticity and possible disease importation"
    )
end 

function plot_sir66(results, population::Real)
    return plot_sir66(
        results, 
        "p6.6.jl: SIR model with τ-leap stochasticity and possible disease importation
        Initial population = $population"
    )
end 

function plot_sir66(results, label::String)
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
    Label(fig[0, 1], label; justification = :left)
    
    return fig
end 

end # module MID_66
