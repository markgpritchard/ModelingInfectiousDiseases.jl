
module MID_63

# SIS model with demographic stochasticity (page 202)
  
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

function run_sis63(; X0::Int, Y0::Int, beta, gamma, duration, kwargs...)
    u0 = [X0, Y0] 
    p = [beta, gamma]
    return run_sis63(u0, p, duration; kwargs...)
end 

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
