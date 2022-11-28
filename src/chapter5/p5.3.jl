
module MID_53

# SIR model with sinusoidal births (page 184)
  
using CairoMakie, DataFrames, DifferentialEquations

export sir53!, run_sir53, dataframe_sir53, bifurcationdata_sir53, plot_sir53, plot_sir53!, bifurcationplot_sir53

function sir53!(du, u, p, t) 
    # compartments 
    S, I, R = u

    # Parameters 
    α0, α1, β, γ, μ, ω = p

    # Seasonal value of α
    α = α0 * (1 + α1 * sin(ω * t))
    
    # the ODEs
    du[1] = α - β * S * I - μ * S       # dS
    du[2] = β * S * I - (γ + μ) * I     # dI
    du[3] = γ * I - μ * R               # dR
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions for a single value of alpha1 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function run_sir53(; S0, I0, R0 = 1 - (S0 + I0), alpha0, alpha1, beta, gamma, mu = alpha0, 
        omega = 2pi / 365, duration, kwargs...
    )
    u0 = [S0, I0, R0]
    p = [alpha0, alpha1, beta, gamma, mu, omega]
    return run_sir53(u0, p, duration; kwargs...)
end 

function run_sir53(u0, p, duration; saveat = 1, kwargs...)
    @assert minimum(u0) >= 0 "Input u0 = $u0: cannot run with negative compartment values"
    @assert sum(u0) ≈ 1 "Input u0 = $u0: compartments are proportions so should sum to 1"
    @assert minimum(p) >= 0 "Input p = $p: cannot run with negative parameters"
    @assert p[2] <= 1 "Input α1 = $(p[2]): if α1 > 1 then at some times the beta parameter will be negative"
    @assert duration > 0 "Input duration = $duration: cannot run with zero or negative duration"

    tspan = ( 0., Float64(duration) )

    prob = ODEProblem(sir53!, u0, tspan, p)
    sol = solve(prob; saveat, kwargs...)

    return sol
end 

function dataframe_sir53(sol)
    result = DataFrame(t = Float64[], S = Float64[], I = Float64[], R = Float64[])
    for i ∈ eachindex(sol.u)
        push!( result, Dict(
            :t => sol.t[i], 
            :S => sol.u[i][1], 
            :I => sol.u[i][2], 
            :R => sol.u[i][3]
        ) )
    end 
    return result 
end 

function plot_sir53(result; kwargs...)
    fig = Figure()
    plot_sir53!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_sir53!(any, sol; kwargs...) = plot_sir53!(any, dataframe_sir53(sol); kwargs...)

function plot_sir53!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_sir53!(gl, result; kwargs...)
end 

function plot_sir53!(gl::GridLayout, result::DataFrame; 
        label = "p5.3.jl: SIR model with sinusoidal births", 
        legend = :right, kwargs...
    )
    ax = Axis(gl[1, 1])
    plot_sir53!(ax, result; kwargs...)

    Label(gl[0, :], label)

    if legend == :right
        leg = Legend(gl[1, 2], ax)
    elseif legend == :below 
        leg = Legend(gl[2, 1], ax)
    elseif legend == :none 
        # no legend 
    else 
        @info "Unrecognised legend position given. Recognised options are `:right`, `:below` and `:none`"
    end
end 

function plot_sir53!(ax::Axis, result::DataFrame; plotR = false)
    xs = result.t ./ 365  # to plot time in weeks
    
    lines!(ax, xs, result.S, label = "Susceptible")
    lines!(ax, xs, result.I, label = "Infectious")
    plotR && lines!(ax, xs, result.R, label = "Recovered")
    ax.xlabel = "Time, years"
    ax.ylabel = "Fraction of population"
end 


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Functions to produce a bifurcation diagram (a vector for alpha1)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

bifurcationdata_sir53(; alpha0, alpha1, beta, gamma, mu = alpha0, omega = 2pi / 365, kwargs...) = 
    bifurcationdata_sir53(alpha0, alpha1, beta, gamma, mu, omega; kwargs...)

function bifurcationdata_sir53(alpha0, alpha1::Vector{<:Real}, beta, gamma, mu = alpha0, omega = 2pi / 365; 
        S0 = .5, I0 = 1e-4, R0 = 1 - (S0 + I0), kwargs...
    )
    u0 = [S0, I0, R0]
    p = [alpha0, alpha1[1], beta, gamma, mu, omega]
    return bifurcationdata_sir53!(alpha1, u0, p; kwargs...)
end

function bifurcationdata_sir53!(alpha1::Vector{<:Real}, u0, p; kwargs...)
    data = DataFrame(alpha1 = Float64[], t = Float64[], I = Float64[])
    bifurcationdata_sir53!(data, alpha1, u0, p; kwargs...) 
    return data 
end

function bifurcationdata_sir53!(data, alpha1::Vector{<:Real}, u0, p; kwargs...) 
    for α1 ∈ alpha1
        bifurcationdata_sir53!(data, α1, u0, p; kwargs...) 
    end
end 

function bifurcationdata_sir53!(data, alpha1::Real, u0, p; runin = 100, plot = 100, kwargs...) 
    p[2] = alpha1
    duration = (runin + plot) * 365
    sol = run_sir53(u0, p, duration; saveat = 365, kwargs...)
    try
        for i ∈ runin+1:runin+plot
            push!(data, Dict(
                :alpha1 => alpha1,
                :t => sol.t[i], 
                :I => sol.u[i][2]
            ) )
        end 
    catch e
        @warn "Error when alpha1 = $alpha1: $e"
    end 
end 

function bifurcationplot_sir53(data; label = "p5.3.jl: SIR model with sinusoidal births")
    fig = Figure()
    ax = Axis(fig[1, 1])
    scatter!(ax, data.alpha1, data.I; marker = :rect, markersize = 3)
    ax.xlabel = "α1 (amplitude of sinuoidal birth rate)"
    ax.ylabel = "Annual infectiveness"
    Label(fig[0, :], label)
    return fig
end

end # module MID_53
