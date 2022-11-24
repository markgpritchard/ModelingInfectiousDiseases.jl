
module MID_51
  
using CairoMakie, DataFrames, DifferentialEquations

export sir51!, run_sir51, dataframe_sir51, bifurcationdata_sir51, plot_sir51, plot_sir51!, bifurcationplot_sir51

function sir51!(du, u, p, t) 
    # compartments 
    S, I, R = u

    # Parameters 
    β0, β1, γ, μ, ω = p

    # Seasonal value of β
    β = β0 * (1 + β1 * sin(ω * t))
    
    # the ODEs
    du[1] = μ - β * S * I - μ * S       # dS
    du[2] = β * S * I - (γ + μ) * I     # dI
    du[3] = γ * I - μ * R               # dR
end 

function run_sir51(u0, p, duration; saveat = 1, kwargs...)
    @assert minimum(u0) >= 0 "Input u0 = $u0: cannot run with negative compartment values"
    @assert sum(u0) ≈ 1 "Input u0 = $u0: compartments are proportions so should sum to 1"
    @assert minimum(p) >= 0 "Input p = $p: cannot run with negative parameters"
    @assert p[2] <= 1 "Input β1 = $(p[2]): if β1 > 1 then at some times the beta parameter will be negative"
    @assert duration > 0 "Input duration = $duration: cannot run with zero or negative duration"

    tspan = ( 0., Float64(duration) )

    prob = ODEProblem(sir51!, u0, tspan, p)
    sol = solve(prob; saveat, kwargs...)

    return sol
end 

function run_sir51(; S0, I0, R0 = 1 - (S0 + I0), beta0, beta1, gamma, mu, omega = 2pi / 365, duration, kwargs...)
    u0 = [S0, I0, R0]
    p = [beta0, beta1, gamma, mu, omega]
    return run_sir51(u0, p, duration; kwargs...)
end 

function dataframe_sir51(sol)
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

function bifurcationdata_sir51(beta0, beta1::Vector{<:Real}, gamma, mu, omega = 2pi / 365; 
        S0 = .5, I0 = 1e-4, R0 = 1 - (S0 + I0), kwargs...
    )
    u0 = [S0, I0, R0]
    p = [beta0, beta1[1], gamma, mu, omega]
    return bifurcationdata_sir51!(beta1, u0, p; kwargs...)
end

bifurcationdata_sir51(; beta0, beta1, gamma, mu, omega = 2pi / 365, kwargs...) = 
    bifurcationdata_sir51(beta0, beta1, gamma, mu, omega; kwargs...)

function bifurcationdata_sir51!(beta1::Vector{<:Real}, u0, p; kwargs...)
    data = DataFrame(beta1 = Float64[], t = Float64[], I = Float64[])
    bifurcationdata_sir51!(data, beta1, u0, p; kwargs...) 
    return data 
end

function bifurcationdata_sir51!(data, beta1::Vector{<:Real}, u0, p; kwargs...) 
    for β1 ∈ beta1
        bifurcationdata_sir51!(data, β1, u0, p; kwargs...) 
    end
end 

function bifurcationdata_sir51!(data, beta1::Real, u0, p; runin = 100, plot = 100, kwargs...) 
    p[2] = beta1
    duration = (runin + plot) * 365
    sol = run_sir51(u0, p, duration; saveat = 365, kwargs...)
    try
        for i ∈ runin+1:runin+plot
            push!(data, Dict(
                :beta1 => beta1,
                :t => sol.t[i], 
                :I => sol.u[i][2]
            ) )
        end 
    catch e
        @warn "Error when beta1 = $beta1: $e"
    end 
end 

function plot_sir51(result; kwargs...)
    fig = Figure()
    plot_sir51!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_sir51!(any, sol; kwargs...) = plot_sir51!(any, dataframe_sir51(sol); kwargs...)

function plot_sir51!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_sir51!(gl, result; kwargs...)
end 

function plot_sir51!(gl::GridLayout, result::DataFrame; 
        label = "p5.1.jl: SIR model with sinusoidal forcing", 
        legend = :right, kwargs...
    )
    ax = Axis(gl[1, 1])
    plot_sir51!(ax, result; kwargs...)

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

function plot_sir51!(ax::Axis, result::DataFrame; plotR = false)
    xs = result.t ./ 365  # to plot time in weeks
    
    lines!(ax, xs, result.S, label = "Susceptible")
    lines!(ax, xs, result.I, label = "Infectious")
    plotR && lines!(ax, xs, result.R, label = "Recovered")
    ax.xlabel = "Time, years"
    ax.ylabel = "Fraction of population"
end 

function bifurcationplot_sir51(data; label = "p5.1.jl: SIR model with sinusoidal forcing")
    fig = Figure()
    ax = Axis(fig[1, 1])
    scatter!(ax, data.beta1, data.I; marker = :rect, markersize = 3)
    ax.xlabel = "β1 (amplitude of sinuoidal forcing)"
    ax.ylabel = "Annual infectiveness"
    Label(fig[0, :], label)
    return fig
end

end # module MID_51
