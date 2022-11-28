
module MID_52
  
using CairoMakie, DataFrames, DifferentialEquations

export sir52!, run_sir52, dataframe_sir52, bifurcationdata_sir52, plot_sir52, plot_sir52!, bifurcationplot_sir52

function sir52!(du, u, p, t) 
    # compartments 
    S, I, R = u

    # Parameters 
    β0, β1, γ, μ, nu, term = p

    # Seasonal value of β
    β = β0 * (1 + β1 * term)
    
    # the ODEs
    du[1] = nu - β * S * I - μ * S      # dS
    du[2] = β * S * I - (γ + μ) * I     # dI
    du[3] = γ * I - μ * R               # dR
end 

function run_sir52!(u0, p, termstarttimes, termendtimes, duration; saveat = 1, kwargs...)
    @assert minimum(u0) >= 0 "Input u0 = $u0: cannot run with negative compartment values"
    @assert sum(u0) ≈ 1 "Input u0 = $u0: compartments are proportions so should sum to 1"
    @assert minimum(p[1:4]) >= 0 "Input p = $(p[1:4]): cannot run with negative parameters"
    @assert p[2] <= 1 "Input β1 = $(p[2]): if β1 > 1 then at some times the beta parameter will be negative"
    @assert duration > 0 "Input duration = $duration: cannot run with zero or negative duration"

    # Expand the lists of term start and end times to cover all years being simulated 
    i = 1 
    while maximum(termstarttimes) < duration 
        push!(termstarttimes, termstarttimes[i] + 365) 
        i +=1 
    end 
    i = 1 
    while maximum(termendtimes) < duration 
        push!(termendtimes, termendtimes[i] + 365) 
        i +=1 
    end 

    # Callbacks to switch between term and holiday
    termaffect!(integrator) = integrator.p[6] = 1
    termcb = PresetTimeCallback(termstarttimes, termaffect!; save_positions=(false, false))
    holidayaffect!(integrator) = integrator.p[6] = -1
    holidaycb = PresetTimeCallback(termendtimes, holidayaffect!; save_positions=(false, false))
    cbs = CallbackSet(termcb, holidaycb)

    tspan = ( 0., Float64(duration) )

    prob = ODEProblem(sir52!, u0, tspan, p)
    sol = solve(prob; saveat, callback = cbs)
    return sol
end 

function run_sir52(u0, p, termstarttimes, termendtimes, duration; kwargs...)
    parms = deepcopy(p)
    ts = deepcopy(termstarttimes)
    te = deepcopy(termendtimes)
    return run_sir52!(u0, parms, ts, te, duration; kwargs...)
end 

function run_sir52(; S0, I0, R0 = 1 - (S0 + I0), beta0, beta1, gamma, mu, nu = mu, 
        termstarttimes, termendtimes, duration, kwargs...
    )
    # Do we start in term time or in holiday?
    term = ifelse(termstarttimes[1] < termendtimes[1], -1, 1)

    u0 = [S0, I0, R0]
    p = [beta0, beta1, gamma, mu, nu, term]
    return run_sir52(u0, p, termstarttimes, termendtimes, duration; kwargs...)
end 

function dataframe_sir52(sol)
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

function bifurcationdata_sir52(beta0, beta1::Vector{<:Real}, gamma, mu, termstarttimes, termendtimes; 
        S0 = .5, I0 = 1e-4, R0 = 1 - (S0 + I0), kwargs...
    )
    # Do we start in term time or in holiday?
    term = ifelse(termstarttimes[1] < termendtimes[1], -1, 1)

    u0 = [S0, I0, R0]
    p = [beta0, beta1[1], gamma, mu, term]
    return bifurcationdata_sir52!(beta1, u0, p, termstarttimes, termendtimes; kwargs...)
end

bifurcationdata_sir52(; beta0, beta1, gamma, mu, termstarttimes, termendtimes, kwargs...) = 
    bifurcationdata_sir52(beta0, beta1, gamma, mu, termstarttimes, termendtimes; kwargs...)

function bifurcationdata_sir52!(beta1::Vector{<:Real}, u0, p, termstarttimes, termendtimes; kwargs...)
    data = DataFrame(beta1 = Float64[], t = Float64[], I = Float64[])
    ts = deepcopy(termstarttimes)
    te = deepcopy(termendtimes)
    bifurcationdata_sir52!(data, beta1, u0, p, ts, te; kwargs...) 
    return data 
end

function bifurcationdata_sir52!(data, beta1::Vector{<:Real}, u0, p, termstarttimes, termendtimes; kwargs...) 
    for β1 ∈ beta1
        bifurcationdata_sir52!(data, β1, u0, p, termstarttimes, termendtimes; kwargs...) 
    end
end 

function bifurcationdata_sir52!(data, beta1::Real, u0, p, termstarttimes, termendtimes; 
        runin = 100, plot = 100, kwargs...
    ) 
    p[2] = beta1
    duration = (runin + plot) * 365
    sol = run_sir52!(u0, p, termstarttimes, termendtimes, duration; saveat = 365, kwargs...)
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

function plot_sir52(result; kwargs...)
    fig = Figure()
    plot_sir52!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_sir52!(any, sol; kwargs...) = plot_sir52!(any, dataframe_sir52(sol); kwargs...)

function plot_sir52!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_sir52!(gl, result; kwargs...)
end 

function plot_sir52!(gl::GridLayout, result::DataFrame; 
        label = "p5.2.jl: SIR model with corrected term-time forcing", 
        legend = :right, kwargs...
    )
    ax = Axis(gl[1, 1])
    plot_sir52!(ax, result; kwargs...)

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

function plot_sir52!(ax::Axis, result::DataFrame; plotR = false)
    xs = result.t ./ 365  # to plot time in weeks
    
    lines!(ax, xs, result.S, label = "Susceptible")
    lines!(ax, xs, result.I, label = "Infectious")
    plotR && lines!(ax, xs, result.R, label = "Recovered")
    ax.xlabel = "Time, years"
    ax.ylabel = "Fraction of population"
end 

function bifurcationplot_sir52(data; label = "p5.2.jl: SIR model with corrected term-time forcing")
    fig = Figure()
    ax = Axis(fig[1, 1])
    scatter!(ax, data.beta1, data.I; marker = :rect, markersize = 3)
    ax.xlabel = "β1 (amplitude of term-time forcing)"
    ax.ylabel = "Annual infectiveness"
    Label(fig[0, :], label)
    return fig
end

end # module MID_52
