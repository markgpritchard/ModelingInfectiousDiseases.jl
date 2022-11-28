
module MID_83

# SIR model with pulsed vaccination (page 302)
  
using CairoMakie, DataFrames, DifferentialEquations

export sir83!, run_sir83, dataframe_sir83, plot_sir83, plot_sir83!

function sir83!(du, u, p, t) 
    # compartments 
    S, I, R = u
    # parameters 
    beta, gamma, mu, nu = p
    
    # ODEs without vaccination (that is added in pulses via callbacks)
    du[1] = nu - beta * S * I - mu * S          # dS 
    du[2] = beta * S * I - (gamma + mu) * I     # dI 
    du[3] = gamma * I - mu * R                  # dR
end 

function run_sir83(; S0, I0, R0 = 1 - (S0 + I0), beta, gamma, mu, nu = mu, duration, 
        vaccinationstarttime, vaccinationfrequency, vaccinationproportion, kwargs...
    )
    u0 = [S0, I0, R0] 
    p = [beta, gamma, mu, nu]
    return run_sir83(u0, p, duration, vaccinationstarttime, vaccinationfrequency, vaccinationproportion; kwargs...)
end 

function run_sir83(u0, p, duration, vaccinationstarttime, vaccinationfrequency, vaccinationproportion; 
        saveat = 1
    )
    @assert minimum(u0) >= 0 "Input u0 = $u0: no compartments can contain negative proportions"
    @assert sum(u0) ≈ 1 "Input sum(u0) = $(sum(u0)): compartments are proportions so should sum to 1"
    @assert minimum(p) >= 0 "Input p = $p: cannot run with negative parameters" 
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"
    @assert vaccinationstarttime >= 0 "Input vaccinationstarttime = $vaccinationstarttime: cannot run with negative vaccination start time"
    @assert vaccinationfrequency > 0 "Input vaccinationfrequency = $vaccinationfrequency: cannot run with negative or zero vaccination frequency"
    @assert vaccinationproportion >= 0 "Input vaccinationproportion = $vaccinationproportion: cannot run with negative vaccination rate"
    
    tspan = ( 0., Float64(duration) )
    prob = ODEProblem(sir83!, u0, tspan, p)
    
    if vaccinationstarttime > duration 
        @info "Vaccination start time is after end of model run" 
        return solve(prob; saveat)
    else 
        # The callback to change the parameter for vaccine rate at vaccinationstarttime
        vaccinetimes = collect(vaccinationstarttime:vaccinationfrequency:duration)
        function affect!(integrator) 
            integrator.u[1] = integrator.u[1] * (1 - vaccinationproportion)
            integrator.u[3] += integrator.u[1] * vaccinationproportion
        end 
        cb = PresetTimeCallback(vaccinetimes, affect!)
        return solve(prob; callback = cb, saveat)
    end
end 

function dataframe_sir83(sol)
    return DataFrame(
        t = sol.t,
        S = [ sol[i][1] for i ∈ axes(sol, 2) ],
        I = [ sol[i][2] for i ∈ axes(sol, 2) ],
        R = [ sol[i][3] for i ∈ axes(sol, 2) ]
    ) 
end 

function plot_sir83(result; kwargs...)
    fig = Figure()
    plot_sir83!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_sir83!(any, sol; kwargs...) = plot_sir83!(any, dataframe_sir83(sol); kwargs...)

function plot_sir83!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_sir83!(gl, result; kwargs...)
end 

function plot_sir83!(gl::GridLayout, result::DataFrame; 
        label = "p8.3.jl: SIR model with pulsed vaccination", kwargs...
    )
    ax = Axis(gl[1, 1])
    plot_sir83!(ax, result; kwargs...)
    leg = Legend(gl[1, 2], ax)
    Label(gl[0, :], label)
end 

function plot_sir83!(ax::Axis, result::DataFrame; plotr = true)
    lines!(ax, result.t / 365, result.S; label = "Susceptible")
    lines!(ax, result.t / 365, result.I; label = "Infectious")
    if plotr lines!(ax, result.t / 365, result.R; label = "Resistant") end
    ax.xlabel = "Time, years"
    ax.ylabel = "Proportion"
end 

end # module MID_83
