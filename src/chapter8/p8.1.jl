
module MID_81
  
using CairoMakie, DataFrames, DifferentialEquations

export sir81!, run_sir81, run_sir81!, dataframe_sir81, plot_sir81, plot_sir81!

function sir81!(du, u, p, t) 
    # compartments 
    S, I, R = u

    # parameters 
    beta, gamma, mu, nu, pr = p
    
    # ODEs
    du[1] = nu * (1 - pr) - beta * I * S - mu * S   # dS 
    du[2] = beta * I * S - (gamma + mu) * I         # dI 
    du[3] = gamma * I + nu * pr - mu * R            # dR
end 

function run_sir81(u0, p, duration, vaccinationstarttime, vaccinationrate; saveat = 1)
    pcopy = deepcopy(p)
    return run_sir81!(u0, pcopy, duration, vaccinationstarttime, vaccinationrate; saveat)
end 

function run_sir81!(u0, p, duration, vaccinationstarttime, vaccinationrate; saveat = 1)
    @assert minimum(u0) >= 0 "Input u0 = $u0: no compartments can contain negative proportions"
    @assert sum(u0) ≈ 1 "Input sum(u0) = $(sum(u0)): compartments are proportions so should sum to 1"
    @assert minimum(p) >= 0 "Input p = $p: cannot run with negative parameters" 
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"
    
    tspan = ( 0., Float64(duration) )
    prob = ODEProblem(sir81!, u0, tspan, p)
    
    if vaccinationstarttime > duration 
        @info "Vaccination start time is after end of model run" 
        sol = solve(prob; saveat)
    else 
        # The callback to change the parameter for vaccine rate at vaccinationstarttime
        affect!(integrator) = integrator.p[5] = vaccinationrate
        cb = PresetTimeCallback(vaccinationstarttime, affect!)
        sol = solve(prob; callback = cb, saveat)
    end

    return sol
end 

function dataframe_sir81(sol)
    result = DataFrame(t = sol.t)
    insertcols!(result, :S => [ sol[i][1] for i ∈ axes(sol, 2) ])
    insertcols!(result, :I => [ sol[i][2] for i ∈ axes(sol, 2) ])
    insertcols!(result, :R => [ sol[i][3] for i ∈ axes(sol, 2) ])
    return result 
end 

function plot_sir81(result; kwargs...)
    fig = Figure()
    plot_sir81!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_sir81!(any, sol; kwargs...) = plot_sir81!(any, dataframe_sir81(sol); kwargs...)

function plot_sir81!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_sir81!(gl, result; kwargs...)
end 

function plot_sir81!(gl::GridLayout, result::DataFrame; 
        label = "p8.1.jl: SIR model with childhood vaccination", kwargs...)

    ax = Axis(gl[1, 1])
    plot_sir81!(ax, result; kwargs...)
    leg = Legend(gl[1, 2], ax)
    Label(gl[0, 1], label)
end 

function plot_sir81!(ax::Axis, result::DataFrame; plotr = true)
    lines!(ax, result.t / 365, result.S; label = "Susceptible")
    lines!(ax, result.t / 365, result.I; label = "Infectious")
    if plotr lines!(ax, result.t / 365, result.R; label = "Resistant") end
    ax.xlabel = "Time, years"
    ax.ylabel = "Proportion"
end 

end # module MID_81
