
module MID_82
  
using CairoMakie, DataFrames, DifferentialEquations

export sir82!, run_sir82, dataframe_sir82, plot_sir82, plot_sir82!

function sir82!(du, u, p, t) 
    # compartments 
    S, I, R = u

    # parameters 
    beta, gamma, mu, nu, v = p
    
    # ODEs
    du[1] = nu - (beta * I + mu + v) * S        # dS 
    du[2] = beta * I * S - (gamma + mu) * I     # dI 
    du[3] = gamma * I + v * S - mu * R          # dR
end 

function run_sir82(u0, p, duration, vaccinationstarttime, vaccinationrate; saveat = 1)
    @assert minimum(u0) >= 0 "Input u0 = $u0: no compartments can contain negative proportions"
    @assert minimum(p) >= 0 "Input p = $p: cannot run with negative parameters" 
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"
    
    tspan = ( 0., Float64(duration) )
    prob = ODEProblem(sir82!, u0, tspan, p)
    
    if vaccinationstarttime > duration 
        @info "Vaccination start time is after end of model run" 
        return solve(prob; saveat)
    else 
        # The callback to change the parameter for vaccine rate at vaccinationstarttime
        affect!(integrator) = integrator.p[5] = vaccinationrate
        cb = PresetTimeCallback(vaccinationstarttime, affect!)
        return solve(prob; callback = cb, saveat)
    end
end 

function dataframe_sir82(sol)
    result = DataFrame(t = sol.t)
    insertcols!(result, :S => [ sol[i][1] for i ∈ axes(sol, 2) ])
    insertcols!(result, :I => [ sol[i][2] for i ∈ axes(sol, 2) ])
    insertcols!(result, :R => [ sol[i][3] for i ∈ axes(sol, 2) ])
    return result 
end 

function plot_sir82(result; kwargs...)
    fig = Figure()
    plot_sir82!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_sir82!(any, sol; kwargs...) = plot_sir82!(any, dataframe_sir82(sol); kwargs...)

function plot_sir82!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_sir82!(gl, result; kwargs...)
end 

function plot_sir82!(gl::GridLayout, result::DataFrame; 
        label = "p8.2.jl: SIR model with random vaccination", kwargs...)

    ax = Axis(gl[1, 1])
    plot_sir82!(ax, result; kwargs...)
    leg = Legend(gl[1, 2], ax)
end 

function plot_sir82!(ax::Axis, result::DataFrame; plotr = true)
    lines!(ax, result.t / 365, result.S; label = "Susceptible")
    lines!(ax, result.t / 365, result.I; label = "Infectious")
    if plotr lines!(ax, result.t / 365, result.R; label = "Resistant") end
    ax.xlabel = "Time, years"
    ax.ylabel = "Proportion"
end 

end # module MID_82
