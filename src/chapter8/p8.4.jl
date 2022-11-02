
module MID_84
  
using CairoMakie, DataFrames, DifferentialEquations
import Base: minimum

export Parameters84, sir84!, run_sir84, dataframe_sir84, plot_sir84, plot_sir84!

mutable struct Parameters84 # beta will be a matrix but gamma is a float, so group them in a structure
    beta    :: Matrix{<:Float64} 
    gamma   :: Float64
    mu      :: Float64
    nu      :: Vector{<:Float64}
    pv      :: Vector{<:Float64}
end

minimum(p::Parameters84) = min(minimum(p.beta), p.gamma, p.mu, minimum(p.nu), minimum(p.pv))

function sir84!(du, u, p, t) 
    # compartments 
    Sh, Ih, Rh, Sℓ, Iℓ, Rℓ = u

    # parameters 
    (βhh, βhl, βlh, βll) = p.beta 
    γ = p.gamma
    μ = p.mu 
    nuh, nul = p.nu 
    ph, pl = p.pv
    
    # the ODEs
    # high risk
    du[1] = nuh * (1 - ph) - (βhh * Ih + βhl * Iℓ) * Sh - μ * Sh    # dSh
    du[2] = (βhh * Ih + βhl * Iℓ) * Sh - (γ + μ) * Ih               # dIh
    du[3] = nuh * ph + γ * Ih - μ * Rh                              # dRh

    # low risk
    du[4] = nul * (1 - pl) - (βlh * Ih + βll * Iℓ) * Sℓ - μ * Sℓ    # dSℓ
    du[5] = (βlh * Ih + βll * Iℓ) * Sℓ - (γ + μ) * Iℓ               # dIℓ
    du[6] = nul * pl + γ * Iℓ - μ * Rℓ                              # dRℓ
end 

function run_sir84(u0, p, duration, vaccinationstarttime, vaccinationrate; saveat = 1)
    @assert minimum(u0) >= 0 "Input u0 = $u0: no compartments can contain negative proportions"
    @assert sum(u0) ≈ 1 "Input sum(u0) = $(sum(u0)): compartments are proportions so should sum to 1"
    @assert minimum(p) >= 0 "Input p = $p: cannot run with negative parameters"
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"
    @assert minimum(vaccinationrate) >= 0 "Input vaccinationrate = $vaccinationrate: cannot run with negative vaccination rate"

    tspan = ( 0., Float64(duration) )
    prob = ODEProblem(sir84!, u0, tspan, p)
    
    if vaccinationstarttime > duration 
        @info "Vaccination start time is after end of model run" 
        return solve(prob; saveat)
    else 
        # The callback to change the parameter for vaccine rate at vaccinationstarttime
        affect!(integrator) = integrator.p.pv = vaccinationrate
        cb = PresetTimeCallback(vaccinationstarttime, affect!)
        return solve(prob; callback = cb, saveat)
    end
end 

function dataframe_sir84(sol)
    return DataFrame(
        t = sol.t,
        Sh = [ sol[i][1, 1] for i ∈ axes(sol, 3) ],
        Ih = [ sol[i][2, 1] for i ∈ axes(sol, 3) ],
        Rh = [ sol[i][3, 1] for i ∈ axes(sol, 3) ],
        Sl = [ sol[i][1, 2] for i ∈ axes(sol, 3) ],
        Il = [ sol[i][2, 2] for i ∈ axes(sol, 3) ],
        Rl = [ sol[i][3, 2] for i ∈ axes(sol, 3) ]
    ) 
end 

function plot_sir84(result; kwargs...)
    fig = Figure()
    plot_sir84!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_sir84!(any, sol; kwargs...) = plot_sir84!(any, dataframe_sir84(sol); kwargs...)

function plot_sir84!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_sir84!(gl, result; kwargs...)
end 

function plot_sir84!(gl::GridLayout, result::DataFrame; 
        label = "p8.4.jl: SIR model with high and low risk groups", plotr = true)

    ax1 = Axis(gl[1, 1])
    lines!(ax1, result.t / 365, result.Sh; label = "Susceptible high risk")
    lines!(ax1, result.t / 365, result.Ih; label = "Infectious high risk")
    if plotr 
        lines!(ax1, result.t / 365, result.Rh; label = "Resistant high risk") 
    end

    ax2 = Axis(gl[2, 1])
    lines!(ax2, result.t / 365, result.Sl; label = "Susceptible low risk")
    lines!(ax2, result.t / 365, result.Il; label = "Infectious low risk")
    if plotr 
        lines!(ax2, result.t / 365, result.Rl; label = "Resistant low risk") 
    end
    linkxaxes!(ax1, ax2)
    hidexdecorations!(ax1; ticks = false)
    ax2.xlabel = "Time, years"
    ylbl = Label(gl[:, 0], "Proportion"; rotation = pi/2)
    leg1 = Legend(gl[1, 2], ax1); leg2 = Legend(gl[2, 2], ax2)
end 

end # module MID_84
