
module MID_44
  
using CairoMakie, DataFrames, DifferentialEquations
import Base: minimum

export Parameters44, sir44!, run_sir44, dataframe_sir44, plot_sir44, plot_sir44!

struct Parameters44
    r       :: Float64
    gamma   :: Vector{<:Float64}
    mu      :: Vector{<:Float64}
    nu      :: Vector{<:Float64}
    Thm     :: Float64
    Tmh     :: Float64
end

minimum(p::Parameters44) = min(p.r, minimum(p.gamma), minimum(p.mu), minimum(p.nu), p.Thm, p.Tmh)

function sir44!(du, u, p, t) 
    # compartments (including parameter results that are integrated)
    (Xh, Yh, Xm, Ym) = u  
    
    # ODEs 
    du[1, 1] = p.nu[1] - p.r * p.Thm * Ym * Xh - p.mu[1] * Xh       # dXh
    du[2, 1] = p.r * p.Thm * Ym * Xh - (p.gamma[1] + p.mu[1]) * Yh  # dYh 
    du[1, 2] = p.nu[2] - p.r * p.Tmh * Yh * Xm - p.mu[2] * Xm       # dXm
    du[2, 2] = p.r * p.Tmh * Yh * Xm - (p.gamma[2] + p.mu[2]) * Ym  # dYm 
    # The example code does not include the γₘ (p.gamma[2] in this code) term in 
    # the final equation but does take a value for γₘ. The term is included here. 
    # We will use p.gamma[2] = 0 but include the term in the equation to allow it 
    # to be generalised to other settings with mosquito recovery or to represent 
    # increased mortality among infected mosquitos
end 

function run_sir44(u0, p, duration; saveat = 1)
    @assert minimum(u0) >= 0 "Input u0 = $u0: cannot run with negative compartment values"
    @assert minimum(p) >= 0 "Input p = $p: cannot run with negative parameters"
    @assert maximum(p.Thm) <= 1 "Input p.Thm = $(p.Thm): is a probability and cannot be > 1"
    @assert maximum(p.Tmh) <= 1 "Input p.Tmh = $(p.Tmh): is a probability and cannot be > 1"
    @assert duration > 0 "Input duration = $duration: cannot run with zero or negative duration"

    tspan = ( 0., Float64(duration) )
    prob = ODEProblem(sir44!, u0, tspan, p)
    sol = solve(prob; saveat)
    return sol
end 

function run_sir44(; Xh, Yh, Xm, Ym, r, gamma, mu, nu, Thm, Tmh, duration, kwargs...)
    u0 = [ Xh Xm
           Yh Ym ]
    p = Parameters44(r, gamma, mu, nu, Thm, Tmh)
    return run_sir44(u0, p, duration; kwargs...)
end 

function dataframe_sir44(sol)
    result = DataFrame(t = Float64[], Xh = Float64[], Yh = Float64[], Xm = Float64[], Ym = Float64[])
    for i ∈ eachindex(sol.u)
        push!( result, Dict(
            :t => sol.t[i], 
            :Xh => sol.u[i][1, 1], 
            :Yh => sol.u[i][2, 1], 
            :Xm => sol.u[i][1, 2],
            :Ym => sol.u[i][2, 2]
        ) )
    end 
    return result 
end 

function plot_sir44(sol; kwargs...)
    fig = Figure()
    plot_sir44!(fig, sol; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_sir44!(any, sol; kwargs...) = plot_sir44!(any, dataframe_sir44(sol); kwargs...)

function plot_sir44!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_sir44!(gl, result; kwargs...)
end 

function plot_sir44!(gl::GridLayout, result::DataFrame; 
        label = "p4.4.jl: SIR model for mosquito vectors", kwargs...
    )
    axs = [ Axis(gl[i, 1]) for i ∈ 1:2 ]
    plot_sir44!(axs[1], select(result, :t, :Xh, :Yh))
    plot_sir44!(axs[2], select(result, :t, :Xm, :Ym))
    legs = [ Legend(gl[i, 2], axs[i]) for i ∈ 1:2 ]
    hidexdecorations!(axs[1]; ticks = false, grid = false)
    axs[2].xlabel = "Time"
    axs[1].ylabel = "Number of people"
    axs[2].ylabel = "Number of mosquitos"

 #   for (i, c) ∈ enumerate([:S, :E, :I])
 #       plot_seicr43!(axs[i], result, c)
 #       leg = Legend(gl[i, 2], axs[i])
 #       i < 3 && hidexdecorations!(axs[i]; ticks = false, grid = false)
 #       axs[i].ylabel = "Proportion $(String(c))" 
 #   end 
    axs[2].xlabel = "Time"
    Label(gl[0, 1], label)
end 

function plot_sir44!(ax::Axis, result::DataFrame)
    xs = result.t
    lines!(ax, xs, result[:, 2], label = names(result)[2])
    lines!(ax, xs, result[:, 3], label = names(result)[3])
end 

end # module MID_44
