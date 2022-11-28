
module MID_41

# SIR model with partial immunity (page 118)
  
using CairoMakie, DataFrames, DifferentialEquations
import Base: minimum

export Parameters41, sir41!, run_sir41, dataframe_sir41, plot_sir41, plot_sir41!

struct Parameters41 
    a       :: Vector{<:Float64} 
    alpha   :: Vector{<:Float64} 
    beta    :: Vector{<:Float64} 
    gamma   :: Vector{<:Float64} 
    mu      :: Float64 
    nu      :: Float64
end

function sir41!(du, u, p, t) 
    # compartments 
    SS, IS, RS, SI, RI, SR, IR, RR = u

    I1 = IS + p.a[1] * IR 
    I2 = SI + p.a[2] * RI
    
    # the ODEs
    du[1] = p.nu - ( p.beta[1] * I1 + p.beta[2] * I2 + p.mu ) * SS          # dSS  
    du[2] = p.beta[1] * I1 * SS - ( p.gamma[1] + p.mu ) * IS                # dIS 
    du[3] = p.gamma[1] * IS - ( p.alpha[2] * p.beta[2] * I2 + p.mu ) * RS   # dRS 
    du[4] = p.beta[2] * I2 * SS - ( p.gamma[2] + p.mu ) * SI                # dSI 
    du[5] = p.alpha[2] * p.beta[2] * RS * I2 - ( p.gamma[2] + p.mu ) * RI   # dRI 
    du[6] = p.gamma[1] * IS - ( p.alpha[1] * p.beta[1] * I1 + p.mu ) * SR   # dSR 
    du[7] = p.alpha[1] * p.beta[1] * SR * I1 - ( p.gamma[1] + p.mu ) * IR   # dRI 
    du[8] = p.gamma[1] * IR + p.gamma[2] * RI - p.mu * RR                   # dRR  
end 

function run_sir41(; SS0, IS0, RS0, SI0, RI0, SR0, IR0, RR0, a, alpha, beta, gamma, mu, nu = mu, duration, kwargs...)
    u0 = [SS0, IS0, RS0, SI0, RI0, SR0, IR0, RR0]
    p = Parameters41(a, alpha, beta, gamma, mu, nu)
    return run_sir41(u0, p, duration; kwargs...)
end 

function run_sir41(u0, p, duration; saveat = 1)
    @assert minimum(u0) >= 0 "Input u0 = $u0: cannot run with negative compartment values"
    @assert sum(u0) ≈ 1 "Input u0 = $u0: compartments are proportions so should sum to 1"
    @assert minimum(p) >= 0 "Input p = $p: cannot run with negative parameters"
    @assert duration > 0 "Input duration = $duration: cannot run with zero or negative duration"

    tspan = ( 0., Float64(duration) )

    prob = ODEProblem(sir41!, u0, tspan, p)
    sol = solve(prob; saveat)
    return sol
end 

minimum(p::Parameters41) = min( minimum(p.a), minimum(p.alpha), minimum(p.beta), 
    minimum(p.gamma), p.mu, p.nu )

function dataframe_sir41(sol)
    result = DataFrame(t = sol.t)
    for (i, lbl) ∈ enumerate([ "SS", "IS", "RS", "SI", "RI", "SR", "IR", "RR" ])
        insertcols!(result, Symbol(lbl) => [ sol[j][i] for j ∈ axes(sol, 2) ])
    end 
    return result 
end 

function plot_sir41(result; kwargs...)
    fig = Figure()
    plot_sir41!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_sir41!(any, sol; kwargs...) = plot_sir41!(any, dataframe_sir41(sol); kwargs...)

function plot_sir41!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_sir41!(gl, result; kwargs...)
end 

function plot_sir41!(gl::GridLayout, result::DataFrame; 
        label = "p4.1.jl: SIR model with partial immunity", kwargs...
    )
    ax = Axis(gl[1, 1])
    plot_sir41!(ax, result; kwargs...)
    leg = Legend(gl[1, 2], ax)
    Label(gl[0, 1], label)
end 

function plot_sir41!(ax::Axis, result::DataFrame)
    for lbl ∈ [ "SS", "IS", "RS", "SI", "RI", "SR", "IR", "RR" ]
        lines!(ax, result.t / 365, result[:, lbl]; label = lbl)
    end 
    ax.xlabel = "Time, years"
    ax.ylabel = "Proportion"
end 

end # module MID_41
