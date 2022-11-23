
module MID_43
  
using CairoMakie, DataFrames, DifferentialEquations
import Base: minimum

export Parameters43, seicr43!, run_seicr43, dataframe_seicr43, plot_seicr43, plot_seicr43!

struct Parameters43
    alpha   :: Vector{<:Float64}
    beta    :: Vector{<:Float64}
    gamma   :: Vector{<:Float64}
    delta   :: Vector{<:Float64}
    mu      :: Float64
    nu      :: Float64
    xi      :: Vector{<:Float64}
    rho     :: Vector{<:Float64}
    sigma   :: Vector{<:Float64}
    phi     :: Vector{<:Float64}
    psi     :: Vector{<:Float64}
end

minimum(p::Parameters43) = min( 
    minimum(p.alpha), 
    minimum(p.beta), 
    minimum(p.gamma), 
    minimum(p.delta), 
    p.mu, 
    p.nu, 
    minimum(p.xi), 
    minimum(p.rho), 
    minimum(p.sigma), 
    minimum(p.phi), 
    minimum(p.psi) 
)

function seicr43!(du, u, p, t) 
    # compartments (including parameter results that are integrated)
    S, E1, E2, I1, I2, C1, C2, R1, R2, R12, ε1, ε2, λ1, λ2 = u  
    
    # ODEs 
    du[1] = p.nu - (λ1 + λ2 + p.mu) * S                                                         # dS 
    du[2] = λ1 * S - (p.phi[2] * λ2 + p.sigma[1] + p.mu) * E1                                   # dE1 
    du[3] = λ2 * S - (p.phi[1] * λ1 + p.sigma[2] + p.mu) * E2                                   # dE2 
    du[4] = p.sigma[1] * E1 - (p.phi[2] * λ2 + p.gamma[1] + p.mu) * I1                          # dI1 
    du[5] = p.sigma[2] * E2 - (p.phi[1] * λ1 + p.gamma[2] + p.mu) * I2                          # dI2 
    du[6] = p.gamma[1] * I1 - (p.xi[2] * p.phi[2] * λ2 + p.delta[1] + p.mu) * C1                # dC1 
    du[7] = p.gamma[2] * I2 - (p.xi[1] * p.phi[1] * λ1 + p.delta[2] + p.mu) * C2                # dC2 
    du[8] = (1 - p.rho[1]) * p.delta[1] * C1 - (p.alpha[2] * λ2 + p.mu) * R1                    # dR1 
    du[9] = (1 - p.rho[2]) * p.delta[2] * C1 - (p.alpha[1] * λ1 + p.mu) * R2                    # dR2 
    du[10] = (1 - p.rho[1]) * (1 - p.rho[2]) * (λ2 * p.phi[2] * E1 + p.phi[2] * I1 + 
        p.xi[2] * p.phi[2] * C1 + λ1 * p.phi[1] * E2 + p.phi[1] * I2 + p.xi[1] * p.phi[1] * C2) + 
        (1 - p.psi[2] * p.rho[2]) * p.alpha[2] * λ2 * R1 + 
        (1 - p.psi[1] * p.rho[1]) * p.alpha[1] * λ1 * R2 - p.mu * R12                           # dR12 
    du[11] = λ1 * S + p.phi[1] * λ1 * (E2 + I2 + p.xi[1] * C2) + p.alpha[1] * λ1 * R2 - 
        (p.sigma[1] + p.mu) * ε1                                                                # dε1
    du[12] = λ2 * S + p.phi[2] * λ2 * (E1 + I1 + p.xi[2] * C1) + p.alpha[2] * λ2 * R1 - 
        (p.sigma[2] + p.mu) * ε2                                                                # dε2
    du[13] = p.beta[1] * p.sigma[1] * ε1 - (p.gamma[1] + p.mu) * λ1                             # dλ1
    du[14] = p.beta[2] * p.sigma[2] * ε2 - (p.gamma[2] + p.mu) * λ2                             # dλ2
end 

function run_seicr43(u0, p, duration; saveat = 1)
    @assert minimum(u0) >= 0 "Input u0 = $u0: cannot run with negative compartment values"
    @assert minimum(p) >= 0 "Input p = $p: cannot run with negative parameters"
    @assert maximum(p.rho) <= 1 "Input p.rho = $(p.rho): is a probability and cannot be > 1"
    @assert maximum(p.phi) <= 1 "Input p.phi = $(p.phi): is a probability and cannot be > 1"
    @assert maximum(p.psi) <= 1 "Input p.psi = $(p.psi): is a probability and cannot be > 1"
    @assert duration > 0 "Input duration = $duration: cannot run with zero or negative duration"

    tspan = ( 0., Float64(duration) )
    prob = ODEProblem(seicr43!, u0, tspan, p)
    sol = solve(prob; saveat)
    return sol
end 

function run_seicr43(; S_0, E1_0, E2_0, I1_0, I2_0, C1_0, C2_0, R1_0, R2_0, R12_0, 
        ε1_0, ε2_0, λ1_0, λ2_0, alpha, beta, gamma, delta, mu, nu, xi, rho, sigma, 
        phi, psi, duration, kwargs...
    )
    u0 = [S_0, E1_0, E2_0, I1_0, I2_0, C1_0, C2_0, R1_0, R2_0, R12_0, ε1_0, ε2_0, λ1_0, λ2_0]
    p = Parameters43(alpha, beta, gamma, delta, mu, nu, xi, rho, sigma, phi, psi)
    return run_seicr43(u0, p, duration; kwargs...)
end 

dataframe_seicr43(sol, p::Parameters43) = dataframe_seicr43(sol, p.beta)

function dataframe_seicr43(sol, beta::Vector{<:Real})
    result = DataFrame(
        t = Float64[], 
        S1 = Float64[], S2 = Float64[], 
        E1 = Float64[], E2 = Float64[], 
        I1 = Float64[], I2 = Float64[]
    )
    for i ∈ eachindex(sol.u)
        push!( result, Dict(
            :t => sol.t[i], 
            :S1 => sol.u[i][1] + sol.u[i][3] + sol.u[i][5] + sol.u[i][7] + sol.u[i][9], 
            :S2 => sol.u[i][1] + sol.u[i][2] + sol.u[i][4] + sol.u[i][6] + sol.u[i][8], 
            :E1 => sol.u[i][11],
            :E2 => sol.u[i][12],
            :I1 => sol.u[i][13] / beta[1],
            :I2 => sol.u[i][14] / beta[2],
        ) )
    end 
    return result 
end 

function plot_seicr43(sol; kwargs...)
    fig = Figure()
    plot_seicr43!(fig, sol; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_seicr43!(any, sol, beta; kwargs...) = 
    plot_seicr43!(any, dataframe_seicr43(sol, beta); kwargs...)

function plot_seicr43!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_seicr43!(gl, result; kwargs...)
end 

function plot_seicr43!(gl::GridLayout, result::DataFrame; 
        label = "p4.2.jl: SIR model with multiple strains giving partial immunity", kwargs...
    )
    axs = [ Axis(gl[i, 1]) for i ∈ 1:3 ]
    for (i, c) ∈ enumerate([:S, :E, :I])
        plot_seicr43!(axs[i], result, c)
        leg = Legend(gl[i, 2], axs[i])
        i < 3 && hidexdecorations!(axs[i]; ticks = false, grid = false)
        axs[i].ylabel = "Proportion $(String(c))" 
    end 
    axs[3].xlabel = "Time"
    Label(gl[0, 1], label)
end 

function plot_seicr43!(ax::Axis, result::DataFrame, compartments)
    xs = result.t
    for i ∈ 1:2
        lines!(ax, xs, result[:, "$compartments$i"], label = "$compartments$i")
    end 
end 

end # module MID_43
