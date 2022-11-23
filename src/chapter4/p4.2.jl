
module MID_42
  
using CairoMakie, DataFrames, DifferentialEquations
import Base: minimum

export Parameters42, spr42!, run_spr42, plot_spr42, plot_spr42!

struct Parameters42
    β   :: Float64
    γ   :: Float64
    μ   :: Float64
    a   :: Float64
    n   :: Int 
    c   :: Matrix{Int64}
end

const compartmentsymbols = Dict(:S => 1, :P => 2, :R => 3, :λ => 4)

minimum(p::Parameters42) = min( p.β, p.γ, p.μ, p.a, p.n, minimum(p.c) )

function spr42!(du, u, p, t) 
    # compartments 
    S = u[1, :] 
    P = u[2, :] 
    R = u[3, :] 
    λ = u[4, :]     

    for i ∈ 1:p.n
        du[1, i] = p.μ * (1 - S[i]) - S[i] * sum([ p.c[i, j] * λ[j] for j ∈ 1:p.n ])                        # dS
        du[2, i] = S[i] * sum([ ifelse(i == j, 0, p.c[i, j]) * λ[j] for j ∈ 1:p.n ]) - P[i] * (λ[i] + p.μ)  # dP
        du[3, i] = (S[i] + P[i]) * λ[i] - p.μ * R[i]                                                        # dR
        du[4, i] = (p.β * (S[i] + p.a * P[i]) - p.γ - p.μ) * λ[i]                                           # dλ
    end 
end 

function run_spr42(u0, p::Parameters42, duration; saveat = 1)
    @assert minimum(u0) >= 0 "Input u0 = $u0: cannot run with negative compartment values"
    @assert minimum(p) >= 0 "Input p = $p: cannot run with negative parameters"
    @assert duration > 0 "Input duration = $duration: cannot run with zero or negative duration"
    for i ∈ axes(u0, 2)
        if !isapprox(sum(u0[1:3, i]), 1)
            @error "Input u0 = $u0: compartments are the sum of S, R, R in each column must equal 1"
        end 
    end 

    tspan = ( 0., Float64(duration) )
    prob = ODEProblem(spr42!, u0, tspan, p)
    sol = solve(prob; saveat)
    return sol
end 

function run_spr42(u0, p::Vector{<:Real}, duration; kwargs...)
    n = Int(p[5])
    c = zeros(Int, n, n)
    for i ∈ 1:n, j ∈ 1:n 
        if i == j || i == j + 1 || i == j - 1
            c[i, j] = 1 
        end
        if i == 1 && j == n c[i, j] = 1 end
        if i == n && j == 1 c[i, j] = 1 end
    end 
    parms = Parameters42(
        p[1],   # β
        p[2],   # γ
        p[3],   # μ
        p[4],   # a
        n,      # n
        c       # c
    )

    return run_spr42(u0, parms, duration; kwargs...)
end 

function run_spr42(; n, S0, P0, R0 = nothing, lambda0, beta, gamma, mu, a, duration, kwargs...)
    if isnothing(R0) R0 = [ 1 - (S0[i] + P0[i]) for i ∈ eachindex(S0) ] end 
    @assert length(S0) == length(P0) == length(R0) == length(lambda0) == n "Must have 
    equal length vectors for `S0`, `P0`, `R0` and `lambda0`, and must be equal to `n`"

    u0 = zeros(4, length(S0))
    u0[1, :] = S0
    u0[2, :] = P0
    u0[3, :] = R0
    u0[4, :] = lambda0
    p = [beta, gamma, mu, a, n]
    return run_spr42(u0, p, duration; kwargs...)
end 

function dataframe_spr42(sol, compartment::Symbol)
    compartmentnumber = compartmentsymbols[compartment]
    compartmentlabel = String(compartment)

    result = DataFrame(t = sol.t)
    for i ∈ axes(sol, 2)
        insertcols!(result, Symbol("$compartmentlabel$i") => [ sol[j][compartmentnumber, i] for j ∈ axes(sol, 3) ])
    end 

    return result 
end 

function plot_spr42(sol; kwargs...)
    fig = Figure()
    plot_spr42!(fig, sol; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

function plot_spr42!(fig::Figure, sol; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_spr42!(gl, sol; kwargs...)
end 

function plot_spr42!(gl::GridLayout, sol; 
        label = "p4.2.jl: SIR model with multiple strains giving partial immunity", kwargs...
    )
    axs = [ Axis(gl[i, 1]) for i ∈ 1:4 ]
    for (i, c) ∈ enumerate([:S, :P, :R, :λ])
        data = dataframe_spr42(sol, c)
        plot_spr42!(axs[i], data)#; kwargs...)
        leg = Legend(gl[i, 2], axs[i])
        i < 4 && hidexdecorations!(axs[i]; ticks = false, grid = false)
        if i <= 3 
            axs[i].ylabel = "Proportion $(String(c))" 
        else 
            axs[i].ylabel = "$(String(c))"
        end 
    end 
    axs[4].xlabel = "Time"
    Label(gl[0, 1], label)
end 

function plot_spr42!(ax::Axis, data)
    xs = data.t
    select!(data, Not(:t))
    for i ∈ axes(data, 2) 
        lines!(ax, xs, data[:, i], label = names(data)[i])
    end 
end 

end # module MID_42
