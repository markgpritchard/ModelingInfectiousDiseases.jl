
module MID_31
  
using CairoMakie, DataFrames, DifferentialEquations
import Base: minimum

export Parameters31, sir31!, run_sir31, dataframe_sir31, plot_sir31, plot_sir31!

struct Parameters31 
    beta    :: Matrix{<:Float64} 
    gamma   :: Float64
end

minimum(p::Parameters31) = min(minimum(p.beta), p.gamma)

function sir31!(du, u, p, t) 
    # compartments 
    Sh, Ih, Sℓ, Iℓ = u

    # parameters 
    (βhh, βhl, βlh, βll) = p.beta 
    γ = p.gamma
    
    # the ODEs
    du[1] = -(βhh * Ih + βhl * Iℓ) * Sh + γ * Ih    # dSh
    du[2] = (βhh * Ih + βhl * Iℓ) * Sh - γ * Ih     # dIh
    du[3] = -(βlh * Ih + βll * Iℓ) * Sℓ + γ * Iℓ    # dSℓ
    du[4] = (βlh * Ih + βll * Iℓ) * Sℓ - γ * Iℓ     # dIℓ
end 

function run_sir31(u0, p, duration; saveat = 1)
    @assert minimum(u0) >= 0 "Input u0 = $u0: Cannot run model with negative starting values in any compartment"
    @assert sum(u0) <= 1 "Input u0 = $u0: Compartment values are proportions so must sum to 1 or less"
    @assert minimum(p) >= 0 "Input p = $p: Cannot run with negative values for any parameters"
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"

    tspan = ( 0., Float64(duration) )

    prob = ODEProblem(sir31!, u0, tspan, p)
    sol = solve(prob; saveat)

    return sol
end 

function run_sir31(; Nh, Nl = 1 - Nh, Ih0, Il0, Sh0 = Nh - Ih0, Sl0 = Nl - Il0, 
        beta_hh = nothing, beta_hl = nothing, beta_lh = nothing, beta_ll = nothing, 
        beta = [beta_hh beta_hl; beta_lh beta_ll], gamma, duration, kwargs...
    )
    @assert isa(beta, Matrix{Float64}) "Either enter `beta` as a matrix of Float64, 
    or enter a Float64 for each of `beta_hh`, `beta_hl`, `beta_lh`, `beta_ll`"

    u0 = [Sh0, Ih0, Sl0, Il0]
    p = Parameters31(beta, gamma)
    return run_sir31(u0, p, duration; kwargs...)
end 

function dataframe_sir31(sol)
    result = DataFrame(t = Float64[], Sh = Float64[], Ih = Float64[], Sl = Float64[], Il = Float64[])
    for i ∈ eachindex(sol.u)
        push!( result, Dict(
            :t => sol.t[i], 
            :Sh => sol.u[i][1], 
            :Ih => sol.u[i][2], 
            :Sl => sol.u[i][3],
            :Il => sol.u[i][4]
        ) )
    end 
    return result 
end 

function plot_sir31(result; kwargs...)
    fig = Figure()
    plot_sir31!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_sir31!(any, sol; kwargs...) = plot_sir31!(any, dataframe_sir31(sol); kwargs...)

function plot_sir31!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_sir31!(gl, result; kwargs...)
end 

function plot_sir31!(gl::GridLayout, result::DataFrame; 
        label = "p3.1.jl: Susceptible--infectious--resistant model with a constant population", 
        legend = :right, kwargs...
    )
    ax1 = Axis(gl[1, 1])
    plot_sir31!(ax1::Axis, result::DataFrame; ylabel = "Fraction infectious\n(linear scale)")
    ax2 = Axis(gl[2, 1], yscale = log10)
    plot_sir31!(ax1::Axis, result::DataFrame; ylabel = "Fraction infectious\n(log scale)")
    ax2.xlabel = "Time, days"
    hidexdecorations!(ax1; grid = false, ticks = false)
    linkxaxes!(ax1, ax2)

    Label(gl[0, :], label)

    if legend == :right
        leg = Legend(gl[1:2, 2], ax)
    elseif legend == :below 
        leg = Legend(gl[3, 1], ax)
    elseif legend == :none 
        # no legend 
    else 
        @info "Unrecognised legend position given. Recognised options are `:right`, `:below` and `:none`"
    end
end 

function plot_sir31!(ax::Axis, result::DataFrame; ylabel = "")
    lines!(ax1, xs, result.Ih, label = "High risk")
    lines!(ax1, xs, result.Il, label = "Low risk")
    ax.ylabel = ylabel 
end 

end # module MID_31
