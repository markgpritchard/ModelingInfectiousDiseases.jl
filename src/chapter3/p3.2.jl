
module MID_32
  
# SIS model with multiple risk groups (page 64)

using CairoMakie, DataFrames, DifferentialEquations
import Base: minimum

export Parameters32, sis32!, run_sis32, dataframe_sis32, plot_sis32, plot_sis32!

struct Parameters32
    beta    :: Matrix{<:Float64} 
    gamma   :: Vector{<:Float64}
end

function sis32!(du, u, p, t) 
    # parameters 
    β = p.beta 
    γ = p.gamma
    
    # the ODEs
    for i ∈ axes(u, 2) # the number of risk groups
        infections = sum( [ β[i, j] * u[1, i] * u[2, j] for j ∈ axes(u, 2) ] )     
                                                # = Σ βᵢⱼ Sᵢ Iⱼ
        recoveries = γ[i] * u[2, i]             # = γᵢ Iᵢ

        du[1, i] = -infections + recoveries     # dSᵢ 
        du[2, i] = infections - recoveries      # dIᵢ 
    end 
end 

function run_sis32(; S0, I0, betavector = nothing, betaconstant = 1, beta = nothing, gamma, duration, kwargs...)
    @assert length(S0) == length(I0) "Must have equal length vectors for S0 and I0"
    u0 = zeros(2, length(S0))
    u0[1, :] = S0 
    u0[2, :] = I0 
    # This function allows you to provide a matrix for `beta`, or a vector and constant 
    # (`betavector` and `betaconstant`) to calculate `beta`. If calculating `beta` 
    # we use an intermediary step when creating the transmission factor so that 
    # all βᵢⱼ == βⱼᵢ 
    if isnothing(beta) beta = betaconstant .* betavector * betavector' end 
    p = Parameters32(beta, gamma)
    run_sis32(u0, p, duration; kwargs...)
end

function run_sis32(u0, p, duration; saveat = 1)
    @assert minimum(u0) >= 0 "Input u0 = $u0: no compartments can contain negative proportions"
    @assert sum(u0) ≈ 1 "Input u0 = $u0: compartment values are proportions and must approximately sum to 1"
    @assert minimum(p) >= 0 "Input p = $p: cannot run with negative parameters" 
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"

    tspan = ( 0., Float64(duration) )

    prob = ODEProblem(sis32!, u0, tspan, p)
    sol = solve(prob; saveat)
    return sol
end 

minimum(p::Parameters32) = min( minimum(p.beta), minimum(p.gamma) )

function dataframe_sis32(sol; type = :both)
    result = DataFrame(t = sol.t)
    if type == :both    # susceptible and infectious to be included in DataFrame 
        dataframe_sis32!(result, sol, :S)
        dataframe_sis32!(result, sol, :I)
    else 
        dataframe_sis32!(result, sol, type)
    end 
    return result 
end 

function dataframe_sis32!(result, sol, type)
    row = dataframerow(type)
    if isnothing(row) return end    # for when dataframerow returns an error
    categories = axes(sol[1], 2)
    for i ∈ categories
        rs = zeros( size(sol, 3) )
        for j ∈ axes(sol, 3) 
            rs[j] = sol[j][row, i] 
        end 
        insertcols!(result, Symbol("$type$i") => rs)
    end 
end 

function dataframerow(type)
    if type == :S 
        return 1 
    elseif type == :I 
        return 2 
    else 
        @error "Unrecognised type entered. Recognised options are :S, :I and :both"
    end 
end 

function plot_sis32(result; kwargs...)
    fig = Figure()
    plot_sis32!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_sis32!(any, sol; kwargs...) = plot_sis32!(any, dataframe_sis32(sol); kwargs...)

function plot_sis32!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_sis32!(gl, result; kwargs...)
end 

function plot_sis32!(gl::GridLayout, result::DataFrame; 
        label = "p3.2.jl: Susceptible--infectious--susceptible model with multiple risk groups", 
        legend = :right
    )
    # Plot on a linear scale 
    ax1 = Axis(gl[1, 1])
    plot_sis32!(ax1, result)
    ax1.ylabel = "Fraction of population\ninfectious (linear scale)"

    # Plot on a log scale
    ax2 = Axis(gl[2, 1], yscale = log10)
    plot_sis32!(ax2, result, plotvalue_nozero)
    ax2.ylabel = "Fraction of population\ninfectious (log scale)"
    ax2.xlabel = "Time, days"

    # And plot proportions within risk group instead of within population 
    # Plot on a linear scale 
    ax3 = Axis(gl[1, 2])
    plot_sis32!(ax3, result, plotvalue_dividen)
    ax3.ylabel = "Fraction of risk group\ninfectious (linear scale)"

    # Plot on a log scale
    ax4 = Axis(gl[2, 2], yscale = log10)
    plot_sis32!(ax4, result, plotvalue_dividen_nozero)
    ax4.ylabel = "Fraction of risk group\ninfectious (log scale)"
    ax4.xlabel = "Time, days"

    linkxaxes!(ax1, ax2); hidexdecorations!(ax1; grid = false, ticks = false)
    linkxaxes!(ax3, ax4); hidexdecorations!(ax3; grid = false, ticks = false)

    Label(gl[0, :], label)

    if legend == :right
        leg = Legend(gl[1:2, 3], ax1)
    elseif legend == :below 
        leg = Legend(gl[3, 1:2], ax1; orientation = :horizontal)
    elseif legend == :none 
        # no legend 
    else 
        @info "Unrecognised legend position given. Recognised options are `:right`, `:below` and `:none`"
    end
end 

function plot_sis32!(ax::Axis, result::DataFrame, f = plotvalue)
    for col ∈ names(result)
        col[1] != 'I' && continue 
        lines!(ax, result.t, f(result, col), label = "$col")
    end 
end 


## Functions to choose which values to plot 

plotvalue(result, col) = result[!, col]
function plotvalue_dividen(result, col) # proportion within group rather than within population
    scol = "S$(col[2:end])"
    n = result[!, col] .+ result[!, scol]
    vals =  result[!, col] ./ n 
    return vals 
end 

# for plots on log scale, need to omit any zeros from the plotting
value_nozero(vals) = [ ifelse(v == 0, NaN, v) for v ∈ vals ]
plotvalue_nozero(result, col) = value_nozero(plotvalue(result, col))
plotvalue_dividen_nozero(result, col) = value_nozero(plotvalue_dividen(result, col))

end # module MID_31
