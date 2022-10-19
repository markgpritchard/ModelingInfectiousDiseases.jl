
module MID_32
  
using CairoMakie, DataFrames, DifferentialEquations
import Base: minimum

export Parameters32, sis32!, run_sis32, dataframe_sis32, plot_sis32, plot_sis32!

struct Parameters32
    beta    :: Matrix{<:Float64} 
    gamma   :: Vector{<:Float64}
end

minimum(p::Parameters32) = min( minimum(p.beta), minimum(p.gamma) )

function sis32!(du, u, p, t) 
    # parameters 
    β = p.beta 
    γ = p.gamma
    
    # the ODEs
    for i ∈ axes(u0, 2) # the number of risk groups
        infections = sum( [ β[i, j] * u[1, i] * u[2, j] for j ∈ axes(u0, 2) ] )     
                                        # = Σ βᵢⱼ Sᵢ Iⱼ
        recoveries = γ[i] * u[2, i]     # γᵢ Iᵢ

        du[1, i] = -infections + recoveries     # dSᵢ 
        du[2, i] = infections - recoveries      # dIᵢ 
    end 
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

function dataframe_sis32(sol; type = :both)
    result = DataFrame(t = sol.t)
    if type == :both    # susceptible and infectious to be included in DataFrame 
        _dataframe_sis32!(result, sol, :S)
        _dataframe_sis32!(result, sol, :I)
    else 
        _dataframe_sis32!(result, sol, type)
    end 
    return result 
end 

function _dataframe_sis32!(result, sol, type)
    row = _dataframerow(type)
    if isnothing(row) return end    # for when _dataframerow returns an error
    categories = axes(sol[1], 2)
    for i ∈ categories
        rs = zeros( size(sol, 3) )
        for j ∈ axes(sol, 3) 
            rs[j] = sol[j][row, i] 
        end 
        insertcols!(result, Symbol("$type$i") => rs)
    end 
end 

function _dataframerow(type)
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
    ax1.ylabel = "Fraction of population infectious\n(linear scale)"

    # Plot on a log scale
    ax2 = Axis(gl[2, 1], yscale = log10)
    plot_sis32!(ax2, result, plotvalue_nozero)
    ax2.ylabel = "Fraction of population infectious\n(log scale)"
    ax2.xlabel = "Time, days"

    # And plot proportions within risk group instead of within population 
    # Plot on a linear scale 
    ax3 = Axis(gl[1, 2])
    plot_sis32!(ax3, result, plotvalue_dividen)
    ax3.ylabel = "Fraction of risk group infectious\n(linear scale)"

    # Plot on a log scale
    ax4 = Axis(gl[2, 2], yscale = log10)
    plot_sis32!(ax4, result, plotvalue_dividen_nozero)
    ax4.ylabel = "Fraction of risk group infectious\n(log scale)"
    ax4.xlabel = "Time, days"


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


## A set of functions to choose which values to plot 

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


## Help text

"""
    run_sis32(u0, p, duration[; saveat])

Run the model `sis32!`, a susceptible--infectious--susceptible model with multiple 
risk group categories.

## Arguments 

`u0` is a `2xn` matrix of starting conditions for the model. There is a column for 
each risk group. The upper value is the proportion of the whole population that 
is susceptible and in that risk group, and the lower value is the proportion of 
the whole population that is infectious and in that risk group. 

`p` is the parameters for the model. This is expected in a `Parameters32` type. 
`p.beta` is an `nxn` matrix of infectiousness parameters between each risk group. 
`p.gamma` is an n-element vector of recovery rates for each risk group.

`duration` is the length of time that the model should run.

## Optional keyword arguments 
* `saveat`: how frequently the ODE solver should save results. Default is `1`

## Example 
```
julia> u0 = [ .06 .31 .52 .08 .02999; .0 .0 .0 .0 1e-5 ]    
2×5 Matrix{Float64}:
 0.06  0.31  0.52  0.08  0.02999
 0.0   0.0   0.0   0.0   1.0e-5

julia> betavector = [ 0, 3, 10, 60, 100 ]
5-element Vector{Int64}:
   0
   3
  10
  60
 100

julia> betamatrix = .0016 .* betavector * betavector'       
5×5 Matrix{Float64}:
 0.0  0.0     0.0    0.0     0.0
 0.0  0.0144  0.048  0.288   0.48
 0.0  0.048   0.16   0.96    1.6
 0.0  0.288   0.96   5.76    9.6
 0.0  0.48    1.6    9.6    16.0

julia> p = Parameters32( betamatrix, [.2, .2, .2, .2, .2] ) 
Parameters32([0.0 0.0 … 0.0 0.0; 0.0 0.0144 … 
0.288 0.48; … ; 0.0 0.288 … 5.76 9.6; 0.0 0.48 … 9.6 16.0], [0.2, 0.2, 0.2, 0.2, 0.2])

julia> run_sis32(u0, p, 8; saveat = 2)
retcode: Success
Interpolation: 1st order linear
t: 5-element Vector{Float64}:
 0.0
 2.0
 4.0
 6.0
 8.0
u: 5-element Vector{Matrix{Float64}}:
 [0.06 0.31 … 0.08 0.02999; 0.0 0.0 … 0.0 1.0e-5]
 [0.06 0.30999 … 0.07997 0.02997; 0.0 6.6124e-6 … 3.4121e-5 2.802e-5]
 [0.06 0.30996 … 0.07980 0.02987; 0.0 3.9027e-5 … 0.0002 0.00013]
 [0.06 0.30979 … 0.07895 0.02934; 0.0 0.00021 … 0.00105 0.00066]
 [0.06 0.30897 … 0.07487 0.02687; 0.0 0.00102 … 0.00513 0.00313]
```
"""
function run_sis32() end

end # module MID_31
