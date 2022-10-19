
module MID_21
  
using CairoMakie, DataFrames, DifferentialEquations

export sir21!, run_sir21, print_sir21, plot_sir21

function sir21!(du, u, p, t) 
    # compartments 
    S, I, R = u
    # parameters
    beta, gamma = p 

    # the ODEs
    du[1] = -beta * S * I               # dS
    du[2] = beta * S * I - gamma * I    # dI
    du[3] = gamma * I                   # dR
end 

function run_sir21(u0, p, duration; saveat = 1)
    @assert minimum(u0) >= 0 "Input u0 = $u0: Cannot run model with negative starting values in any compartment"
    @assert sum(u0) == 1 "Input u0 = $u0: Compartment values are proportions so must sum to 1"
    @assert minimum(p) >= 0 "Input p = $p: Cannot run with negative values for any parameters"
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"

    tspan = ( 0., Float64(duration) )

    prob = ODEProblem(sir21!, u0, tspan, p)
    sol = solve(prob; saveat)

    return sol
end 

function dataframe_sir21(sol)
    result = DataFrame(t = Float64[], S = Float64[], I = Float64[], R = Float64[])
    for i ∈ eachindex(sol.u)
        push!( result, Dict(
            :t => sol.t[i], 
            :S => sol.u[i][1], 
            :I => sol.u[i][2], 
            :R => sol.u[i][3]
        ) )
    end 
    return result 
end 

function plot_sir21(result; kwargs...)
    fig = Figure()
    plot_sir21!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_sir21!(any, sol; kwargs...) = plot_sir21!(any, dataframe_sir21(sol); kwargs...)

function plot_sir21!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_sir21!(gl, result; kwargs...)
end 

function plot_sir21!(gl::GridLayout, result::DataFrame; legend = :right)
    ax = Axis(gl[1, 1])
    plot_sir21!(ax, result)
    Label(
        gl[0, :], 
        "p2.1.jl: Susceptible--infectious--resistant model with a constant population"
    )
    if legend == :right
        leg = Legend(gl[1, 2], ax)
    elseif legend == :below 
        leg = Legend(gl[2, 1], ax)
    elseif legend == :none 
        # no legend 
    else 
        @info "Unrecognised legend position given. Recognised options are `:right`, `:below` and `:none`"
    end
end 

function plot_sir21!(ax::Axis, result::DataFrame)
    xs = result.t ./ 7  # to plot time in weeks
    
    lines!(ax, xs, result.S, label = "Susceptible")
    lines!(ax, xs, result.I, label = "Infectious")
    lines!(ax, xs, result.R, label = "Recovered")
    ax.xlabel = "Time, weeks"
    ax.ylabel = "Fraction of population"
end 


## Help text 

"""
    run_sir21(u0, p, duration[; saveat])

Run the model `sir21!`, Susceptible--infectious--resistant model with a constant 
population

# Arguments 

`u0` is a vector of the initial proportions in each compartment. The vector should 
have three values, in the order `susceptible`, `infectious`, `recovered`. These 
values are proportions so should sum to 1.

`p` is a vector of the parameters of the model. The vector should have two values, 
in the order `beta` (infectiousness parameter), `gamma` (recovery rate).

`duration` is the time the model should run for (time units are interpretted as 
days).

# Optional keyword arguments 

* `saveat`: How frequently the model should save values. Default is 1 (day).

# Examples 
```
julia> u0 = [1 - 1e-6, 1e-6, 0]
3-element Vector{Float64}:
 0.999999
 1.0e-6
 0.0

julia> p = [520 / 365, 1 / 7]
2-element Vector{Float64}:
 1.4246575342465753
 0.14285714285714285

julia> run_sir21(u0, p, 10; saveat = 2)
retcode: Success
Interpolation: 1st order linear
t: 6-element Vector{Float64}:
  0.0
  2.0
  4.0
  6.0
  8.0
 10.0
u: 6-element Vector{Vector{Float64}}:
 [0.999999, 1.0e-6, 0.0]
 [0.9999857036962815, 1.2963010421217267e-5, 1.3332932973141789e-6]
 [0.9998132460495519, 0.0001681257322247459, 1.8628218223384836e-5]
 [0.9975816100019375, 0.0021756898099265973, 0.0002427001881360446]
 [0.9695595143428245, 0.027340832126607505, 0.0030996535305681827]
 [0.7144249430461295, 0.2518577279442249, 0.033717329009645705] 
```
"""
function run_sir21() end 

"""
    plot_sir21(result[; legend])

Plot the results of the model `sir21!`.

`result` can be either the ODE solver output or a DataFrame produced with `dataframe_sir21`.

`legend` indicates where the plot legend should be positioned relative to the plot. 
The recognised options are `:right`, `:below` and `:none`.
"""
function plot_sir21() end 

"""
    plot_sir21!(fig::Figure, result[; legend])
    plot_sir21!(gl::GridLayout, result[; legend])
    plot_sir21!(ax::Axis, result)

Plot the results of the model `sir21!`.

This function is provided to give more flexibility than only using `plot_sir21`. 
It allows the results of the model to be added to an existing figure. The function 
accepts a `Figure`, `GridLayout` or `Axis` to plot on to.

`result` can be either the ODE solver output or a DataFrame produced with `dataframe_sir21`.

`legend` indicates where the plot legend should be positioned relative to the plot. 
The recognised options are `:right`, `:below` and `:none`. Note that this option 
is not available if an `Axis` was provided to the function.
"""
function plot_sir21!() end 

end # module MID_21
