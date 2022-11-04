
module MID_33
  
using CairoMakie, DataFrames, DifferentialEquations
import Base: minimum

export Parameters33, sir33!, run_sir33, dataframe_sir33, plot_sir33, plot_sir33!

struct Parameters33
    beta    :: Matrix{<:Float64} 
    gamma   :: Float64
    lambda  :: Float64
    mu      :: Vector{<:Float64}
    nu      :: Float64
end

minimum(p::Parameters33) = min( minimum(p.beta), p.gamma, p.lambda, minimum(p.mu), p.nu )

function sir33!(du, u, p, t) 
    # compartments 
    (Sc, Ic, Rc, Sa, Ia, Ra) = u

    # parameters 
    (βcc, βca, βac, βaa) = p.beta 
    γ = p.gamma
    λ = p.lambda 
    μc, μa = p.mu 
    nu = p.nu 
    
    # events 
    infection_c = Sc * ( βcc * Ic + βca * Ia )
    recovery_c = γ * Ic 
    infection_a = Sa * ( βac * Ic + βaa * Ia )
    recovery_a = γ * Ia 

    # the ODEs
    du[1, 1] = nu - infection_c - (μc + λ) * Sc             # dSc
    du[2, 1] = infection_c - recovery_c - (μc + λ) * Ic     # dIc
    du[3, 1] = recovery_c - (μc + λ) * Rc                   # dRc
    du[1, 2] = λ * Sc - infection_a - μa * Sa               # dSa
    du[2, 2] = λ * Ic + infection_a - recovery_a - μa * Ia  # dIa
    du[3, 2] = λ * Rc + recovery_a - μa * Ra                # dRa
end 

function run_sir33(u0, p, duration; saveat = 1)
    @assert minimum(u0) >= 0 "Input u0 = $u0: no compartments can contain negative proportions"
    @assert sum(u0) ≈ 1 "Input u0 = $u0: compartment values are proportions and must sum to 1"
    @assert minimum(p) >= 0 "Input p = $p: cannot run with negative parameters" 
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"

    tspan = ( 0., Float64(duration) )

    prob = ODEProblem(sir33!, u0, tspan, p)
    sol = solve(prob; saveat)

    return sol
end 

function dataframe_sir33(sol)
    result = DataFrame(t = sol.t)
    _dataframe_sir33!(result, sol, :c, 1)
    _dataframe_sir33!(result, sol, :a, 2)
    return result 
end 

function _dataframe_sir33!(result, sol, type, row)
    categories = [ "S", "I", "R" ]
    for (i, cat) ∈ enumerate(categories)
        rs = zeros( size(sol, 3) )
        for j ∈ axes(sol, 3) 
            rs[j] = sol[j][i, row] 
        end 
        insertcols!(result, Symbol("$cat$type") => rs)
    end 
end 

function plot_sir33(result; kwargs...)
    fig = Figure()
    plot_sir33!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_sir33!(any, sol; kwargs...) = plot_sir33!(any, dataframe_sir33(sol); kwargs...)

function plot_sir33!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_sir33!(gl, result; kwargs...)
end 

function plot_sir33!(gl::GridLayout, result::DataFrame; 
        label = "p3.3.jl: Susceptible--infectious--recovered model with two age groups", 
        legend = :right
    )

    axs = [ Axis(gl[i, 1]) for i ∈ 1:2 ]
    plot_sir33!(axs[1], select(result, :t, :Sc, :Ic, :Rc))
    axs[1].ylabel = "Fraction children"

    plot_sir33!(axs[2], select(result, :t, :Sa, :Ia, :Ra))
    axs[2].ylabel = "Fraction adults"

    linkxaxes!(axs...)
    hidexdecorations!(axs[1]; grid = false, ticks = false)
    axs[2].xlabel = "Time"
    Label(gl[0, :], label)

    if legend == :right
        leg = Legend(gl[1:2, 2], axs[1])
    elseif legend == :below 
        leg = Legend(gl[3, 1], axs[1]; orientation = :horizontal)
    elseif legend == :none 
        # no legend 
    else 
        @info "Unrecognised legend position given. Recognised options are `:right`, `:below` and `:none`"
    end
end 

function plot_sir33!(ax::Axis, result::DataFrame)
    lbls = [ "Susceptible", "Infectious", "Recovered"] 
    for (i, lbl) ∈ enumerate(lbls)
        lines!(ax, result.t, result[!, names(result)[i+1]]; label = lbl)
    end 
end 


## Help text

"""
    run_sir33(u0, p, duration[; saveat])

Run the model `sir33!`, a susceptible--infectious--recovered model with two age  
groups.

## Arguments 

`u0` is a `3x2` matrix of starting conditions for the model. There is a column for 
each age group with values from top to bottom in the order susceptible, infectious, 
recovered. 

`p` is the parameters for the model. This is expected in a `Parameters33` type. 
`p.beta` is an `2x2` matrix of infectiousness parameters between each age group. 
`p.gamma` is the recovery rate. `p.lambda` is the rate at which children become 
adults. `p.mu` is a vector of mortality rates. `p.nu` is the birth rate.

`duration` is the length of time that the model should run.

## Optional keyword arguments 
* `saveat`: how frequently the ODE solver should save results. Default is `1`

## Example 
```
julia> u0 = [.1 .1; .0001 .0001; .0999 .6999]
3×2 Matrix{Float64}:
 0.1     0.1
 0.0001  0.0001
 0.0999  0.6999

julia> p = Parameters33( [100. 10.; 10. 20.], 10., 1/15, 
[0., 1/60], 1/60 )
Parameters33([100.0 10.0; 10.0 20.0], 10.0, 0.06666666666666667, [0.0, 0.016666666666666666], 0.016666666666666666)

julia> run_sir33(u0, p, 8; saveat = 2)
retcode: Success
Interpolation: 1st order linear
t: 5-element Vector{Float64}:
 0.0
 2.0
 4.0
 6.0
 8.0
u: 5-element Vector{Matrix{Float64}}:
 [0.1 0.1; 0.0001 0.0001; 0.0999 0.6999]
 [0.11369619644458144 0.11022606092659487; 0.0006128138490162058 7.716636864826869e-5; 0.09193232375421641 0.6901173671791364]
 [0.09148067115054088 0.1147344368340036; 0.002279884006076226 0.00038815171076258696; 0.11794302792518478 0.6864707713377334]
 [0.09430041595582973 0.11967324985168493; 0.0003866986851766598 6.839905516128923e-5; 0.1217968830571779 0.6836563760193186]
 [0.10631969678297161 0.1277247071891025; 0.0005598435665349642 9.138146304083705e-5; 0.11378814867496992 0.6779141075704351]
```
"""
function run_sir33() end

end # module MID_33
