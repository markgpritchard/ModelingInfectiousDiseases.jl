
module MID_34
  
using CairoMakie, DataFrames, DifferentialEquations
import Base: minimum

export Parameters34, seir34!, run_seir34, plot_seir34, plot_seir34!

struct Parameters34
    beta    :: Matrix{<:Float64} 
    sigma   :: Float64
    gamma   :: Float64
    mu      :: Float64  
    nu      :: Float64
end

minimum(p::Parameters34) = min( minimum(p.beta), p.sigma, p.gamma, p.mu, p.nu )

function seir34!(du, u, p, t)  # model runs for one year at a time
    # compartments 
    S = u[1, :]
    E = u[2, :]
    I = u[3, :]
    R = u[4, :]

    # parameters 
    β = p.beta 
    σ = p.sigma
    γ = p.gamma
    μ = zeros(4); μ[4] = p.mu   # make it a vector with a mu only for the adult group
    nu = p.nu 

    births = zeros(4); births[1] = nu * sum(u[:, 4])   # births only into youngest group

    for i ∈ 1:4     # 4 age groups
        du[1, i] = births[i] - sum([ β[i, j] * I[j] * S[i] for j ∈ 1:4 ]) - μ[i] * S[i] # dSi
        du[2, i] = sum([ β[i, j] * I[j] * S[i] for j ∈ 1:4 ]) - σ * E[i] - μ[i] * E[i]  # dEi 
        du[3, i] = σ * E[i] - γ * I[i] - μ[i] * I[i]                                    # dIi 
        du[4, i] = γ * I[i] - μ[i] * R[i]                                               # dRi
    end 
end 

function run_seir34(u0, p, duration; saveat = 1)
    @assert minimum(u0) >= 0 "Input u0 = $u0: no compartments can contain negative proportions"
    @assert sum(u0) ≈ 1 "Input u0 = $u0: compartment values are proportions and must sum to 1"
    @assert minimum(p) >= 0 "Input p = $p: cannot run with negative parameters" 
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"

    # initial u 
    u = u0

    # dataframe for results
    result = DataFrame(t = Float64[]) 
    columnnames = ["S1", "E1", "I1", "R1", "S2", "E2", "I2", "R2", "S3", "E3", 
        "I3", "R3", "S4", "E4", "I4", "R4"]
    for col ∈ columnnames 
        insertcols!(result, Symbol(col) => Float64[])
    end 

    τs = 1. # starting time

    while τs < duration    # for each year
        τe = min(duration, τs + 365)

        # run the ODE 
        tspan = (τs, τe)
        prob = ODEProblem(seir34!, u0, tspan, p)
        sol = solve(prob; saveat)

        # put these results into the dataframe
        dataframe_seir34!(result, sol)
        # and remove the last row (which is the same day as the first row for the next year)
        filter!(:t => x -> x != τe, result)

        # u for the new year 
        u = last(sol)

        # age individuals 
        for q ∈ 1:4     # q = { S, E, I, R }
            u[q, 1] = u[q, 1] * 5 / 6                       # aging of youngest group 
            u[q, 2] = u[q, 1] * 1 / 6 + u[q, 2] * 3 / 4     # next group, etc. 
            u[q, 3] = u[q, 2] * 1 / 4 + u[q, 3] * 9 / 10
            u[q, 4] = u[q, 3] * 1 / 10 + u[q, 4]  
        end 

        # time for new year 
        τs = τe
    end 

    return result
end 

function dataframe_seir34!(df, sol)
    df2 = dataframe_seir34(sol)
    append!(df, df2)
end 

function dataframe_seir34(sol)
    result = DataFrame(t = sol.t)
    for i ∈ 1:4 
        _dataframe_seir34!(result, sol, i)
    end 
    return result 
end 

function _dataframe_seir34!(result, sol, row)
    categories = [ "S", "E", "I", "R" ]
    for (i, cat) ∈ enumerate(categories)
        rs = zeros( size(sol, 3) )
        for j ∈ axes(sol, 3) 
            rs[j] = sol[j][i, row] 
        end 
        insertcols!(result, Symbol("$cat$row") => rs)
    end 
end 

function plot_seir34(result; kwargs...)
    fig = Figure()
    plot_seir34!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

function plot_seir34!(fig::Figure, result; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_seir34!(gl, result; kwargs...)
end 

function plot_seir34!(gl::GridLayout, result; 
        label = "p3.3.jl: Susceptible--exposed--infectious--recovered model with four age groups", 
        legend = :right
    )

    axs = [ Axis(gl[i, 1]) for i ∈ 1:4 ]
    for i ∈ 1:4 
        plot_seir34!(
            axs[i], 
            select(result, :t, Symbol("S$i"), Symbol("E$i"), Symbol("I$i"), Symbol("R$i"))
        )
        i < 4 && hidexdecorations!(axs[i]; grid = false, ticks = false)
        Label(gl[i, 2], "age $i"; rotation = -π/2)
    end 
    linkxaxes!(axs...)
    Label(gl[:, 0], "Proportion of the population"; rotation = π/2)
    Label(gl[5, :], "Time, years")

    Label(gl[0, :], label)

    if legend == :right
        leg = Legend(gl[:, 3], axs[1])
    elseif legend == :below 
        leg = Legend(gl[6, 1], axs[1]; orientation = :horizontal)
    elseif legend == :none 
        # no legend 
    else 
        @info "Unrecognised legend position given. Recognised options are `:right`, `:below` and `:none`"
    end
end 

function plot_seir34!(ax::Axis, result)
    lbls = [ "Susceptible", "Exposed", "Infectious", "Recovered"] 
    for (i, lbl) ∈ enumerate(lbls)
        lines!(ax, result.t ./ 365, result[!, names(result)[i+1]]; label = lbl)
    end 
end 


## Help text

"""
    run_sir33(u0, p, duration[; saveat])
    run_seir34(u0, p, duration[; saveat])

Run the model `seir34!`, a susceptible--exposed--infectious--recovered model with 
four age  groups, and include annual aging between runs of the model.

## Arguments 

`u0` is a `4x4` matrix of starting conditions for the model. There is a column for 
each age group with values from top to bottom in the order susceptible, exposed, 
infectious, recovered. 

`p` is the parameters for the model. This is expected in a `Parameters34` type. 
`p.beta` is an `4x4` matrix of infectiousness parameters between each age group. 
`p.sigma`is the rate at which exposed individuals become infectious. `p.gamma` is 
the recovery rate. `p.mu` is the mortality rate. `p.nu` is the birth rate.

`duration` is the length of time that the model should run.

## Optional keyword arguments 
* `saveat`: how frequently the ODE solver should save results. Default is `1`
"""
function run_seir34() end

end # module MID_34
