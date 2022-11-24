
module MID_72
  
using CairoMakie, DataFrames, DifferentialEquations
import Base: minimum

export Parameters72, sir72!, run_sir72, plot_sir72, plot_sir72!

struct Parameters72
    beta    :: Vector{<:Float64}
    gamma   :: Vector{<:Float64}
    mu      :: Vector{<:Float64}
    nu      :: Matrix{<:Float64}
    l       :: Matrix{<:Float64}
    r       :: Matrix{<:Float64}
end

# Define the Base.minimum for Parameters71
minimum(p::Parameters72) = min( 
    minimum(p.beta), 
    minimum(p.gamma), 
    minimum(p.mu), 
    minimum(p.nu), 
    minimum(p.l),
    minimum(p.r)
)

function sir72!(du, u, p, t) 
    # compartments 
    X = u[:, :, 1]
    Y = u[:, :, 2]
    Z = u[:, :, 3]

    N = .+(X, Y, Z)
    
    # ODEs
    for i ∈ axes(u, 1), j ∈ axes(u, 2) 
        # some actions will appear multiple times in the equations 
        births = p.nu[i, j]
        infections = p.beta[i] * X[i, j] * sum(Y[i, :]) / sum(N[i, :])
        recoveries = p.gamma[i] * Y[i, j]
        deaths_x = p.mu[i] * X[i, j]
        deaths_y = p.mu[i] * Y[i, j]
        deaths_z = p.mu[i] * Z[i, j]

        if i == j 
            # dS
            du[i, j, 1] = births - infections - deaths_x - 
                sum(p.l[:, i]) * X[i, j] +                          # leaving 
                sum([ p.r[a, i] * X[a, i] for a ∈ axes(u, 2) ]) +   # returning
                p.l[i, j] * X[j, j] - p.r[i, j] * X[i, j]           # undo terms 
                        # for [i, i] in case p.r[i, i] != 0 or p.l[i, i] != 0

            # dI 
            du[i, j, 2] = infections - recoveries - deaths_y -
                sum(p.l[:, i]) * Y[i, j] +                          # leaving 
                sum([ p.r[a, i] * Y[a, i] for a ∈ axes(u, 2) ]) +   # returning
                p.l[i, j] * Y[j, j] - p.r[i, j] * Y[i, j]           # undo terms for [i, i]

            # dZ 
            du[i, j, 3] = recoveries - deaths_z -
                sum(p.l[:, i]) * Z[i, j] +                          # leaving 
                sum([ p.r[a, i] * Z[a, i] for a ∈ axes(u, 2) ]) +   # returning
                p.l[i, j] * Z[j, j] - p.r[i, j] * Z[i, j]           # undo terms for [i, i]
        
        else 
            # dS
            du[i, j, 1] = births - infections - deaths_x + 
                p.l[i, j] * X[j, j] -                               # arriving
                p.r[i, j] * X[i, j]                                 # returning

            # dI 
            du[i, j, 2] = infections - recoveries - deaths_y +
                p.l[i, j] * Y[j, j] -                               # arriving
                p.r[i, j] * Y[i, j]                                 # returning

            # dZ 
            du[i, j, 3] = recoveries - deaths_z +   
                p.l[i, j] * Z[j, j] -                               # arriving
                p.r[i, j] * Z[i, j]                                 # returning
        
        end 
    end
end 

function run_sir72(u0, p, duration; saveat = 1)
    @assert minimum(u0) >= 0 "Input u0 = $u0: no compartments can contain negative proportions"
    @assert minimum(p) >= 0 "Input p = $p: cannot run with negative parameters" 
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"

    tspan = ( 0., Float64(duration) )

    prob = ODEProblem(sir72!, u0, tspan, p)
    sol = solve(prob; saveat)

    return sol
end 

function dataframe_sir72(sol, i, j)
    result = DataFrame(t = sol.t)
    for (k, c) ∈ enumerate([ "Susceptible", "Infectious", "Recovered" ])
        rs = zeros( size(sol, 4) )
        for t ∈ axes(sol, 4) 
            rs[t] = sol[t][i, j, k] 
        end 
        insertcols!(result, Symbol("$c") => rs)
    end 
    return result 
end 

function plot_sir72(sol; kwargs...)
    fig = Figure()
    plot_sir72!(fig, sol; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

function plot_sir72!(fig::Figure, sol; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_sir72!(gl, sol; kwargs...)
end 

function plot_sir72!(gl::GridLayout, sol; 
        label = "p7.2.jl: SIR metapopulation model", legend = :below, kwargs...
    )

    @assert axes(sol, 1) == axes(sol, 2) 
    solaxs = axes(sol, 1)
    solsz = size(sol, 1)

    axs = [ Axis(gl[i, j]) for i ∈ solaxs, j ∈ solaxs ]
    for i ∈ solaxs, j ∈ solaxs
        plot_sir72!( axs[i, j], sol, i, j; hidex = i < solsz, hidey = j != 1, kwargs... )
    end 
    Label(gl[solsz+1, :], "Time")
    Label(gl[:, 0], "Counts"; rotation = π/2)
    Label(gl[0, :], label)
    linkaxes!(axs...)

    if legend == :right
        leg = Legend(gl[:, solsz+1], axs[1])
    elseif legend == :below 
        leg = Legend(gl[solsz+2, :], axs[1]; orientation = :horizontal)
    elseif legend == :none 
        # no legend 
    else 
        @info "Unrecognised legend position given. Recognised options are `:right`, `:below` and `:none`"
    end 
end 

function plot_sir72!(ax::Axis, sol, i, j; 
        hidex = false, hidey = false, lbls = [ "Susceptible", "Infectious", "Recovered" ]
    )
    
    data = dataframe_sir72(sol, i, j)
    for lbl ∈ lbls
        lines!(ax, data.t, data[:, lbl]; label = lbl)
    end 
    hidex && hidexdecorations!(ax; grid = false, ticks = false)
    hidey && hideydecorations!(ax; grid = false, ticks = false)
end 

end # module MID_72
