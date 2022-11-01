
module MID_71
  
using CairoMakie, DataFrames, DifferentialEquations
import Base: minimum

export Parameters71, sir71!, run_sir71, dataframe_sir71, plot_sir71, plot_sir71!

struct Parameters71 
    beta    :: Vector{<:Float64}
    gamma   :: Vector{<:Float64}
    mu      :: Vector{<:Float64}
    nu      :: Vector{<:Float64}
    m       :: Matrix{<:Float64}
end

# Define the Base.minimum for Parameters71
minimum(p::Parameters71) = min( 
    minimum(p.beta), 
    minimum(p.gamma), 
    minimum(p.mu), 
    minimum(p.nu), 
    minimum(p.m) 
)

function sir71!(du, u, p, t) 
    # compartments 
    X = u[1, :]
    Y = u[2, :]
    Z = u[3, :]
    
    # ODEs
    for i ∈ axes(u, 2)
        du[1, i] = p.nu[i] - p.beta[i] * X[i] * Y[i] - p.mu[i] * X[i] + 
            sum( [ p.m[i, j] * X[j] for j ∈ axes(u, 2) ] ) - 
            sum( [ p.m[j, i] * X[i] for j ∈ axes(u, 2) ] )              # dX 
        du[2, i] = p.beta[i] * X[i] * Y[i] - p.gamma[i] * Y[i] - p.mu[i] * Y[i] + 
            sum( [ p.m[i, j] * Y[j] for j ∈ axes(u, 2) ] ) - 
            sum( [ p.m[j, i] * Y[i] for j ∈ axes(u, 2) ] )              # dY
        du[3, i] = p.gamma[i] * Y[i] - p.mu[i] * Z[i] + 
            sum( [ p.m[i, j] * Z[j] for j ∈ axes(u, 2) ] ) - 
            sum( [ p.m[j, i] * Z[i] for j ∈ axes(u, 2) ] )              # dZ
    end 
end 

function run_sir71(u0, p, duration; saveat = 1)
    @assert minimum(u0) >= 0 "Input u0 = $u0: no compartments can contain negative proportions"
    @assert minimum(p) >= 0 "Input p = $p: cannot run with negative parameters" 
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"

    tspan = ( 0., Float64(duration) )

    prob = ODEProblem(sir71!, u0, tspan, p)
    sol = solve(prob; saveat)

    return sol
end 

function dataframe_sir71(sol)
    result = DataFrame(t = sol.t)
    for i ∈ axes(sol, 2), (j, c) ∈ enumerate(["S", "I", "R"])
        rs = zeros( size(sol, 3) )
        for k ∈ axes(sol, 3) 
            rs[k] = sol[k][j, i] 
        end 
        insertcols!(result, Symbol("$c$i") => rs)
    end 
    return result 
end 

function plot_sir71(result; kwargs...)
    fig = Figure()
    plot_sir71!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_sir71!(any, sol; kwargs...) = plot_sir71!(any, dataframe_sir71(sol); kwargs...)

function plot_sir71!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_sir71!(gl, result; kwargs...)
end 

function plot_sir71!(gl::GridLayout, result::DataFrame; 
        label = "p7.1.jl: SIR metapopulation model", legend = :belowsol72_2 = run_sir72(u0_2, p_2, duration_2; saveat = .125)
        plot_sir72(sol72_2) )

    n = Int((size(result, 2) - 1) / 3)

    axs = [ Axis(gl[i, 1]) for i ∈ 1:n ]
    for i ∈ 1:n 
        plot_sir71!( axs[i], select(result, :t, Symbol("S$i"), Symbol("I$i"), Symbol("R$i")) )
    end 
    axs[n].xlabel = "Time, years" 
    Label(gl[:, 0], "Counts"; rotation = π/2)
    Label(gl[0, :], label)

    if legend == :right
        leg = Legend(gl[:, 2], axs[1])
    elseif legend == :below 
        leg = Legend(gl[n+1, 1], axs[1]; orientation = :horizontal)
    elseif legend == :none 
        # no legend 
    else 
        @info "Unrecognised legend position given. Recognised options are `:right`, `:below` and `:none`"
    end 
end 

function plot_sir71!(ax::Axis, result::DataFrame)
    lbls = [ "Susceptible", "Infectious", "Recovered" ]
    for (i, lbl) ∈ enumerate(lbls)
        lines!(ax, result.t ./ 365, result[!, i+1]; label = lbl)
    end 
end 

end # module MID_6_1
