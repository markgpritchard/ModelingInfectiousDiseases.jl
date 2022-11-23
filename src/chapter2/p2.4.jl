
module MID_24
  
using CairoMakie, DataFrames, DifferentialEquations

export sir24!, run_sir24, dataframe_sir24, plot_sir24, plot_sir24!

function sir24!(du, u, p, t)
    # compartments 
    X, Y, Z = u # following convention in book, {S, I, R} refer to proportions and {X, Y, Z} refer to numbers
    # parameters 
    beta, gamma, mu, nu, rho = p 

    # total population size 
    N = X + Y + Z

    # the ODEs
    du[1] = nu - beta * X * Y / N - mu * X                      # dX
    du[2] = beta * X * Y / N - (gamma + mu) * Y / (1 - rho)     # dY
    du[3] = gamma * Y - mu * Z                                  # dZ
end 

function run_sir24(u0, p, duration; reltol = 1e-12, saveat = 1)
    # set low tolerance for solver, otherwise oscillations don't converge appropriately
    @assert minimum(u0) >= 0 "Input u0 = $u0: Cannot run model with negative starting values in any compartment"
    @assert minimum(p) >= 0 "Input p = $p: Cannot run with negative values for any parameters"
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"

    tspan = ( 0., Float64(duration) )

    prob = ODEProblem(sir24!, u0, tspan, p)
    sol = solve(prob; reltol, saveat)

    return sol
end 

function run_sir24(; N0 = 0, X0, Y0, Z0 = N0 - (X0 + Y0), beta, gamma, mu, nu, rho, duration, kwargs...)
    u0 = [X0, Y0, Z0]
    p = [beta, gamma, mu, nu, rho]
    return run_sir24(u0, p, duration; kwargs...)
end 

function dataframe_sir24(sol)
    result = DataFrame(t = Float64[], X = Float64[], Y = Float64[], Z = Float64[])
    for i ∈ eachindex(sol.u)
        push!( result, Dict(
            :t => sol.t[i], 
            :X => sol.u[i][1], 
            :Y => sol.u[i][2], 
            :Z => sol.u[i][3]
        ) )
    end 
    insertcols!(result, :N => +(result.X, result.Y, result.Z) )
    return result 
end 

function plot_sir24(result; kwargs...)
    fig = Figure()
    plot_sir24!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_sir24!(any, sol; kwargs...) = plot_sir24!(any, dataframe_sir24(sol); kwargs...)

function plot_sir24!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_sir24!(gl, result; kwargs...)
end 

function plot_sir24!(gl::GridLayout, result::DataFrame; 
        label = "p2.4.jl: SIR model with infection-induced mortality and frequency-dependent transmission", 
        legend = :right, kwargs...
    )
    ax = Axis(gl[1, 1])
    plot_sir24!(ax, result; kwargs...)

    Label(gl[0, :], label)

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

function plot_sir24!(ax::Axis, result::DataFrame; plotZ = true, plotN = true)
    xs = result.t ./ 365  # to plot time in years
    
    lines!(ax, xs, result.X, label = "Susceptible")
    lines!(ax, xs, result.Y, label = "Infectious")
    plotZ && lines!(ax, xs, result.Z, label = "Recovered")
    plotN && lines!(ax, xs, result.N, label = "Total population")
    ax.xlabel = "Time, years"
    ax.ylabel = "Numbers"
end 

end # module MID_24
