
module MID_23
  
using CairoMakie, DataFrames, DifferentialEquations

export sir23!, run_sir23, dataframe_sir23, plot_sir23, plot_sir23!

function sir23!(du, u, p, t)
    # compartments 
    X, Y, Z = u # following convention in book, {S, I, R} refer to proportions and {X, Y, Z} refer to numbers
    # parameters 
    beta, gamma, mu, nu, rho = p 

    # the ODEs
    du[1] = nu - beta * X * Y - mu * X                      # dX
    du[2] = beta * X * Y - (gamma + mu) * Y / (1 - rho)     # dY
    du[3] = gamma * Y - mu * Z                              # dZ
end 

function run_sir23(u0, p, duration; reltol = 1e-12, saveat = 1)
    # set low tolerance for solver, otherwise oscillations don't converge appropriately
    @assert minimum(u0) >= 0 "Input u0 = $u0: Cannot run model with negative starting values in any compartment"
    @assert minimum(p) >= 0 "Input p = $p: Cannot run with negative values for any parameters"
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"

    tspan = ( 0., Float64(duration) )

    prob = ODEProblem(sir23!, u0, tspan, p)
    sol = solve(prob; reltol, saveat)

    return sol
end 

function run_sir23(; N0 = 0, X0, Y0, Z0 = N0 - (X0 + Y0), beta, gamma, mu, nu = mu, rho, duration, kwargs...)
    u0 = [X0, Y0, Z0]
    p = [beta, gamma, mu, nu, rho]
    return run_sir23(u0, p, duration; kwargs...)
end 

function dataframe_sir23(sol)
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

function plot_sir23(result; kwargs...)
    fig = Figure()
    plot_sir23!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_sir23!(any, sol; kwargs...) = plot_sir23!(any, dataframe_sir23(sol); kwargs...)

function plot_sir23!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_sir23!(gl, result; kwargs...)
end 

function plot_sir23!(gl::GridLayout, result::DataFrame; 
        label = "p2.3.jl: SIR model with infection-induced mortality and density-dependent transmission", 
        legend = :right, kwargs...
    )
    ax = Axis(gl[1, 1])
    plot_sir23!(ax, result; kwargs...)

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

function plot_sir23!(ax::Axis, result::DataFrame; plotZ = true, plotN = true)
    xs = result.t ./ 365  # to plot time in years
    
    lines!(ax, xs, result.X, label = "Susceptible")
    lines!(ax, xs, result.Y, label = "Infectious")
    plotZ && lines!(ax, xs, result.Z, label = "Recovered")
    plotN && lines!(ax, xs, result.N, label = "Total population")
    ax.xlabel = "Time, years"
    ax.ylabel = "Numbers"
end 

#="""
    run_sir23([; beta, gamma, mu, nu, rho, X0, Y0, N0, duration, saveat])

Run the model `sir23`

# Keyword arguments 

All keyword arguments are optional with default values supplied for each. 

* `beta`: The beta parameter in the model (infectiousness of infectives). Default is 520 / 365 (520 per year).
* `gamma`: The gamma parameter in the model (recovery rate). Default is 1 / 7.
* `mu`: The model's mortality rate not due to the pathogen. Default is 1 / (70 * 365) 
    (i.e. mean life duration 70 years).
* `nu`: The model's birth rate. Default is 1 / (70 * 365) (i.e. match the mortality in the absence of the pathogen).
* `rho`: The mortality probability for infecteds. Default is 0.5.
* `X0`: Number susceptible at time = 0. Default is 0.2. Note, we follow the convention 
    used in the book that {S, I, R} represent proportions and {X, Y, Z} represent 
    numbers susceptible, infectious and resistant.
* `Y0`: Number infectious at time = 0. Default is 1e-6.
* `N0`: The initial total population size. Default is 1.
* `duration`: How long the model will run (time units are interpretted as days). 
    Default is 100 years. (The default differs between the programmes in different 
    languages -- 100 years appears to give a reasonable illustration of the results.)
* `saveat`: How frequently the model should save values. Default is 1 (day).
""" =#

end # module MID_23
