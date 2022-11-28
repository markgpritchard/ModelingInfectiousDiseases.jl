
module MID_27

# SIR model with carrier state (page 44)

using CairoMakie, DataFrames, DifferentialEquations

export sir27!, run_sir27, dataframe_sir27, plot_sir27, plot_sir27!

function sir27!(du, u, p, t)
    # compartments 
    S, I, C = u 
    # parameters 
    beta, gamma_i, gamma_c, epsilon, mu, nu, q = p 
    # Note that gamma (γ) and Gamma (Γ) are described as different parameters for 
    # this model. Within this code we use gamma_i for gamma and gamma_c for Gamma

    # the ODEs
    du[1] = nu - beta * S * (I + epsilon * C) - mu * S              # dS
    du[2] = beta * S * (I + epsilon * C) - gamma_i * I - mu * I     # dI
    du[3] = gamma_i * q * I - gamma_c * C - mu * C                  # dC
end 

function run_sir27(; S0, I0, C0, beta, gamma_i, gamma_c, epsilon, mu, nu = mu, q, duration, kwargs...)
    u0 = [S0, I0, C0]
    p = [beta, gamma_i, gamma_c, epsilon, mu, nu, q]
    return run_sir27(u0, p, duration; kwargs...)
end 

function run_sir27(u0, p, duration; saveat = 1)
    @assert minimum(u0) >= 0 "Input u0 = $u0: Cannot run model with negative starting values in any compartment"
    @assert sum(u0) <= 1 "Input u0 = $u0: Compartment values are proportions so must sum to 1 or less"
    @assert minimum(p) >= 0 "Input p = $p: Cannot run with negative values for any parameters"
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"

    tspan = ( 0., Float64(duration) )

    prob = ODEProblem(sir27!, u0, tspan, p)
    sol = solve(prob; saveat)
    return sol
end 

function dataframe_sir27(sol)
    result = DataFrame(t = Float64[], S = Float64[], I = Float64[], C = Float64[])
    for i ∈ eachindex(sol.u)
        push!( result, Dict(
            :t => sol.t[i], 
            :S => sol.u[i][1], 
            :I => sol.u[i][2], 
            :C => sol.u[i][3]
        ) )
    end 
    return result 
end 

function plot_sir27(result; kwargs...)
    fig = Figure()
    plot_sir27!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_sir27!(any, sol; kwargs...) = plot_sir27!(any, dataframe_sir27(sol); kwargs...)

function plot_sir27!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_sir27!(gl, result; kwargs...)
end 

function plot_sir27!(gl::GridLayout, result::DataFrame; 
        label = "p2.7.jl: Susceptible--infectious--resistant model with a carrier state", 
        legend = :right, kwargs...
    )
    ax = Axis(gl[1, 1])
    plot_sir27!(ax, result; kwargs...)

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

function plot_sir27!(ax::Axis, result::DataFrame)
    xs = result.t ./ 365  # to plot time in years
    
    lines!(ax, xs, result.S, label = "Susceptible")
    lines!(ax, xs, result.I, label = "Infectious")
    lines!(ax, xs, result.C, label = "Carriers")
    ax.xlabel = "Time, years"
    ax.ylabel = "Fraction of population"
end 

end # module MID_27
