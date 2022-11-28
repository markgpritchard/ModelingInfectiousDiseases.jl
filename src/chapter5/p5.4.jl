
module MID_54

# Rabbit Hemorrhagic Disease model (page 186)
  
using CairoMakie, DataFrames, DifferentialEquations

export sir54!, run_sir54, dataframe_sir54, plot_sir54, plot_sir54!

function sir54!(du, u, p, t) 
    # compartments 
    X, Y, N = u
    # Parameters 
    α0, α1, β0, β1, γ, μ, ω, m, K = p

    # Seasonal value of α and β
    α = α0 * (1 + α1 * sin(ω * t))
    β = β0 * (1 + β1 * sin(ω * t))
    
    # the ODEs
    du[1] = α * N - β * X * Y - (μ + N / K) * X     # dX
    du[2] = β * X * Y - (γ + m + μ + N / K) * Y     # dY
    du[3] = (α - μ - N / K) * N - m * Y             # dN 
    # -γ term removed from dN as these are not deaths and do not reduce the population total
end 

function run_sir54(; X0, Y0, N0, alpha0, alpha1, beta0, beta1, gamma, mu, omega = 2pi / 365, m, K, duration, kwargs...)
    u0 = [X0, Y0, N0]
    p = [alpha0, alpha1, beta0, beta1, gamma, mu, omega, m, K]
    return run_sir54(u0, p, duration; kwargs...)
end 

function run_sir54(u0, p, duration; saveat = 1, kwargs...)
    @assert minimum(u0) >= 0 "Input u0 = $u0: cannot run with negative compartment values"
    @assert u0[1] + u0[2] <= u0[3] "Input u0 = $u0: N0 must be ≥ X0 + Y0"
    @assert minimum(p) >= 0 "Input p = $p: cannot run with negative parameters"
    @assert p[2] <= 1 "Input α1 = $(p[2]): if α1 > 1 then at some times the beta parameter will be negative"
    @assert p[4] <= 1 "Input β1 = $(p[2]): if β1 > 1 then at some times the beta parameter will be negative"
    @assert duration > 0 "Input duration = $duration: cannot run with zero or negative duration"

    tspan = ( 0., Float64(duration) )

    prob = ODEProblem(sir54!, u0, tspan, p)
    sol = solve(prob; saveat, kwargs...)

    return sol
end 

function dataframe_sir54(sol)
    result = DataFrame(t = Float64[], X = Float64[], Y = Float64[], N = Float64[])
    for i ∈ eachindex(sol.u)
        push!( result, Dict(
            :t => sol.t[i], 
            :X => sol.u[i][1], 
            :Y => sol.u[i][2], 
            :N => sol.u[i][3]
        ) )
    end 
    return result 
end 

function plot_sir54(result; kwargs...)
    fig = Figure()
    plot_sir54!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

plot_sir54!(any, sol; kwargs...) = plot_sir54!(any, dataframe_sir54(sol); kwargs...)

function plot_sir54!(fig::Figure, result::DataFrame; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_sir54!(gl, result; kwargs...)
end 

function plot_sir54!(gl::GridLayout, result::DataFrame; 
        label = "p5.3.jl: Rabbit Hemorrhagic Disease model", 
        legend = :right, kwargs...
    )
    ax = Axis(gl[1, 1])
    plot_sir54!(ax, result; kwargs...)

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

function plot_sir54!(ax::Axis, result::DataFrame; plotN = true)
    xs = result.t ./ 365  # to plot time in weeks
    
    lines!(ax, xs, result.X, label = "Susceptible")
    lines!(ax, xs, result.Y, label = "Infectious")
    plotN && lines!(ax, xs, result.N, label = "Population")
    ax.xlabel = "Time, years"
    ax.ylabel = "Number"
end 

end # module MID_54
