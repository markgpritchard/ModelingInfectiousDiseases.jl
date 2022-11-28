
module MID_33

# SIR model with 2 age classes (page 79)
  
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

function run_sir33(; S0, I0, R0, beta, gamma, lambda, mu, nu, duration, kwargs...)
    @assert length(S0) == length(I0) == length(R0) "Must have equal length vectors for S0, I0 and R0"
    u0 = zeros(3, length(S0))
    u0[1, :] = S0 
    u0[2, :] = I0 
    u0[3, :] = R0 
    p = Parameters33(beta, gamma, lambda, mu, nu)
    return run_sir33(u0, p, duration; kwargs...)
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

minimum(p::Parameters33) = min( minimum(p.beta), p.gamma, p.lambda, minimum(p.mu), p.nu )

function dataframe_sir33(sol)
    result = DataFrame(t = sol.t)
    dataframe_sir33!(result, sol, :c, 1)
    dataframe_sir33!(result, sol, :a, 2)
    return result 
end 

function dataframe_sir33!(result, sol, typesymbol, row)
    categories = [ "S", "I", "R" ]
    for (i, cat) ∈ enumerate(categories)
        rs = zeros( size(sol, 3) )
        for j ∈ axes(sol, 3) 
            rs[j] = sol[j][i, row] 
        end 
        insertcols!(result, Symbol("$cat$typesymbol") => rs)
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

end # module MID_33
