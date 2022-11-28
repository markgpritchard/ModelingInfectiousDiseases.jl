
module MID_35

# SEIR model with n stages (page 94)

# I am deviating from the example in having a seperate sigma parameter for the exposed 
# compartments
  
using CairoMakie, DataFrames, DifferentialEquations
import Base: minimum

export Parameters35, seir35!, run_seir35, dataframe_seir35, separateddataframe_seir35, 
    plot_seir35, plot_seir35!, seir35_u0

struct Parameters35
    beta    :: Float64
    sigma   :: Float64
    gamma   :: Float64
    mu      :: Float64  
    nu      :: Float64
    m       :: Int 
    n       :: Int
end

function seir35!(du, u, p, t)  
    # Compartments will be referred to only with their location in the vector  

    # Number of compartments 
    m = p.m; n = p.n

    # Parameters 
    β = p.beta
    σ = p.sigma * m         # to keep total time in exposed compartments 1 / p.sigma 
    γ = p.gamma * (n - m)   # to keep total time in infectious compartments 1 / p.gamma 
    μ = p.mu    
    nu = p.nu  

    # ODEs 
    du[1] = nu - β * u[1] * sum(u[m+2:n+1]) - μ * u[1]                  # dS 
    if m > 0        # if there is at least one E compartment
        du[2] = β * u[1] * sum(u[m+2:n+1]) - (σ + μ) * u[2]             # dE1 
        if m > 1    # if there is more than one E compartment 
            du[3:m+1] = [ σ * u[i-1] - (σ + μ) * u[i] for i ∈ 3:m+1 ]   # dE2:dEm
        end 
                    # if there is at least one E compartment
        du[m+2] = σ * u[m+1] - (γ + μ) * u[m+2]                         # dI1 
    else            # if there are no E compartments 
        du[2] = β * u[1] * sum(u[m+2:n+1]) - (γ + μ) * u[2]             # dI1 
    end 
    if n - m > 1
        du[m+3:n+1] = [ γ * u[i-1] - (γ + μ) * u[i] for i ∈ m+3:n+1 ]   # dI2:dIn
    end 
    du[n+2] = γ * u[n+1] - μ * u[n+2]                                   # dR
end 

function run_seir35(; m, n, S0, E0, I0, R0 = 1 - (S0 + E0 + I0), beta, sigma, gamma, mu, nu = mu, duration, kwargs...)
    u0 = seir35_u0(S0, E0, I0, R0, m, n)
    p = Parameters35(beta, sigma, gamma, mu, nu, m, n)
    return run_seir35(u0, p, duration; kwargs...)
end

function run_seir35(u0, p, duration; saveat = 1)
    @assert minimum(u0) >= 0 "Input u0 = $u0: no compartments can contain negative proportions"
    @assert sum(u0) ≈ 1 "Input u0 = $u0: compartment values are proportions and must sum to 1"
    @assert minimum(p) >= 0 "Input p = $p: cannot run with negative parameters" 
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"

    tspan = ( 0., Float64(duration) )

    prob = ODEProblem(seir35!, u0, tspan, p)
    sol = solve(prob; saveat)

    return sol
end 

minimum(p::Parameters35) = min( p.beta, p.sigma, p.gamma, p.mu, p.nu, p.m, p.n )

# u0 needs to be the right size for the function, so use a function to create it 
function seir35_u0(S, E, I, R, m::Int, n::Int)
    @assert m < n "`m` is exposed and `n` is exposed + infectious, no need n < m"
    @assert +(S, E, I, R) ≈ 1 "Compartments are proportions so must sum to 1"
    @assert min(S, E, I, R, m, n) >= 0 "Cannot take any negative values"

    u = zeros(2 + n)
    
    u[1] = S 

    if m == 0 
        @assert E == 0 "When m = 0 there are no exposed compartments so need E == 0" 
    else 
        e_prop = E / m
        for i ∈ 2:m+1 
            u[i] = e_prop 
        end 
    end 

    i_prop = I / (n - m)  
    for i ∈ m+2:n+1 
        u[i] = i_prop
    end 

    u[n+2] = R
    
    return u
end 

seir35_u0(S, E, I, R, p::Parameters35) = seir35_u0(S, E, I, R, p.m, p.n)
seir35_u0(S, E, I, m::Int, n::Int) = seir35_u0(S, E, I, 1 - (S + E + I), m, n) 
seir35_u0(S, E, I, p::Parameters35) = seir35_u0(S, E, I, 1 - (S + E + I), p)

dataframe_seir35(sol, p) = dataframe_seir35(sol, p.m, p.n)

function dataframe_seir35(sol, m, n)
    rs = zeros( size(sol, 2) )
    re = zeros( size(sol, 2) )
    ri = zeros( size(sol, 2) )
    rr = zeros( size(sol, 2) )

    for i ∈ axes(sol, 1)
        if i == 1           # susceptible
            for j ∈ axes(sol, 2) 
                rs[j] = sol[j][i] 
            end 
        elseif i <= m + 1   # exposed
            for j ∈ axes(sol, 2) 
                re[j] = sol[j][i] 
            end 
        elseif i <= n + 1   # infectious
            for j ∈ axes(sol, 2) 
                ri[j] = sol[j][i] 
            end 
        else                # recovered
            for j ∈ axes(sol, 2) 
                rr[j] = sol[j][i] 
            end 
        end 
    end 
    result = DataFrame(t = sol.t, Susceptible = rs, Exposed = re, Infectious = ri, Recovered = rr)
    return result 
end 

separateddataframe_seir35(sol, p) = separateddataframe_seir35(sol, p.m, p.n)

function separateddataframe_seir35(sol, m, n)
    result = DataFrame(t = sol.t)
    for i ∈ axes(sol, 1)
        if i == 1
            name = "S" 
        elseif i <= m + 1 
            name = "E$(i-1)" 
        elseif i <= n + 1 
            name = "I$(i-m-1)" 
        else 
            name = "R" 
        end 
        rs = zeros( size(sol, 2) )
        for j ∈ axes(sol, 2) 
            rs[j] = sol[j][i] 
        end 
        insertcols!(result, Symbol(name) => rs, makeunique=true)
    end 
    return result 
end 

function plot_seir35(result; kwargs...)
    fig = Figure()
    plot_seir35!(fig, result; kwargs...)
    resize_to_layout!(fig)
    return fig 
end 

function plot_seir35!(fig::Figure, result; kwargs...)
    gl = GridLayout(fig[1, 1])
    plot_seir35!(gl, result; kwargs...)
end 

function plot_seir35!(gl::GridLayout, result; 
        label = "p3.5.jl: Susceptible--exposed--infectious--recovered model with multiple stages", 
        legend = :right
    )
    ax = Axis(gl[1, 1])
    plot_seir35!(ax, result)
    ax.xlabel = "Time, years" 
    ax.ylabel = "Proportion"
    Label(gl[0, :], label)

    if legend == :right
        leg = Legend(gl[1, 2], ax)
    elseif legend == :below 
        leg = Legend(gl[2, 1], ax; orientation = :horizontal)
    elseif legend == :none 
        # no legend 
    else 
        @info "Unrecognised legend position given. Recognised options are `:right`, `:below` and `:none`"
    end
end 

function plot_seir35!(ax::Axis, result)
    lbls = names(result)
    popfirst!(lbls)
    for (i, lbl) ∈ enumerate(lbls)
        lines!(ax, result.t ./ 365, result[!, i+1]; label = lbl)
    end 
end 

end # module MID_35
