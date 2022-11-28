
module MID_73

# Coupled lattice model with commuter-like coupling (page 256)
  
using CairoMakie, DifferentialEquations, Random
using StatsBase: sample

export sir73_u0, sir73!, run_sir73, plot_sir73, video_sir73
 
# Neighbours can be in one of four directions 
@enum Directions North East South West

function sir73!(du, u, p, t) 
    # compartments 
    X = u[:, :, 1]
    Y = u[:, :, 2]
    Z = u[:, :, 3]
    # parameters
    beta, gamma, mu, rho = p

    N = .+(X, Y, Z)

    for i ∈ axes(u, 1), j ∈ axes(u, 2) 
        nu = mu * N[i, j]   # so that births equal deaths

        # ODEs

        # dX
        du[i, j, 1] = nu - beta * X[i, j] * ( 
                ( 1 - rho * countneighbours(Y, i, j) ) * Y[i, j] + rho * sumneighbours(Y, i, j) 
            ) / ( 
                ( 1 - rho * countneighbours(N, i, j) ) * N[i, j] + rho * sumneighbours(N, i, j) 
            ) - 
            mu * X[i, j]

        # dY 
        du[i, j, 2] = beta * X[i, j] * ( 
                ( 1 - rho * countneighbours(Y, i, j) ) * Y[i, j] + rho * sumneighbours(Y, i, j) 
            ) / ( 
                ( 1 - rho * countneighbours(N, i, j) ) * N[i, j] + rho * sumneighbours(N, i, j) 
            ) - 
            (gamma + mu) * Y[i, j]

        # dZ 
        du[i, j, 3] = gamma * Y[i, j] - mu * Z[i, j]
    end
end 

function run_sir73(; n, x0, yvector = nothing, ni = yvector, y0, n0, beta, gamma, 
        mu, rho, duration, seed = nothing, kwargs...
    )
    u0 = sir73_u0(n, x0, ni, y0, n0; seed)
    p = [beta, gamma, mu, rho]
    return run_sir73(u0, p, duration; kwargs...)
end 

function run_sir73(u0, p, duration; saveat = 1)
    @assert minimum(u0) >= 0 "Input u0 = $u0: no compartments can contain negative proportions"
    @assert minimum(p) >= 0 "Input p = $p: cannot run with negative parameters" 
    @assert duration > 0 "Input duration = $duration: cannot run with negative or zero duration"

    tspan = ( 0., Float64(duration) )

    prob = ODEProblem(sir73!, u0, tspan, p)
    sol = solve(prob; saveat)
    return sol
end 

sir73_u0(; n, x0, ni, y0, n0, seed = nothing) = sir73_u0(n, x0, ni, y0, n0; seed)
sir73_u0(n, x0, ni, y0, n0; seed = nothing) = _sir73_u0(n, x0, ni, y0, n0, seed)

function _sir73_u0(n, x0, ni, y0, n0, seed)
    @error "Input `n` must be an integer and `ni` must be an integer or vector"
end 

function _sir73_u0(n::Int, x0, ni, y0, n0, seed::Real)
    Random.seed!(seed)
    return _sir73_u0(n, x0, ni, y0, n0, nothing)
end 

function _sir73_u0(n::Int, x0, ni::Int, y0, n0, seed::Nothing)
    yvector = sample(collect(1:1:n^2), ni; replace = false)
    return _sir73_u0(n, x0, yvector, y0, n0, seed)
end 

function _sir73_u0(n::Int, x0, yvector::Vector{<:Int}, y0, n0, seed::Nothing)
    @assert n0 >= x0 + y0 "x0 + y0 > n0 will lead to some negative compartments"

    X0 = x0 * ones(n, n)
    Y0 = zeros(n, n) 
    for y ∈ yvector 
        Y0[y] = y0 
    end 
    Z0 = zeros(n, n) 
    for i ∈ 1:n, j ∈ 1:n 
        Z0[i, j] = n0 - (X0[i, j] + Y0[i, j])
    end 
    u0 = zeros(n, n, 3) 
    u0[:, :, 1] = X0; u0[:, :, 2] = Y0; u0[:, :, 3] = Z0
    return u0
end 

sumneighbours(x, i, j) = sum( findneighbours(x, i, j) ) # sums the values of all neighbouring cells

# return the values of all neighbouring cells:
findneighbours(x, i, j) = [ findneighbour(x, i, j, Directions(d)) for d ∈ 0:3 ]

function findneighbour(x, i, j, direction) # returns the value of a neighbouring cell
    if direction == North 
        if i == 1 
            return .0 
        else 
            return x[i-1, j] 
        end
    elseif direction == East 
        if j == size(x, 2)
            return .0 
        else 
            return x[i, j+1] 
        end 
    elseif direction == South 
        if i == size(x, 1)
            return .0 
        else 
            return x[i+1, j] 
        end 
    else # direction == West 
        if j == 1 
            return .0 
        else 
            return x[i, j-1]
        end 
    end 
end 

function countneighbours(x, i, j) 
    # number of neighbours is usually 4 but not for cells on an edge or in a corner
    n = 0 
    if i > 1 n += 1 end 
    if j < size(x, 2) n += 1 end 
    if i < size(x, 1) n += 1 end 
    if j > 1 n += 1 end  
    return n 
end 

function plot_sir73(sol, i::Int; forcepositive = false, kwargs...)
    data = values_sir73(sol, i; forcepositive)
    return plot_sir73(data, sol; forcepositive, kwargs...) 
end 

function plot_sir73(data, sol; 
        fixmax = false, colormap = :viridis, 
        label = "p7.3.jl: Coupled lattice model with commuter-like coupling"
    )
    fig = Figure()
    if fixmax 
        maxval = maximum([ maximum(sol[i][:, :, 2]) for i ∈ axes(sol, 4) ])
        colorrange = (0, maxval)
        ax, hm = heatmap(fig[1, 1][1, 1], data; colorrange, colormap)
    else 
        ax, hm = heatmap(fig[1, 1][1, 1], data; colormap)
    end   
    Colorbar(fig[1, 1][1, 2], hm)
    hidexdecorations!(ax)
    hideydecorations!(ax)
    #Label(fig[0, :], label)
    return fig 
end 

function values_sir73(sol, t; forcepositive = false) 
    d = sol[t][:, :, 2]
    if forcepositive 
        # an uncomfortable fix for situations where some compartments have very slightly negative values
        for i ∈ axes(d, 1), j ∈ axes(d, 2)
            d[i, j] = max(.0, d[i, j])
        end 
    end 
    return d
end 

function video_sir73(sol; folder = "outputvideos", filename = "video73.mp4", forcepositive = false, kwargs...)
    data = values_sir73(sol, 1; forcepositive)
    values = Observable(data)
    fig = plot_sir73(values, sol; kwargs...)

    record(fig, "$folder/$filename") do io
        for i ∈ axes(sol, 4)
            values[] = values_sir73(sol, i; forcepositive)  # animate scene by changing values
            recordframe!(io)                                # record a new frame
        end
    end
end 
 
end # module MID_73
