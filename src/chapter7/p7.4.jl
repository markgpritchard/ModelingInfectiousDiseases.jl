
module MID_74
  
using CairoMakie, Random
using StatsBase: sample

export sir74!, u0_sir74, run_sir74, plot_sir74, video_sir74

# Each cell can have one of three states  
@enum States S I R 

u0_sir74(n, I0 = 0, R0 = 0; seed = nothing) = _u0_sir74(n, I0, R0, seed)

function _u0_sir74(n, I0, R0, seed::Real)
    Random.seed!(seed)
    return _u0_sir74(n, I0, R0, nothing)
end 

function _u0_sir74(n, I0, R0, seed::Nothing)
    @assert I0 + R0 <= n^2 "Cannot have more than all individuals susceptible or recovered"

    # make a square grid of individuals 
    grid = fill(S, (n, n))

    # which ones are not susceptible 
    notsusceptible = sample(collect(1:1:n^2), I0 + R0; replace = false)

    # which ones are infectious and which are recovered 
    if I0 > 0
        infectious = notsusceptible[1: I0]
        for i ∈ infectious grid[i] = I end 
    end  
    if R0 > 0 
        recovered  = notsusceptible[I0+1: end]
        for i ∈ recovered grid[i] = R end 
    end
    
    return grid
end 

function sir74!(u, p, t, neighbours, rates) 
    tau, gamma, nu, epsilon = p 

    # Each individual can only do one thing. Therefore, only one rate applies to 
    # each cell.
    # How many infectious neighbours does each cell have? 
    for i ∈ axes(u, 1), j ∈ axes(u, 2) 
        neighbours[i, j] = infectiousneighbours(u, i, j) 
    end 
    # Rate of relevant event for each cell 
    for i ∈ axes(u, 1), j ∈ axes(u, 2) 
        if u[i, j] == S 
            rates[i, j] = neighbours[i, j] * tau + epsilon
        elseif u[i, j] == I 
            rates[i, j] = gamma 
        else # u[i, j] == R 
            rates[i, j] = nu
        end  
    end 

    # Random numbers to decide what happens 
    r1 = rand(); r2 = rand() 

    # How soon is something going to happen 
    timestep = -log(r1) / sum(rates)
    newt = last(t) + timestep
    push!(t, newt)

    # What is going to happen 
    cumsumratesratio = cumsummat(rates) / sum(rates)
    i = searchsortedfirstmat(cumsumratesratio, r2)
    if u[i] == S 
        u[i] = I
    elseif u[i] == I 
        u[i] = R
    else # u[i] == R 
        u[i] = S
    end  
end 

function infectiousneighbours(x, i, j) 
    if x[i, j] != S return 0 end # only need to know about neighbours of susceptibles
    n = 0 
    # Look North
    if i > 1 && x[i-1, j] == I n += 1 end
    # Look East 
    if j < size(x, 2) && x[i, j+1] == I n += 1 end
    # Look South 
    if i < size(x, 1) && x[i+1, j] == I n += 1 end
    # Look West 
    if j > 1 && x[i, j-1] == I n += 1 end

    return n 
end 

function cumsummat(mat) # function cumsum not doing what I want, so making own version 
    tot = size(mat, 1) * size(mat, 2)
    result = zeros(size(mat, 1), size(mat, 2))
    result[1] = mat[1]
    for i ∈ 2:tot 
        result[i] = mat[i] + result[i-1]
    end 
    return result 
end 

function searchsortedfirstmat(mat, v)
    tot = size(mat, 1) * size(mat, 2)
    i = 1
    while 1 < tot + 1 
        if mat[i] >= v return i end 
        i += 1 
    end 
end 

function run_sir74(u0, p, duration; seed = nothing)
    ucopy = deepcopy(u0) # run_sir74! mutates the u0 input
    return run_sir74!(ucopy, p, duration, seed)
end 

function run_sir74(; n, I0 = 0, R0 = 0, tau, gamma, nu, epsilon, duration, seed = nothing, kwargs...)
    u0 = u0_sir74(n, I0, R0; seed)
    p = [tau, gamma, nu, epsilon]
    # u0_sir74 has reset the global random number generator seed so does not also 
    # need to be passed to run_sir74!. Pass directly to run_sir74! rather than 
    # run_sir74(u0, p, duration; seed = nothing) as u0 was created by this function 
    # so can be safely mutated.
    return run_sir74!(u0, p, duration, nothing; kwargs...)
end 

function run_sir74!(u0, p, duration, seed::Int)
    Random.seed!(seed)
    return run_sir74!(u0, p, duration, nothing)
end 

function run_sir74!(u0, p, duration, seed::Nothing)
    @assert minimum(p) >= 0 "Model cannot run with negative parameters. Running with p = $p."
    @assert duration > 0 "Model needs duration > 0. Model supplied duration = $duration."

    u = u0
    t = .0
    neighbours = zeros(Int, size(u0, 1), size(u0, 2))
    rates = zeros(size(u0, 1), size(u0, 2))

    tvector = [t]
    ulatest = copy(u)
    uvector = [ulatest]
    
    while t <= duration 
        sir74!(u, p, tvector, neighbours, rates)
        ulatest = copy(u)
        push!(uvector, ulatest)
        t = last(tvector)
    end 

    return tvector, uvector 
end 

function plot_sir74(uvector, i::Int; kwargs...)
    data = values_sir74(uvector, i)
    return plot_sir74(data; kwargs...) 
end 

function plot_sir74(data; colormap = :viridis)
    fig = Figure()
    ax, hm = heatmap(fig[1, 1][1, 1], data; colormap, colorrange = (0, 2))
    Colorbar(fig[1, 1][1, 2], hm, ticks = ([0, 1, 2], ["S", "I", "R"]))
    hidexdecorations!(ax)
    hideydecorations!(ax)
    return fig 
end 

function values_sir74(uvector, i) 
    u = uvector[i]
    d = zeros(Int, size(u, 1), size(u, 2))
    for i ∈ axes(u, 1), j ∈ axes(u, 2) 
        if u[i, j] == I
            d[i, j] = 1 
        elseif u[i, j] == R
            d[i, j] = 2
        end 
    end 
    return d 
end 

function video_sir74(uvector, tvector; 
        step = 1/24, folder = "outputvideos", filename = "video74.mp4", kwargs...)

    data = values_sir74(uvector, 1)
    values = Observable(data)
    fig = plot_sir74(values; kwargs...)

    frametimes = collect(0:step:last(tvector))

    record(fig, "$folder/$filename") do io
        for i ∈ frametimes 
            u = max(1, searchsortedfirst(tvector, i) - 1)
            values[] = values_sir74(uvector, u)     # animate scene by changing values
            recordframe!(io)                        # record a new frame
        end
    end
end 

end # module MID_74
