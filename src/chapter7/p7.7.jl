
module MID_77
  
using CairoMakie, DataFrames, Distributions, GraphMakie, Graphs, Random, StatsBase

export sis77!, u0_sis77, run_sis77, rundf_sis77, dataframe_sis77, plot_sis77, plot_sis77!, 
    graphplot_sis77!, video_sis77

mutable struct Environment 
    g           :: SimpleGraph{Int64} 
    Y           :: Vector{Int64}
    historyY    :: Vector{Vector{Int}}
end

struct SpatialPosition 
    x           :: Float64 
    y           :: Float64 
end 

function u0_sis77(N, averageconnections, Y0, networktype; kwargs...)
    @assert Y0 <= N "Cannot have more than all individuals infectious"
    @assert averageconnections <= N - 1 "Need averageconnections ≤ N - 1: 
        Cannot connect with more than every other node"

    if networktype == :random || networktype == :Random 
        g = makenetwork_random(N, averageconnections; kwargs...)
    elseif networktype == :lattice || networktype == :Lattice 
        g = makenetwork_lattice(N, averageconnections; kwargs...)
    elseif networktype == :smallworld || networktype == :Smallworld || networktype == :SmallWorld
        g = makenetwork_smallworld(N, averageconnections; kwargs...)
    elseif networktype == :spatial || networktype == :Spatial
        g = makenetwork_spatial(N, averageconnections; kwargs...)
    else 
        @error "`networktype` not recognised"
    end 
    Y = sample(vertices(g), Y0; replace = false)

    return Environment(g, Y, [Y])
end 

function makenetwork_random(N, averageconnections; seed = nothing)
    distances = ones(N, N) # all probabilities are equal
    return __makenetwork(N, averageconnections, distances, seed) 
end 

makenetwork_lattice(N, averageconnections; kwargs...) = watts_strogatz(N, averageconnections, 0)
    
makenetwork_smallworld(N, averageconnections; seed = nothing, beta = .05) = 
    watts_strogatz(N, averageconnections, beta; seed)

function makenetwork_spatial(N, averageconnections; seed = nothing, spacesize = 1) 
    if spacesize > 100 
        @warn "Nodes further apart than ≈ 150 will have a probability of connecting of 0" 
    end

    # Make a matrix of positions and identify which nodes link to each other
    positions = [ randomposition(spacesize) for _ ∈ 1:N ]
    distances = [ dists(positions[i], positions[j]) for i ∈ 1:N, j ∈ 1:N ]
    return __makenetwork(N, averageconnections, distances, seed) 
end 

function __makenetwork(N, averageconnections, distances, seed::Real) 
    Random.seed!(seed)
    return __makenetwork(N, averageconnections, distances, nothing) 
end 

function __makenetwork(N, averageconnections, distances, seed::Nothing) 
    probabilities = [ calcprobs(distances, i, j) for i ∈ 1:N, j ∈ 1:N ]
    wts = ProbabilityWeights(vec(probabilities))
    identities = vec([ (i, j) for i ∈ 1:N, j ∈ 1:N ])
    contacts = round(Int, N * averageconnections / 2) # total number of connections
    connecteds = sample(identities, wts, contacts; replace = false) # vector of connections

    g = SimpleGraph(N)
    for c ∈ connecteds
        add_edge!(g, c...)
    end 
    
    return g
end 

randomposition(spacesize) = SpatialPosition(rand(Uniform(0, spacesize), 2)...)

dists(A, B) = sqrt( (A.x - B.x)^2 + (A.y - B.y)^2 )

function calcprobs(distances, i, j)
    if i >= j # connections are bidirectional so non-zero values in only half of matrix
        return .0 
    else 
        return exp(-5 * distances[i, j])
    end 
end 

function sis77!(u, p, t; tstep = nothing)
    # Parameters 
    gamma, tau = p 
    
    # Rates of events 
    rates = [
        ifelse(
            i ∈ u.Y, # i.e. is infected 
            gamma,
            forceofinfection(i, u, tau)
        ) for i ∈ vertices(u.g)
    ]

    return _sis77!(u, t, rates, tstep)
end 

function _sis77!(u, t, rates, tstep::Nothing)
    # No tstep is provided so model goes to the time when the next event happens 
    sumrates = sum(rates)

    # Random numbers to decide what happens 
    r1 = rand(); r2 = rand() 

    # How soon is something going to happen 
    timestep = -log(r1) / sumrates
    newt = t + timestep

    # Who is it going to happen to  
    cumsumratesratio = cumsum(rates) / sumrates
    i = searchsortedfirst(cumsumratesratio, r2)

    # Calculate new vector of infected individuals 
    Y = Int[] 
    for j ∈ vertices(u.g) 
        if j ∈ u.Y # i.e. is already infected 
            if j != i push!(Y, j) end # not changing, remains infectious 
        else 
            if j == i push!(Y, j) end # changing, becomes infectious
        end 
    end 
    u.Y = Y 
    push!(u.historyY, Y)

    return u, newt
end 

function _sis77!(u, t, rates, tstep::Real)
    # Calculate how many events occur during period tstep 
    probs = [ 1 - exp(-rate * tstep) for rate in rates ]
    rands = rand(length(rates))

    # Calculate new vector of infected individuals 
    Y = Int[] 
    for i ∈ vertices(u.g) 
        if i ∈ u.Y  # i.e. is already infected 
            if rands[i] >= probs[i] push!(Y, i) end # not changing, remains infectious 
        else        # not already infected
            if rands[i] < probs[i] push!(Y, i) end # changing, becomes infectious
        end 
    end 

    u.Y = Y 
    push!(u.historyY, Y)

    return u, t + tstep
end 

function forceofinfection(i, u, tau)
    lambda = .0 
    for j ∈ neighbors(u.g, i) 
        if j ∈ u.Y 
            lambda += tau 
        end 
    end 
    return lambda  
end 

run_sis77(u0, p, duration; seed = nothing, tstep = nothing, kwargs...) = 
    _run_sis77(u0, p, duration, seed; tstep)

function run_sis77(; N, averageconnections, Y0, networktype, gamma, tau, duration, 
        seed = nothing, tstep = nothing, kwargs...
    )
    u0 = u0_sis77(N, averageconnections, Y0, networktype; seed, kwargs...)
    p = [gamma, tau]
    return run_sis77(u0, p, duration; seed, tstep)
end

function rundf_sis77(args...; N, kwargs...)
    u, times = run_sis77(args...; N, kwargs...)
    df = dataframe_sis77(u, times, N)
    return u, times, df
end 

function _run_sis77(u0, p, duration, seed::Int; tstep)
    Random.seed!(seed)
    return _run_sis77(u0, p, duration, nothing; tstep)
end 

function _run_sis77(u0, p, duration, seed::Nothing; tstep)
    @assert minimum(p) >= 0 "Input p = $p. Model cannot run with negative parameters."
    @assert duration > 0 "Input duration = $duration. Model needs duration > 0."

    u = u0
    t = .0
    times = [t]
    
    while t < duration 
        u, t = sis77!(u, p, t; tstep)
        push!(times, t)
    end 

    return u, times
end 

function dataframe_sis77(u, times, N)
    df = DataFrame(t = Float64[], Susceptible = Int[], Infectious = Int[])
    for (i, t) ∈ enumerate(times)
        infec = length(u.historyY[i])
        susce = N - infec
        append!(
            df, 
            DataFrame(t = t, Susceptible = susce, Infectious = infec)
        )
    end 
    return df 
end 

function plot_sis77(data, t = nothing)
    fig = Figure()
    plot_sis77!(fig, data, t)
    return fig 
end 

function plot_sis77!(fig::Figure, data, t = nothing)
    gl = GridLayout(fig[1, 1])
    plot_sis77!(gl, data, t)
end 

function plot_sis77!(gl::GridLayout, data, t = nothing)
    ax = Axis(gl[1, 1])
    plot_sis77!(ax, data, t)
    leg = Legend(gl[2, 1], ax; orientation = :horizontal)
end 

function plot_sis77!(ax::Axis, data, t = nothing)
    copydata = deepcopy(data)
    if last(data.t) == Inf pop!(copydata) end 
    lines!(ax, copydata.t, copydata.Susceptible; label = "Susceptibles")
    lines!(ax, copydata.t, copydata.Infectious; label = "Infectious")
    plott_sis77!(ax, t)
end 

plott_sis77!(ax, t::Nothing) = nothing 
plott_sis77!(ax, t) = vlines!(ax, t)

function graphplot_sis77!(ax, g, yhistory, t)
    colours = [ ifelse(i ∈ yhistory[t], :red, :black) for i ∈ vertices(g) ]
    graphplot!(ax, g; node_color = colours )
end 

function videoplot_sis77(g, df, t, colours)
    fig = Figure()
    ga = GridLayout(fig[1, 1])
    ax1 = Axis(ga[1, 1])
    graphplot!(ax1, g; node_color = colours )
    hidexdecorations!(ax1); hideydecorations!(ax1)
    gb = GridLayout(fig[1, 2])
    ax2 = Axis(gb[1, 1])
    plot_sis77!(ax2, df, t)
    leg = Legend(gb[2, 1], ax2)
    return fig
end 

function video_sis77(u, times, df; 
        step = 1/24, folder = "outputvideos", filename = "video77.mp4", kwargs...)

    copytimes = deepcopy(times) 
    if last(times) == Inf pop!(copytimes) end 

    g = u.g
    yhistory = u.historyY
    tstart = 0.
    colours = [ ifelse(i ∈ yhistory[1], :red, :black) for i ∈ vertices(g) ]
    ot = Observable(tstart)
    cols = Observable(colours)

    fig = videoplot_sis77(g, df, ot, cols)

    frametimes = collect(0:step:last(copytimes))

    record(fig, "$folder/$filename") do io
        for t ∈ frametimes 
            # animate scene by changing values:
            ot[] = t
            k = max(1, searchsortedfirst(copytimes, t) - 1)
            cols[] = [ ifelse(i ∈ yhistory[k], :red, :black) for i ∈ vertices(g) ]
            recordframe!(io)    # record a new frame
        end
    end
end 

end # module MID_77
