
module MID_75
  
using CairoMakie, Distributions, Random

export sis75, u0_sis75, run_sis75, plot_sis75, video_sis75

# Each individual can have one of two states  
@enum State S I 

mutable struct Individual
    x               :: Float64
    y               :: Float64
    state           :: State
    rate            :: Float64
end

struct Environment
    t               :: Float64 
    individuals     :: Vector{Individual}
end

function setindividual(size, state)
    return Individual(
        rand(Uniform(0, size)),
        rand(Uniform(0, size)),
        state,
        .0
    )
end 

u0_sis75(n, Y0, size = 10; seed = nothing) = _u0_sis75(n, Y0, size, seed)

function _u0_sis75(n, Y0, size, seed::Real)
    Random.seed!(seed)
    return _u0_sis75(n, Y0, size, nothing)
end 

function _u0_sis75(n, Y0, size, seed::Nothing)
    @assert Y0 <= n "Cannot have more than all individuals infectious"

    individuals = Individual[]
    for i ∈ 1:n 
        if i <= Y0 
            push!(individuals, setindividual(size, I))
        else 
            push!(individuals, setindividual(size, S))
        end 
    end 

    return Environment(.0, individuals)
end 

function sis75(u, p; tstep = nothing)
    newu = deepcopy(u) # so that _sis75! does not mutate the original 
    return _sis75!(newu, p, tstep)
end 

function _sis75!(u, p, tstep)
    # Parameters 
    alpha, beta, gamma = p 

    inds = u.individuals

    ## Calculated rates of each possible event 

    # First make a vector of infected individuals 
    infecteds = Individual[] 
    for X ∈ inds 
        if X.state == I push!(infecteds, X) end 
    end 

    # Calculate a rate for each individual 
    for X ∈ inds
        if X.state == S 
            X.rate = forceofinfection(X, infecteds, alpha, beta)
        else # X.state == I 
            X.rate = gamma
        end 
    end 

    # Calculate sum of all rates 
    rates = [ ind.rate for ind ∈ inds ]
    
    if sum(rates) == 0 # there will be no more events 
        return Environment(Inf, inds)
    else 
        return __sis75(inds, rates, u.t, tstep)
    end 
end 

function __sis75(inds, rates, t, tstep::Nothing)
    # No tstep is provided so model goes to the time when the next event happens 
    sumrates = sum(rates)

    # Random numbers to decide what happens 
    r1 = rand(); r2 = rand() 

    # How soon is something going to happen 
    timestep = -log(r2) / sumrates
    newt = t + timestep

    # Who is it going to happen to  
    cumsumratesratio = cumsum(rates) / sumrates
    i = searchsortedfirst(cumsumratesratio, r2)

    # Make the change 
    if inds[i].state == S 
        inds[i].state = I
    else # inds[i].state == I
        inds[i].state = S 
    end 

    return Environment(newt, inds)
end 

function __sis75(inds, rates, t, tstep::Real)
    # Calculate how many events occur during period tstep 
    probs = [ 1 - exp(-rate * tstep) for rate in rates ]
    rands = rand(length(inds))

    # Make the changes 
    for (i, ind) ∈ enumerate(inds)
        if rands[i] < probs[i] # a change occurs 
            if ind.state == S 
                ind.state = I 
            else # ind.state == I 
                ind.state = S 
            end 
        end 
    end 

    newt = t + tstep

    return Environment(newt, inds)
end 

function transmissionkernel(A, B, alpha)
    distance = sqrt( (A.x - B.x)^2 + (A.y - B.y)^2 )
    k = distance^(-alpha) 
    return k 
end 

function forceofinfection(A, infecteds, alpha, beta)
    ksum = .0 
    for J ∈ infecteds 
        ksum += transmissionkernel(A, J, alpha)
    end 
    lambda = beta * ksum
    return lambda 
end 

run_sis75(u0, p, duration; seed = nothing, tstep = nothing) = 
    _run_sis75(u0, p, duration, seed; tstep)

function _run_sis75(u0, p, duration, seed::Int; tstep)
    Random.seed!(seed)
    return _run_sis75(u0, p, duration, nothing; tstep)
end 

function _run_sis75(u0, p, duration, seed::Nothing; tstep)
    @assert minimum(p) >= 0 "Input p = $p. Model cannot run with negative parameters."
    @assert duration > 0 "Input duration = $duration. Model needs duration > 0."

    u = u0
    t = u.t
    results = [u]
    
    while t <= duration 
        nextu = sis75(u, p; tstep)
        push!(results, nextu)
        t = nextu.t
        u = nextu
    end 

    return results
end 

function plot_sis75(s_xs, s_ys, i_xs, i_ys)
    fig = Figure()
    ax = Axis(fig[1, 1])
    scatter!(ax, s_xs, s_ys; label = "Susceptibles")
    scatter!(ax, i_xs, i_ys; label = "Infecteds")
    leg = Legend(fig[2, 1], ax; orientation = :horizontal)
    hidexdecorations!(ax)
    hideydecorations!(ax)
    return fig 
end 

values_sis75(env::Environment) = values_sis75(env.individuals)

function values_sis75(individuals::Vector{Individual})
    s_xs, s_ys = values_sis75(individuals, S)
    i_xs, i_ys = values_sis75(individuals, I)
    return [s_xs, s_ys, i_xs, i_ys]
end 

function values_sis75(individuals, state)
    xs = Float64[]; ys = Float64[] 
    for ind ∈ individuals
        if ind.state == state 
            push!(xs, ind.x); push!(ys, ind.y)
        else 
            push!(xs, NaN); push!(ys, NaN)
        end 
    end 
    return xs, ys
end 

function video_sis75(results; 
        step = 1/24, folder = "outputvideos", filename = "video75.mp4", kwargs...)

    tvector = [ result.t for result ∈ results]
    if last(tvector) == Inf
        pop!(tvector) 
        @info "Model reached point where no more events could occur. Video terminates at that time"
    end 
    s_xs, s_ys, i_xs, i_ys = values_sis75(results[1])
    sxvalues = Observable(s_xs)
    syvalues = Observable(s_ys)
    ixvalues = Observable(i_xs)
    iyvalues = Observable(i_ys)

    fig = plot_sis75(sxvalues, syvalues, ixvalues, iyvalues)

    frametimes = collect(0:step:last(tvector))

    record(fig, "$folder/$filename") do io
        for i ∈ frametimes 
            x = max(1, searchsortedfirst(tvector, i) - 1)
            s_xs, s_ys, i_xs, i_ys = values_sis75(results[x])
            # animate scene by changing values:
            sxvalues[] = s_xs
            syvalues[] = s_ys
            ixvalues[] = i_xs
            iyvalues[] = i_ys
            recordframe!(io)    # record a new frame
        end
    end
end 

end # module MID_75
