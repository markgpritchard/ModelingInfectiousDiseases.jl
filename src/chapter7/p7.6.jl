
module MID_76
  
using CairoMakie, DataFrames, Distributions, Random
import Base: minimum

export Parameters76, seirc76, u0_seirc76, run_seirc76, dataframe_seirc76, plot_seirc76, 
    plot_seirc76!, video_seirc76

# Each individual can have one of two states  
@enum State Susceptible Exposed Infectious Reported Culled RingCulled 

mutable struct Farm
    x               :: Float64
    y               :: Float64
    sheep           :: Float64
    cows            :: Float64
    state           :: State
    timeinstate     :: Int
    sus             :: Real 
    trans           :: Real
end

struct Environment
    t               :: Int 
    farms           :: Vector{Farm}
end

struct Parameters76 
    s               :: Vector{<:Float64}
    t               :: Vector{<:Float64}
    Ring            :: Float64
end

minimum(p::Parameters76) = min(minimum(p.s), minimum(p.t), p.Ring)
minimum(u::Environment) = min(u.t, minimum(u.farms))
function minimum(v::Vector{Farm}) 
    m = Inf
    for f ∈ v
        m = min(m, minimum(f))
    end 
    return m 
end 

minimum(f::Farm) = min(f.x, f.y, f.sheep, f.cows, f.timeinstate, f.sus, f.trans)

function setfarm(size, state, p)
    @assert minimum(p) >= 0 

    sheep = 350 * exp(rand())
    cows = 70 * exp(rand())

    return Farm(
        rand(Uniform(0, size)),
        rand(Uniform(0, size)),
        sheep,
        cows,
        state,
        0,
        p.s[1] * sheep + p.s[2] * cows,
        p.t[1] * sheep + p.t[2] * cows
    )
end 

u0_seirc76(n, Y0, size, p; seed = nothing) = _u0_seirc76(n, Y0, size, p, seed)

function _u0_seirc76(n, Y0, size, p, seed::Real)
    Random.seed!(seed)
    return _u0_seirc76(n, Y0, size, p, nothing)
end 

function _u0_seirc76(n, Y0, size, p, seed::Nothing)
    @assert Y0 <= n "Cannot have more than all individuals infectious"

    farms = Farm[]
    for i ∈ 1:n 
        if i <= Y0 
            push!(farms, setfarm(size, Infectious, p))
        else 
            push!(farms, setfarm(size, Susceptible, p))
        end 
    end 

    return Environment(0, farms)
end 

function seirc76(u, p; tstep = 1)
    newu = deepcopy(u) # so that _seirc76! does not mutate the original 
    return _seirc76!(newu, p, tstep)
end 

function _seirc76!(u, p, tstep)
    farms = u.farms

    # First make a vector of infectious farms 
    infecteds = Farm[] 
    for X ∈ farms 
        if X.state == Infectious || X.state == Reported
            push!(infecteds, X) 
        end 
    end 

    ## Calculated rates of possible infections 
    infectionrates = [
        ifelse(
            farm.state == Susceptible, 
            forceofinfection(farm, infecteds),
            .0
        )
        for farm ∈ farms
    ]

    # Calculate sum of all rates 
    sumrates = sum(infectionrates)

    # Calculate how many events occur during period tstep 
    probs = [ 1 - exp(-rate * tstep) for rate in infectionrates ]
    rands = rand(length(farms))

    # Make the changes 
    for (i, farm) ∈ enumerate(farms) 
        if farm.state == Susceptible 
            if rands[i] < probs[i] 
                farm.state = Exposed
                farm.timeinstate = 0 
            end 
        elseif farm.state == Exposed 
            if farm.timeinstate == 5
                farm.state = Infectious 
                farm.timeinstate = 0 
            end
        elseif farm.state == Infectious 
            if farm.timeinstate == 5
                farm.state = Reported 
                farm.timeinstate = 0 
            end
        elseif farm.state == Reported 
            if farm.timeinstate == 2
                farm.state = Culled 
                farm.timeinstate = 0 
            end
        end  
    end 

    # Perform ring culling 
    # First make a vector of farms that could be at the centre of a ring cull 
    cullcentre = Farm[] 
    for X ∈ farms 
        if X.state == Culled && X.timeinstate == 2
            push!(cullcentre, X) 
        end 
    end 

    if length(cullcentre) > 0 
        # Find farms within Ring distance of cullcentre and cull them 
        for A ∈ cullcentre 
            for B ∈ farms
                dist = sqrt( (A.x - B.x)^2 + (A.y - B.y)^2 ) 
                if dist <= p.Ring && B.state != Culled 
                    B.state = RingCulled 
                    B.timeinstate = 0 
                end 
            end 
        end 
    end 

    # Add one on to all time in state values 

    for farm ∈ farms 
        farm.timeinstate += 1 
    end 

    # increase time of model 
    newt = u.t + tstep

    return Environment(newt, farms)
end 

function transmissionkernel(A, B)
    d2 = (A.x - B.x)^2 + (A.y - B.y)^2 
    if d2 < 0.0138 
        k = .3093 
    elseif d2 > 60^2 
        k = .0 
    else 
        k = exp(-9.2123e-5 * d2^6 + 9.5628e-4 * d2^5 + 3.3966e-3 * d2^4 - 3.3687e-2 * d2^3 - 
            1.30519e-1 * d2^2 - 0.609262 * d2 - 3.231772) 
    end 
    return k 
end 

function forceofinfection(A, infecteds)
    ksum = .0 
    for J ∈ infecteds 
        ksum += transmissionkernel(A, J) * J.trans
    end 
    lambda = A.sus * ksum
    return lambda 
end 

run_seirc76(u0, p, duration; seed = nothing, tstep = 1) = 
    _run_seirc76(u0, p, duration, seed; tstep)

function _run_seirc76(u0, p, duration, seed::Int; tstep)
    Random.seed!(seed)
    return _run_seirc76(u0, p, duration, nothing; tstep)
end 

function _run_seirc76(u0, p, duration, seed::Nothing; tstep)
    @assert minimum(u0) >= 0 "Input u0 = $u0. Model cannot run with negative compartments."
    @assert minimum(p) >= 0 "Input p = $p. Model cannot run with negative parameters."
    @assert duration > 0 "Input duration = $duration. Model needs duration > 0."

    u = u0
    t = u.t
    results = [u]
    
    while t < duration 
        nextu = seirc76(u, p; tstep)
        push!(results, nextu)
        t = nextu.t
        u = nextu
    end 

    return results
end 

function dataframe_seirc76(result)
    df = DataFrame(t = Int[], Susceptible = Int[], Exposed = Int[], Infectious = Int[], 
        Reported = Int[], Culled = Int[], RingCulled = Int[])
    for i ∈ axes(result, 1)
        susce = 0; expos = 0; infec = 0; repor = 0; culle = 0; ringc = 0 
        for farm ∈ result[i].farms
            if farm.state == Susceptible 
                susce += 1 
            elseif farm.state == Exposed 
                expos += 1 
            elseif farm.state == Infectious 
                infec += 1 
            elseif farm.state == Reported 
                repor += 1 
            elseif farm.state == Culled 
                culle += 1 
            else # farm.state == RingCulled 
                ringc += 1 
            end 
        end 
        d = DataFrame(t = result[i].t, Susceptible = susce, Exposed = expos, Infectious = infec, 
            Reported = repor, Culled = culle, RingCulled = ringc)
        append!(df, d)
    end 
    return df 
end 

function plot_seirc76(data, t = nothing)
    fig = Figure()
    plot_seirc76!(fig, data, t)
    return fig 
end 

function plot_seirc76!(fig::Figure, data, t = nothing)
    gl = GridLayout(fig[1, 1])
    plot_seirc76!(gl, data, t)
end 

function plot_seirc76!(gl::GridLayout, data, t = nothing)
    ax = Axis(gl[1, 1])
    plot_seirc76!(ax, data, t)
    leg = Legend(gl[2, 1], ax; orientation = :horizontal)
end 

function plot_seirc76!(ax::Axis, data, t = nothing)
    lines!(ax, data.t, data.Susceptible; label = "Susceptibles")
    lines!(ax, data.t, data.Exposed; label = "Exposed")
    lines!(ax, data.t, data.Infectious; label = "Infectious")
    lines!(ax, data.t, data.Reported; label = "Reported")
    lines!(ax, data.t, data.Culled; label = "Culled")
    lines!(ax, data.t, data.RingCulled; label = "RingCulled")
    plott_seirc76!(ax, t)
end 

plott_seirc76!(ax, t::Nothing) = nothing 
plott_seirc76!(ax, t) = vlines!(ax, t)

function videoplot_seirc76(s_xs, s_ys, e_xs, e_ys, i_xs, i_ys, r_xs, r_ys, c_xs, c_ys, rc_xs, rc_ys, data, t)
    fig = Figure()
    ga = GridLayout(fig[1, 1])
    ax = Axis(ga[1, 1])
    scatter!(ax, s_xs, s_ys; label = "Susceptibles", markersize = 3)
    scatter!(ax, e_xs, e_ys; label = "Exposed", markersize = 3)
    scatter!(ax, i_xs, i_ys; label = "Infectious", markersize = 3)
    scatter!(ax, r_xs, r_ys; label = "Reported", markersize = 3)
    scatter!(ax, c_xs, c_ys; label = "Culled", markersize = 3)
    scatter!(ax, rc_xs, rc_ys; label = "RingCulled", markersize = 3)
    gb = GridLayout(fig[1, 2])
    ax2 = Axis(gb[1, 1])
    plot_seirc76!(ax2, data, t)
    leg = Legend(fig[2, :], ax; orientation = :horizontal)
    hidexdecorations!(ax)
    hideydecorations!(ax)
    return fig 
end 

values_seirc76(env::Environment) = values_seirc76(env.farms)

function values_seirc76(farms::Vector{Farm})
    s_xs, s_ys = values_seirc76(farms, Susceptible)
    e_xs, e_ys = values_seirc76(farms, Exposed)
    i_xs, i_ys = values_seirc76(farms, Infectious)
    r_xs, r_ys = values_seirc76(farms, Reported)
    c_xs, c_ys = values_seirc76(farms, Culled)
    rc_xs, rc_ys = values_seirc76(farms, RingCulled)
    return [s_xs, s_ys, e_xs, e_ys, i_xs, i_ys, r_xs, r_ys, c_xs, c_ys, rc_xs, rc_ys]
end 

function values_seirc76(farms, state)
    xs = Float64[]; ys = Float64[] 
    for farm ∈ farms
        if farm.state == state 
            push!(xs, farm.x); push!(ys, farm.y)
        else 
            push!(xs, NaN); push!(ys, NaN)
        end 
    end 
    return xs, ys
end 

function video_seirc76(results, data; 
        step = 1/24, folder = "outputvideos", filename = "video76.mp4", kwargs...)

    tvector = [ result.t for result ∈ results]
    s_xs, s_ys, e_xs, e_ys, i_xs, i_ys, r_xs, r_ys, c_xs, c_ys, rc_xs, rc_ys = values_seirc76(results[1])
    oxs = Observable(s_xs)
    oys = Observable(s_ys)
    oxe = Observable(e_xs)
    oye = Observable(e_ys)
    oxi = Observable(i_xs)
    oyi = Observable(i_ys)
    oxr = Observable(r_xs)
    oyr = Observable(r_ys)
    oxc = Observable(c_xs)
    oyc = Observable(c_ys)
    oxrc = Observable(rc_xs)
    oyrc = Observable(rc_ys)
    ot = Observable(tvector[1])

    fig = videoplot_seirc76(oxs, oys, oxe, oye, oxi, oyi, oxr, oyr, oxc, oyc, oxrc, oyrc, data, ot)

    frametimes = tvector

    record(fig, "$folder/$filename") do io
        for t ∈ frametimes 
            s_xs, s_ys, e_xs, e_ys, i_xs, i_ys, r_xs, r_ys, c_xs, c_ys, rc_xs, rc_ys = values_seirc76(results[t+1])
            # animate scene by changing values:
            oxs[] = s_xs
            oys[] = s_ys
            oxe[] = e_xs
            oye[] = e_ys
            oxi[] = i_xs
            oyi[] = i_ys
            oxr[] = r_xs
            oyr[] = r_ys
            oxc[] = c_xs
            oyc[] = c_ys
            oxrc[] = rc_xs
            oyrc[] = rc_ys
            ot[] = t
            recordframe!(io)    # record a new frame
        end
    end
end 

end # module MID_76
