
using ModelingInfectiousDiseases.MID_2_6
using Suppressor, Test
using DifferentialEquations

# macro isinferred(ex) suggested by Andras Niedermayer on 
# https://discourse.julialang.org/t/more-informative-inferred-especially-for-unit-testing/5481
macro isinferred(ex)
    quote try
        @inferred $ex
            true
        catch err
            println(err)
            false
        end
    end
end

macro noerrors(ex)
    quote try
            @suppress $ex
            true
        catch err
            println(err)
            false
        end
    end
end

const u01 = [.3, .4, .2]
const p1 = [.2, .33, .1, .4]
const tspan = (0., 70.)
du01 = deepcopy(u01)
@test @isinferred sir_26!(du01, u01, p1, tspan)
const prob1 = ODEProblem(sir_26!, u01, tspan, p1)
const sol1 = solve(prob1, abstol = 1e-6, reltol = 1e-6)
@test minimum(sol1) >= -1e-5 # allowing some tolerance on this function
@test maximum(sol1) <= 1.00001
const u02 = [.1, .1, .7] 
const p2 = [.9, .1, .7, .1] 
du02 = deepcopy(u02)
@test @isinferred sir_26!(du02, u02, p2, tspan)
const prob2 = ODEProblem(sir_26!, u02, tspan, p2)
const sol2 = solve(prob2, abstol = 1e-6, reltol = 1e-6)
@test minimum(sol2) >= -1e-5
@test maximum(sol2) <= 1.00001

# not using @test @isinferred run_sir_21()
# as currently know to return type Any
const sol3 = run_sir_26() 
@test minimum(sol3) >= -1e-5
@test maximum(sol3) <= 1.00001
const sol4 = run_sir_26(beta = .9, gamma = .1, mu = .7, sigma = .1, 
    S0 = .1, E0 = .1, I0 = .7, duration = 70)
@test minimum(sol4) >= -1e-5
@test maximum(sol4) <= 1.00001
@test last(sol4) ≈ last(sol2)

# these bottom tests are only to make sure that the functions run without an error 
@test @noerrors print_sir_26()
@test @noerrors plot_sir_26() 
