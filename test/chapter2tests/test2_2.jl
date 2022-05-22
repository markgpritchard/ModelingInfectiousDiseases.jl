
using ModelingInfectiousDiseases.MID_2_2
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

const u01 = [.3, .4, .3]
const p1 = [.2, .33, 1/12_000]
const tspan = (0., 70.)
du01 = deepcopy(u01)
@test @isinferred sir_22!(du01, u01, p1, tspan)
const prob1 = ODEProblem(sir_22!, u01, tspan, p1)
const sol1 = solve(prob1)
@test sum( last(sol1) ) ≈ 1
@test minimum(sol1) >= 0
@test maximum(sol1) <= 1
const u02 = [.1, .1, .8] 
const p2 = [.9, .1, 1e-6] 
du02 = deepcopy(u02)
@test @isinferred sir_22!(du02, u02, p2, tspan)
const prob2 = ODEProblem(sir_22!, u02, tspan, p2)
const sol2 = solve(prob2; abstol = 1e-12, reltol = 1e-12)
@test sum( last(sol2) ) ≈ 1
@test minimum(sol2) >= 0
@test maximum(sol2) <= 1

# not using @test @isinferred run_sir_22()
# as currently know to return type Any
const sol3 = run_sir_22() 
@test sum( last(sol3) ) ≈ 1
@test minimum(sol3) >= 0
@test maximum(sol3) <= 1
const sol4 = run_sir_22(beta = .9, gamma = .1, mu = 1e-6, S0 = .1, I0 = .1, duration = 70)
@test sum( last(sol4) ) ≈ 1
@test minimum(sol4) >= 0
@test maximum(sol4) <= 1
@test last(sol4) ≈ last(sol2)

# these bottom tests are only to make sure that the functions run without an error 
@test @noerrors print_sir_22()
@test @noerrors plot_sir_22() 
