
using ModelingInfectiousDiseases.MID_2_7
using Test
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

const u01 = [.3, .4, .2]
const p1 = [.2, .01, .001, .1, .001, .3]
const tspan = (0., 70.)
du01 = deepcopy(u01)
@test @isinferred sir_27!(du01, u01, p1, tspan)
const prob1 = ODEProblem(sir_27!, u01, tspan, p1)
const sol1 = solve(prob1, abstol = 1e-6, reltol = 1e-6)
@test minimum(sol1) >= 0
@test maximum(sol1) <= 1
const u02 = [.1, .1, .7] 
const p2 = [.6, .01, .001, .1, .001, .4] 
du02 = deepcopy(u02)
@test @isinferred sir_27!(du02, u02, p2, tspan)
const prob2 = ODEProblem(sir_27!, u02, tspan, p2)
const sol2 = solve(prob2)
@test minimum(sol2) >= 0
@test maximum(sol2) <= 1

# not using @test @isinferred run_sir_21()
# as currently know to return type Any
const sol3 = run_sir_27() 
@test minimum(sol3) >= 0
@test maximum(sol3) <= 1
const sol4 = run_sir_27(beta = .6, gamma = .01, Gamma = .001, epsilon = .1, mu = .001, q = .4, 
    S0 = .1, I0 = .1, C0 = .7, duration = 70)
@test minimum(sol4) >= 0
@test maximum(sol4) <= 1
@test last(sol4) ≈ last(sol2)

# no tests of print_sir_21

@test @isinferred plot_sir_27() 
