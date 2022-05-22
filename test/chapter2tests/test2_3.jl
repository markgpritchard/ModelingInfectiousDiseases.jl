
using ModelingInfectiousDiseases.MID_2_3
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

# N0 = 1
const u01 = [.3, .4, .3]
const p1 = [.9, .2, 1/12_000, 1/10_000, .4]
const tspan = (0., 70.)
du01 = deepcopy(u01)
@test @isinferred sir_23!(du01, u01, p1, tspan)
const prob1 = ODEProblem(sir_23!, u01, tspan, p1)
const sol1 = solve(prob1)
@test minimum(sol1) >= 0
# N0 > 1; u0 all integer
const u02 = [100_000, 10_000, 390_000] 
const p2 = [.33, .1, 1/12_000, 1/10_000, .8]
du02 = deepcopy(u01) # du will never be integer so not testing that here
@test @isinferred sir_23!(du02, u02, p2, tspan)
const prob2 = ODEProblem(sir_23!, u02, tspan, p2)
const sol2 = solve(prob2; abstol = 1e-12, reltol = 1e-12)
@test minimum(sol2) >= 0

# not using @test @isinferred run_sir_23()
# as currently know to return type Any
const sol3 = run_sir_23() 
@test minimum(sol3) >= 0
const sol4 = run_sir_23(beta = .33, gamma = .1, mu = 1/12_000, nu = 1/10_000, rho = .8, 
    X0 = 100_000, Y0 = 10_000, N0 = 500_000, duration = 70)
@test minimum(sol4) >= 0
@test last(sol4) ≈ last(sol2)

# these bottom tests are only to make sure that the functions run without an error 
@test @noerrors print_sir_23()
@test @noerrors plot_sir_23() 
