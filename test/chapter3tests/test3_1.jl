
using ModelingInfectiousDiseases.MID_3_1
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

const u01 = [.3, .2, .3, .2]
const p1 = Parameters_31([5 .2; .2 1], .8)
const tspan = (0., 20.)
du01 = deepcopy(u01)
@test @isinferred sir_31!(du01, u01, p1, tspan)
const prob1 = ODEProblem(sir_31!, u01, tspan, p1)
const sol1 = solve(prob1)
@test sum( last(sol1) ) ≈ 1
@test minimum(sol1) >= 0
@test maximum(sol1) <= 1
const u02 = [.1, .1, .0, .8] 
const p2 = Parameters_31([.1 .2; .2 1], 2)
du02 = deepcopy(u02)
@test @isinferred sir_31!(du02, u02, p2, tspan)
const prob2 = ODEProblem(sir_31!, u02, tspan, p2)
const sol2 = solve(prob2)
@test sum( last(sol2) ) ≈ 1
@test minimum(sol2) >= 0
@test maximum(sol2) <= 1

# not using @test @isinferred run_sir_31()
# as currently know to return type Any
const sol3 = run_sir_31() 
@test sum( last(sol3) ) ≈ 1
@test minimum(sol3) >= 0
@test maximum(sol3) <= 1
const sol4 = run_sir_31(beta = [.1 .2; .2 1], gamma = 2, nh = .2, Ih = .1, Il = .8, duration = 20)
@test sum( last(sol4) ) ≈ 1
@test minimum(sol4) >= 0
@test maximum(sol4) <= 1
@test last(sol4) ≈ last(sol2)

# these bottom tests are only to make sure that the functions run without an error 
@test @noerrors print_sir_31()
@test @noerrors plot_sir_31() 
