
using Test
using ModelingInfectiousDiseases.MID_2_2
using DifferentialEquations, StaticArrays

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

# not with a static array
const u01 = [.3, .4, .3]
const p1 = [.2, .33, 1/12_000]
const tspan = (0., 70.)
@test @isinferred sir_22(u01, p1, tspan)
prob1 = ODEProblem(sir_22, u01, tspan, p1)
sol1 = solve(prob1)
@test sum( last(sol1) ) ≈ 1
# with a static array
const u02 = @SVector [.1, .1, .8] 
const p2 = [.9, .1, 1e-6] 
@test @isinferred sir_22(u02, p2, tspan)
prob2 = ODEProblem(sir_22, u02, tspan, p2)
sol2 = solve(prob2; abstol = 1e-12, reltol = 1e-12)
@test sum( last(sol2) ) ≈ 1

# not using @test @isinferred run_sir_22()
# as currently know to return type Any
sol3 = run_sir_22() 
@test sum( last(sol3) ) ≈ 1
sol4 = run_sir_22(beta = .9, gamma = .1, mu = 1e-6, S0 = .1, I0 = .1, duration = 70)
@test last(sol4) ≈ last(sol2)

# no tests of print_sir_22

@test @isinferred plot_sir_22() 
