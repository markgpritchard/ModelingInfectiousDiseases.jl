
using Test, SafeTestsets

@safetestset "Chapter 2" begin include("chapter2tests/testchapter2.jl") end 
@safetestset "Chapter 3" begin include("chapter3tests/testchapter3.jl") end

# To test: using ModelingInfectiousDiseases, Pkg; Pkg.test("ModelingInfectiousDiseases")
