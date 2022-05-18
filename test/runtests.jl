
using Test, SafeTestsets

@safetestset "Chapter 2" begin include("chapter2tests/testchapter2.jl") end 

# To test: using ModelingInfectiousDiseases, Pkg; Pkg.test("ModelingInfectiousDiseases")
