
#using ModelingInfectiousDiseases
using Test, SafeTestsets

@safetestset "Chapter2" begin 
    @safetestset "Programme 2.1" begin include("test2_1.jl") end
    @safetestset "Programme 2.2" begin include("test2_2.jl") end
end 

# To test: Pkg.test("ModelingInfectiousDiseases")
