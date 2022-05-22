
module ModelingInfectiousDiseases

using DifferentialEquations, Reexport

@reexport using DifferentialEquations

include("chapter2/Chapter2.jl")

using   # each programme is in its own module  
    # chapter 2:
    .MID_2_1, .MID_2_2, .MID_2_3, .MID_2_4, .MID_2_5, .MID_2_6, .MID_2_7, .Chapter2Additions 

@reexport using 
    .MID_2_1, .MID_2_2, .MID_2_3, .MID_2_4, .MID_2_5, .MID_2_6, .MID_2_7, .Chapter2Additions 

end # module ModelingInfectiousDiseases
