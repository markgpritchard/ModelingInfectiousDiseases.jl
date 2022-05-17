
module ModelingInfectiousDiseases

using DifferentialEquations, Reexport

@reexport using DifferentialEquations

# Programmes are grouped by chapter. Each chapter includes a @reexport for each 
# component programme.

include("chapter2/Chapter2.jl")

end # module ModelingInfectiousDiseases
