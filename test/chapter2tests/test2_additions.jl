
using ModelingInfectiousDiseases.Chapter2Additions
using Test

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

@test @isinferred plot_sirs_23_24() 
