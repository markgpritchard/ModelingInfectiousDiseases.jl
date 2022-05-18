# ModelingInfectiousDiseases.jl

This code recreates the programmes from http://www.modelinginfectiousdiseases.org/ into Julia. These are from the book Keeling, M.J. and Rohani, P., *Modeling Infectious Diseases in Humans and Animals*, Princeton University Press (2007). The code is being developed as a learning exercise. It is not associated with the authors or publishers of the book, and is saved here in 'as is' form.

The code is presented in a package. You may install the whole package and use functions from it, or copy individual programmes from the `src` folder. Functions in that folder are grouped by book chapter.  

Plots of model outputs are produced within this package in `CairoMakie`. Numerical outputs are available for the user to produce plots using any package of their choice. 

So far, programmes up to the end of chapter 2 (programme 2.7) have been added.

## Installation 
``` julia 
julia> ]
pkg> add https://github.com/markgpritchard/ModelingInfectiousDiseases.jl
julia> using ModelingInfectiousDiseases
```

## Use
###Chapter 2
Each programme takes the same form, with the last two digits indicating the programme number. For example, for Programme 2.2, the function `sir_22!` holds the ordinary differential equations. `run_sir_21` runs the model, `print_sir_21` displays the output in the REPL, and `plot_sir_21` plots the output.

Programmes 2.3 and 2.4 both have pathogen-related mortality. 2.3 uses density-dependent transmission and 2.4 uses frequency dependent transmission. Their outputs can be seen side-by-side using the function `plot_sirs_23_24`.

Documentation is not yet complete but each function has help text available:
``` julia 
julia> ?
help?> run_sir_21
search: run_sir_21 run_sir_27 run_sir_26 run_sir_25 run_sir_24 run_sir_23 run_sir_22

  run_sir_21([; beta, gamma, S0, I0, duration, saveat])

  Run the model sir_21

  Keyword arguments
  ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

  All keyword arguments are optional with default values supplied for each.

    •  `beta`: The beta parameter in the model (infectiousness of infectives). Default is 520 / 365
       (520 per year).

    •  gamma: The gamma parameter in the model (recovery rate). Default is 1 / 7.

    •  S0: Proportion of the population susceptible at time = 0. Default is 1 - 1e-6.

    •  I0: Proportion of the population infectious at time = 0. Default is 1e-6.

    •  R_at_time0: Proportion of the population resistant at time = 0. Has a long name to avoid
       confusion with the basic reproduction number R₀. Default value is 0.

    •  duration: How long the model will run (time units are interpretted as days). Default is 70.

    •  saveat: How frequently the model should save values. Default is 1 (day).

  Examples
  ≡≡≡≡≡≡≡≡≡≡

  julia> run_sir_21()
  
  julia> run_sir_21(beta = .8, gamma = .6, duration = 1000)
```
