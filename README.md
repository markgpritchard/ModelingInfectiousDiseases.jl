# ModelingInfectiousDiseases.jl

This code recreates the programmes from http://www.modelinginfectiousdiseases.org/ into Julia. These are from the book Keeling, M.J. and Rohani, P., *Modeling Infectious Diseases in Humans and Animals*, Princeton University Press (2007). The code is being developed as a learning exercise. It is not associated with the authors or publishers of the book, and is saved here in 'as is' form.

The code is presented in a package. You may install the whole package and use functions from it, or copy individual programmes from the `src` folder. Functions in that folder are grouped by book chapter.  

Plots of model outputs are produced within this package in `CairoMakie`. Numerical outputs are available for the user to produce plots using any package of their choice. 

So far, programmes up to the end of chapter 2 (programme 2.7) have been added.

Documentation is not yet complete but each function has help text available, which can be accessed with `?` then the function name.

## Installation 
``` julia 
julia> ]
pkg> add https://github.com/markgpritchard/ModelingInfectiousDiseases.jl
julia> using ModelingInfectiousDiseases
```

## Use
### Chapter 2
Each programme takes the same form, with the last two digits indicating the programme number. For example, for Programme 2.2, the function `sir_22!` holds the ordinary differential equations, `run_sir_21` runs the model, `print_sir_21` displays the output in the REPL, and `plot_sir_21` plots the output. (Note that for simplicity, all functions are named `sir` even if they are not strictly SIR [susceptible--infectious--recovered] models.)

Programmes 2.3 and 2.4 both present pathogen-related mortality. 2.3 uses density-dependent transmission and 2.4 uses frequency dependent transmission. Their outputs can be seen side-by-side using the function `plot_sirs_23_24`.

![Outputs of programmes 2.3 and 2.4 side-by-side, showing damped oscilation toward an equilibrium for numbers infectious, but a difference in numbers susceptible, with density-dependent transmission showing an equilibrium and frequency-dependent transmission showing a moving trend](https://github.com/markgpritchard/ModelingInfectiousDiseases.jl/blob/main/assets/plt_2324.png)
