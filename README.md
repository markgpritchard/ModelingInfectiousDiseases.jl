# ModelingInfectiousDiseases.jl

This code recreates the programmes from http://www.modelinginfectiousdiseases.org/ into Julia. These are from the book Keeling, M.J. and Rohani, P., *Modeling Infectious Diseases in Humans and Animals*, Princeton University Press (2007). The code is being developed as a learning exercise. It is not associated with the authors or publishers of the book, and is saved here in 'as is' form.

## Use

The code has been written in Julia, version 1.8.2. 

To find your current version, use the simple command, 
``` julia 
julia> VERSION
v"1.8.2"
```
If you do not have version 1.8.2 or newer, you can [download it now](https://julialang.org/downloads/). 

You are recommended to _clone_ this git repository (see [instructions here](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository)). Then open the file `RunModelingInfectiousDiseases.jl`. You should replace the string on line 11 with the path to your local copy,

``` julia 
###############################################################################
# Enter the path to this file here 
loc = "C:\\Users\\username\\Documents\\GitHub\\ModelingInfectiousDiseases.jl"
###############################################################################
```

The following lines of that file allow you to activate and instantiate the project so you will be using the same versions of each package as the programmes were written with. (The first time you do this it may take a while, depending on how many packages need to be downloaded and precompiled.)

None of the programmes is intended or expected to be the quickest or most efficient solution. Rather, they are intended to make clear how the model can be coded in *Julia*.

## Layout of code

It is intended that you will be able to run all the programmes from within the main `RunModelingInfectiousDiseases.jl` file. Example parameters and initial conditions are supplied for each model, and you can explore the effects of changing these.

The code for each programme is located in the `src` folder, subdivided by chapter. You can read (and modify) the code there. Each programme is stored in its own module so any functions or constants defined for one programme will not interfere with any others. 

Within each file, the code is laid out with any structures, types and constants listed first. Then is the main function for the model (in most cases, this is the function that will be passed to `DifferentialEquations.solve`). We then provide a function that will run the model. In each case the first version of this function takes each initial condition, parameter and duration as a keyword argument which makes it clear exactly which input is which. Functions that support this are then provided. The programmes generally finish with functions to convert the model output into a dataframe, to plot the output, and where applicable (Chapter 7) to produce a video output.

No formal help text is provided. Comments within the functions and in `RunModelingInfectiousDiseases.jl` should make clear what each function does. Dividing the code and parameters into different files allows them to be viewed side-by-side in a split screen, and also makes clear the difference between changing the input and changing the model.
