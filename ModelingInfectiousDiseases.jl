
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set up 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Run all the code in this set up section to ensure you have the packages expected 
# by the programmes and your working directory includes all the required code

###############################################################################
# Enter the path to this file here 
loc = "C:\\Users\\yourname\\Documents\\GitHub\\ModelingInfectiousDiseases.jl"
###############################################################################

cd(loc)

using Pkg 
Pkg.activate(loc)
Pkg.instantiate()

using CairoMakie

# There is an `include` statement for each programme so you do not need to wait 
# for Julia to read code that you will not be using in this session


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Chapter 2 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Programme 2.1

include("src/chapter2/p2.1.jl")
using .MID_21

beta = 520 / 365        # infectiousness parameter
gamma = 1 / 7           # recovery rate
S0 = 1 - 1e-6           # initial proportion susceptible 
I0 = 1e-6               # initial proportion resistant
duration = 70           # duration of model

sol21 = run_sir21(; beta, gamma, S0, I0, duration, saveat = .125) # frequent saveat to give smooth plot
plot_sir21(sol21)


## Programme 2.2

include("src/chapter2/p2.2.jl")
using .MID_22

beta = 520 / 365        # infectiousness parameter 
gamma = 1 / 7           # recovery rate 
mu = 1 / (70 * 365)     # birth and mortality rate 
S0 = .1                 # initial proportion susceptible
I0 = 1e-4               # initial proportion infectious
duration = 60 * 365     # duration of model

sol22 = run_sir22(; beta, gamma, mu, S0, I0, duration)
plot_sir22(sol22)


## Programme 2.3

include("src/chapter2/p2.3.jl")
using .MID_23

beta = 520 / 365        # infectiousness parameter 
gamma = 1 / 7           # recovery rate 
mu = 1 / (70 * 365)     # mortality rate not due to the pathogen
nu = 1 / (70 * 365)     # birth rate 
rho = .5                # mortality probability for infecteds 
X0 = .2                 # initial number susceptible 
Y0 = 1e-6               # initial number infectious 
N0 = 1                  # initial population (NB when N0 = 1 then X, Y and Z are proportions)
duration = 100 * 365    # duration of model

sol23 = run_sir23(; beta, gamma, mu, nu, rho, X0, Y0, N0, duration)
plot_sir23(sol23)


## Programme 2.4

include("src/chapter2/p2.4.jl")
using .MID_24

beta = 520 / 365        # infectiousness parameter 
gamma = 1 / 7           # recovery rate 
mu = 1 / (70 * 365)     # mortality rate not due to the pathogen
nu = 1 / (70 * 365)     # birth rate 
rho = .5                # mortality probability for infecteds
X0 = .2                 # initial number susceptible 
Y0 = 1e-6               # initial number infectious 
N0 = 1                  # initial population (NB when N0 = 1 then X, Y and Z are proportions)
duration = 1e5          # duration of model

sol24 = run_sir24(; beta, gamma, mu, nu, rho = .5, X0 = .2, Y0, N0, duration)
plot_sir24(sol24)


## Additional function to view outputs of programme 2.3 and 2.4 side-by-side 

include("src/chapter2/additions.jl")
using .Chapter2Additions

beta = 520 / 365        # infectiousness parameter 
gamma = 1 / 7           # recovery rate 
mu = 1 / (70 * 365)     # mortality rate not due to the pathogen
nu = 1 / (70 * 365)     # birth rate 
rho = .5                # mortality probability for infecteds
X0 = .2                 # initial number susceptible 
Y0 = 1e-6               # initial number infectious 
N0 = 1                  # initial population (NB when N0 = 1 then X, Y and Z are proportions)
duration = 100 * 365    # duration of model

plot_sirs23_24(; beta, gamma, mu, nu, rho, X0, Y0, N0, duration)


## Programme 2.5

include("src/chapter2/p2.5.jl")
using .MID_25

beta = 520 / 365        # infectiousness parameter 
gamma = 1 / 7           # recovery rate 
I0 = 1e-6               # initial proportion infectious
duration = 70           # duration of model

sol25 = run_sis25(; beta, gamma, I0, duration, saveat = .125) # frequent saveat to give a smooth plot
plot_sis25(sol25)


## Programme 2.6 

include("src/chapter2/p2.6.jl")
using .MID_26

beta = 520 / 365        # infectiousness parameter 
gamma = 1 / 7           # recovery rate 
mu = 1 / (70 * 365)     # birth and mortality rate 
sigma = 1 / 14          # rate at which exposed individuals become infectious
S0 = .1                 # initial proportion susceptible
E0 = 1e-4               # initial proportion exposed
I0 = 1e-4               # initial proportion infectious
duration = 60 * 365     # duration of model

sol26 = run_seir26(; beta, gamma, mu, sigma, S0, E0, I0, duration)
plot_seir26(sol26)


## Programme 2.7 

include("src/chapter2/p2.7.jl")
using .MID_27

beta = 0.2              # infectiousness parameter 
gamma = 1 / 100         # recovery rate from infectious 
Gamma = 1 / 1000        # recovery rate of carriers 
epsilon = 0.1           # proportion reduction in transmission from carriers compared to infecteds 
mu = 1 / (50 * 365)     # birth and mortality rate 
q = .4                  # proportion of infected who become carriers
S0 = .1                 # initial proportion susceptible
I0 = 1e-4               # initial proportion infectious
C0 = 1e-3               # initial proportion carriers
duration = 60 * 365     # duration of model

sol27 = run_sir27(; beta, gamma, Gamma, epsilon, mu, q, S0, I0, C0, duration)
plot_sir27(sol27)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Chapter 3
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Programme 3.1

include("src/chapter3/p3.1.jl")
using .MID_31

beta = [10 .1; .1 1]    # matrix of infectiousness parameters [β_hh, β_hl; β_lh β_ll]
gamma = 1               # recovery rate 
nh = .2                 # proportion of high-risk individuals 
Ih = 1e-5               # initial proportion of infectious high-risk individuals
Il = 1e-3               # initial proportion of infectious low-risk individuals
duration = 15           # duration of model

sol31 = run_sir31(; beta, gamma, nh, Ih, Il, duration, saveat = .025) # frequent saveat to give a smooth plot
plot_sir31(sol31)


## Programme 3.2

include("src/chapter3/p3.2.jl")
using .MID_32 

# For this example we have five risk groups. 

u0 = [                  # Initial conditions for the model
    .06 .31 .52 .08 .02999  # susceptibles 
    .0  .0  .0  .0  1e-5    # infectious
]

# We use an intermediary step when creating the transmission factor so that all βᵢⱼ == βⱼᵢ 
betavector = [ 0, 3, 10, 60, 100]
betamatrix = .0016 .* betavector * betavector'
p = Parameters32(       # Model parameters
    betamatrix,             # beta = matrix of infectiousness parameters 
    [.2, .2, .2, .2, .2]    # gamma = vector of recovery rates
)

duration = 30           # Duration

# u0 must sum to 1, but because there are many small values there may be rounding 
# errors, which we allow for
sum(u0) ≈ 1  

sol32 = run_sis32(u0, p, duration; saveat = .025)
result32 = dataframe_sis32(sol32; type = :both)
plot_sis32(sol32; legend = :below)


## Programme 3.3

include("src/chapter3/p3.3.jl")
using .MID_33

u0 = [                  # Initial conditions for the model
    .1      .1              # susceptibles 
    .0001   .0001           # infectious
    .0999   .6999           # recovered
]
p = Parameters33(       # Model parameters
    [   100.    10.         # beta = matrix of infectiousness parameters 
        10.     20. ],             
    10.,                    # gamma = recovery rate
    1 / 15,                 # lambda = rate that children become adults 
    [0., 1 / 60],           # mu = vector of mortality rates
    1 / 60                  # nu = birth rate
)
duration = 100          # Duration

sol33 = run_sir33(u0, p, duration; saveat = .001)
result33 = dataframe_sir33(sol33)
plot_sir33(sol33; legend = :below)


## Programme 3.4

include("src/chapter3/p3.4.jl")
using .MID_34

u0 = [                  # Initial conditions for the model
    .05     .01         .01         .008        # susceptibles 
    .0001   .0001       .0001       .0001       # exposed
    .0001   .0001       .0001       .0001       # infectious
    .0298   .04313333   .12313333   .72513333   # recovered
]
p = Parameters34(       # Model parameters
    [2.089 2.089 2.086 2.037    # beta = matrix of infectiousness parameters 
    2.089 9.336 2.086 2.037
    2.086 2.086 2.086 2.037
    2.037 2.037 2.037 2.037],
    1 / 8,                      # sigma = rate of movement into infectious class 
    1 / 5,                      # gamma = recovery rate 
    1 / (55 * 365),             # mu = mortality rates 
    1 / (55 * 365)              # nu = birth rate
)
duration = 100 * 365    # Duration

result34 = run_seir34(u0, p, duration)
plot_seir34(result34; legend = :below)


## Programme 3.5

include("src/chapter3/p3.5.jl")
using .MID_35

p_1 = Parameters35(     # Model parameters
    17 / 5,                 # beta = infectiousness parameter 
    1 / 13,                 # sigma = rate of movement into infectious class 
    1 / 13,                 # gamma = recovery rate 
    1 / (55 * 365),         # mu = mortality rates 
    1 / (55 * 365),         # nu = birth rate
    8,                      # m = number of exposed compartments 
    13                      # n = total number of infected compartments (E + I)
)
u0_1 = seir35_u0(.05, 0, .00001, p_1)   # Initial conditions for the model
duration_1 = 30 * 365   # Duration

sol35_1 = run_seir35(u0_1, p_1, duration_1)
result35_1 = dataframe_seir35(sol35_1, p_1)
plot_seir35(result35_1; legend = :below)

### Alternative sets of parameters

p_2 = Parameters35(     # Model parameters
    1.,                     # beta = infectiousness parameter 
    .0,                     # sigma = rate of movement into infectious class 
    .1,                     # gamma = recovery rate 
    .0,                     # mu = mortality rates 
    .0,                     # nu = birth rate
    0,                      # m = number of exposed compartments 
    10                      # n = total number of infected compartments (E + I)
)
u0_2 = seir35_u0(.5, 0, 1e-6, p_2)  # Initial conditions for the model
duration_2 = 60         # Duration

sol35_2 = run_seir35(u0_2, p_2, duration_2; saveat = .1)
result35_2 = dataframe_seir35(sol35_2, p_2)

p_3 = Parameters35(     # Model parameters
    1.,                     # beta = infectiousness parameter 
    .0,                     # sigma = rate of movement into infectious class 
    .1,                     # gamma = recovery rate 
    .0,                     # mu = mortality rates 
    .0,                     # nu = birth rate
    0,                      # m = number of exposed compartments 
    1                       # n = total number of infected compartments (E + I)
)
u0_3 = seir35_u0(.5, 0, 1e-6, p_3)  # Initial conditions for the model
duration_3 = 60         # Duration

sol35_3 = run_seir35(u0_3, p_3, duration_3; saveat = .1)
result35_3 = dataframe_seir35(sol35_3, p_3)

fig35_1 = Figure()
ga = GridLayout(fig35_1[1, 1])
plot_seir35!(ga, result35_2; label = "SIR with 10 I compartments", legend = :none)
gb = GridLayout(fig35_1[1, 2])
plot_seir35!(gb, result35_3; label = "SIR with 1 I compartment", legend = :right)
fig35_1

p_4 = Parameters35(     # Model parameters
    1.,                     # beta = infectiousness parameter 
    .1,                     # sigma = rate of movement into infectious class 
    .1,                     # gamma = recovery rate 
    .0,                     # mu = mortality rates 
    .0,                     # nu = birth rate
    5,                      # m = number of exposed compartments 
    10                      # n = total number of infected compartments (E + I)
)
u0_4 = seir35_u0(.5, 0, 1e-4, p_4)  # Initial conditions for the model
duration_4 = 150        # Duration

sol35_4 = run_seir35(u0_4, p_4, duration_4; saveat = .2)
result35_4 = dataframe_seir35(sol35_4, p_4)

p_5 = Parameters35(     # Model parameters
    1.,                     # beta = infectiousness parameter 
    .1,                     # sigma = rate of movement into infectious class 
    .1,                     # gamma = recovery rate 
    .0,                     # mu = mortality rates 
    .0,                     # nu = birth rate
    1,                      # m = number of exposed compartments 
    2                       # n = total number of infected compartments (E + I)
)
u0_5 = seir35_u0(.5, 0, 1e-4, p_5)  # Initial conditions for the model
duration_5 = 150        # Duration

sol35_5 = run_seir35(u0_5, p_5, duration_5; saveat = .2)
result35_5 = dataframe_seir35(sol35_5, p_5)

fig35_2 = Figure()
ga = GridLayout(fig35_2[1, 1])
plot_seir35!(ga, result35_4; label = "SEIR with 5 E and 5 I", legend = :none)
gb = GridLayout(fig35_2[1, 2])
plot_seir35!(gb, result35_5; label = "SEIR with 1 E and 1 I", legend = :right)
fig35_2


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Chapter 6 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Programme 6.1 

include("src/chapter6/p6.1.jl")
using .MID_61

u0 = [                  # 
    1e5,                # X0 -- initial number susceptible
    500,                # Y0 -- initial number infectious
    1e6 - (1e5 + 500)   # Z0 -- initial number recovered
]

### Run with no noise

p_nonoise = Parameters61(# Model parameters
    1.,                 # beta -- infection parameter 
    .1,                 # gamma -- recovery rate 
    1 / (50 * 365),     # mu -- birth rate 
    1 / (50 * 365),     # nu -- death rate (recommended to equal birth rate)
    0.                  # xi -- magnitude of the noise that will be added 
)

results61_nonoise = run_sir61(u0, p_nonoise, 5 * 365)
plot_sir61(results61_nonoise, p_nonoise)

### Run with noise parameter

p = Parameters61(
    1.,                 # beta -- infection parameter  
    .1,                 # gamma -- recovery rate 
    1 / (50 * 365),     # mu -- birth rate 
    1 / (50 * 365),     # nu -- death rate 
    10.                 # xi -- magnitude of the noise that will be added 
)

results61 = run_sir61(u0, p, 5 * 365; seed = 61)
plot_sir61(results61, p)


## Programme 6.2 

include("src/chapter6/p6.2.jl")
using .MID_62

u0 = [
    1e5,                # X0 -- initial number susceptible
    500,                # Y0 -- initial number infectious
    1e6 - (1e5 + 500)   # Z0 -- initial number recovered
]

### Run with no noise

p_nonoise = Parameters62(
    1.,                 # beta -- infection parameter 
    .1,                 # gamma -- recovery rate 
    1 / (50 * 365),     # mu -- birth rate 
    1 / (50 * 365),     # nu -- death rate (recommended to equal birth rate)
    0.                  # xi -- magnitude of the noise that will be added 
)

results62_nonoise = run_sir62(u0, p_nonoise, 5 * 365)
plot_sir62(results62_nonoise, p_nonoise)

### Run with noise parameter

p = Parameters62(
    1.,                 # beta -- infection parameter 
    .1,                 # gamma -- recovery rate 
    1 / (50 * 365),     # mu -- birth rate 
    1 / (50 * 365),     # nu -- death rate
    1.                  # xi -- magnitude of the noise that will be added 
)

results62 = run_sir62(u0, p, 5 * 365; seed = 62)
plot_sir62(results62, p)

### Run with a large noise parameter

p_bignoise = Parameters62(
    1.,                 # beta -- infection parameter 
    .1,                 # gamma -- recovery rate 
    1 / (50 * 365),     # mu -- birth rate 
    1 / (50 * 365),     # nu -- death rate
    10.                 # xi -- magnitude of the noise that will be added 
)

results62_bignoise = run_sir62(u0, p_bignoise, 5 * 365; seed = 62)
plot_sir62(results62_bignoise, p_bignoise)


## Programme 6.3 

include("src/chapter6/p6.3.jl")
using .MID_63

u0 = [
    30,                 # X0 -- initial number susceptible 
    70                  # Y0 -- initial number infectious
]
p = [
    .03,                # beta -- infection parameter  
    .01                 # gamma -- recovery rate  
]

results63 = run_sis63(u0, p, 10 * 365; seed = 63)
plot_sis63(results63)


## Programme 6.4 

include("src/chapter6/p6.4.jl")
using .MID_64

p = [
    1.,                 # beta -- infection parameter  
    .1,                 # gamma -- recovery rate  
    5e-4                # mu -- birth and death rate
]

### Examine model with small population

u0_50 = u0_sir64(50, p)
results64_50 = run_sir64(u0_50, p, 2 * 365; seed = 64)
plot_sir64(results64_50, 50)

### and with a larger population

u0 = u0_sir64(5000, p)
results64 = run_sir64(u0, p, 2 * 365; seed = 64)
plot_sir64(results64, 5000)


## Programme 6.5

include("src/chapter6/p6.5.jl")
using .MID_65

p = [
    1.,                 # beta -- infection parameter  
    .1,                 # gamma -- recovery rate  
    5e-4                # mu -- birth and death rate
]

### Examine model with small population

u0_50 = u0_sir65(50, p)
results65_50 = run_sir65(u0_50, p, 2 * 365; seed = 65)
plot_sir65(results65_50, 50)

### and with a larger population

u0 = u0_sir65(5000, p)
results65 = run_sir65(u0, p, 2 * 365; seed = 65)
plot_sir65(results65, 5000)


## Programme 6.6

include("src/chapter6/p6.6.jl")
using .MID_66

### Examine model with small population

# Recommended that ε parameter is adjusted with inverse population size
p_50 = [
    1.,                 # beta -- infection parameter  
    .1,                 # gamma -- recovery rate  
    .0001,              # delta -- rate of infectious immigration 
    .002,               # epsilon -- force of external infection
    5e-4                # mu -- birth and death rate
]
u0_50 = u0_sir66(50, p_50)
results66_50 = run_sir66(u0_50, p_50, 2 * 365; seed = 66)
plot_sir66(results66_50, 50)

### and with a larger population

p = [
    1.,                 # beta -- infection parameter  
    .1,                 # gamma -- recovery rate  
    .01,                # delta -- rate of infectious immigration 
    .00002,             # epsilon -- force of external infection
    5e-4                # mu -- birth and death rate
]
u0 = u0_sir66(5000, p)
results66 = run_sir66(u0, p, 2 * 365; seed = 66)
plot_sir66(results66, 5000)

# This plot looks interesting -- what happens over 10 years?
results66_10y = run_sir66(u0, p, 3650; seed = 66)
plot_sir66(results66_10y, 5000)

# And what happens if we change δt? 

results66_10y_d10 = run_sir66(u0, p, 10 * 365; seed = 66, δt = 10)
plot_sir66(results66_10y_d10, "p6.6.jl: SIR model with τ-leap stochasticity\nδt = 10")

results66_10y_d01 = run_sir66(u0, p, 10 * 365; seed = 66, δt = .1)
plot_sir66(results66_10y_d01, "p6.6.jl: SIR model with τ-leap stochasticity\nδt = 0.1")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Chapter 7 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Programme 7.1

include("src/chapter7/p7.1.jl")
using .MID_71

u0 = [                  # Initial conditions for the model
    .1      .1  .1  .1  .1  # susceptibles 
    .0001   .0  .0  .0  .0  # infectious
    .0      .0  .0  .0  .0  # recovered
]
p = Parameters71(       # Model parameters
    ones(5),                # beta = vector of infection parameters
    .1 * ones(5),           # gamma = vector of recovery rates 
    .0001 * ones(5),        # mu = vector of mortality rates
    .0001 * ones(5),        # nu = vector of birth rates
    .001 * ones(5, 5)       # m = matrix of migrations between subpopulations
)
duration = 2910         # Duration

sol71 = run_sir71(u0, p, duration)
result71 = dataframe_sir71(sol71)
plot_sir71(result71) 


## Programme 7.2

include("src/chapter7/p7.2.jl")
using .MID_72

# Initial conditions for the model:
X0 = [
    800 0   0   0   0
    0   800 0   0   0
    0   0   800 0   0
    0   0   0   800 0
    0   0   0   0   800.
]
Y0 = zeros(5, 5); Y0[1, 1] = 1 
Z0 = [
    199 0   0   0   0
    0   200 0   0   0
    0   0   200 0   0
    0   0   0   200 0
    0   0   0   0   200.
]
u0 = zeros(5, 5, 3) 
u0[:, :, 1] = X0; u0[:, :, 2] = Y0; u0[:, :, 3] = Z0

p = Parameters72(       # Model parameters
    ones(5),                # beta = vector of infection parameters
    .3 * ones(5),           # gamma = vector of recovery rates 
    zeros(5),               # mu = vector of mortality rates
    zeros(5, 5),            # nu = matrix of birth rates for every population--location combination
    zeros(5, 5),            # l = matrix of movements from home subpopulation 
    2 * ones(5, 5)          # r = matrix of movements back to home subpopulation 
)
# Formula for p.l from the Python code 
for i ∈ 1:5, j ∈ 1:5 
    if abs(i - j) == 1 p.l[i, j] = .1 end 
end 

duration = 60           # Duration 

sol72 = run_sir72(u0, p, duration; saveat = .125)
# Note that there are so many compartments in this model we do not display a DataFrame
# of results
plot_sir72(sol72) 

### repeat with more varied parameters 

X0_2 = [
    999 0   0   0   0
    0  1000 0   0   0
    0   0  1000 0   0
    0   0   0  1000 0
    0   0   0   0  1000.
]
Y0_2 = zeros(5, 5); Y0_2[1, 1] = 1 
Z0_2 = zeros(5, 5)
u0_2 = zeros(5, 5, 3) 
u0_2[:, :, 1] = X0_2; u0_2[:, :, 2] = Y0_2; u0_2[:, :, 3] = Z0_2

p_2 = Parameters72(     # Model parameters
    .5 * ones(5),               # beta = vector of infection parameters
    .2 * ones(5),               # gamma = vector of recovery rates 
    1 / (70 * 365) * ones(5),   # mu = vector of mortality rates
    zeros(5, 5),                # nu = matrix of birth rates for every population--location combination
    [                           # l = matrix of movements from home subpopulation 
        0   .15 0   0   0
        .15 0   0   0   0
        0   0   0   0   0 
        0   .15 .15 0   .15 
        0   0   0   0   0
    ],                
    3 * ones(5, 5)              # r = matrix of movements back to home subpopulation 
)
for i ∈ 1:5, j ∈ 1:5 
    if i == j p_2.nu[i, j] = 1000 / (70 * 365) end
end 

duration_2 = 100        # Duration 

sol72_2 = run_sir72(u0_2, p_2, duration_2; saveat = .125)
plot_sir72(sol72_2) 


## Programme 7.3

include("src/chapter7/p7.3.jl")
using .MID_73

u0 = sir73_u0(      # Initial conditions for the model
    25,                 # size of grid 
    .1,                 # value of x0 in each cell 
    4,                  # number of cells with non-zero Y0 
    .001,               # value of Y0 in cells with non-zero Y0
    1;                  # population size in each cell
    seed = 73           # seed for random number generator 
)
p = [               # Model parameters
    1.,                 # beta = infectiousness parameter 
    .1,                 # gamma = recovery rate 
    .0001,              # mu = mortality rate
    .1                  # rho = rate at which individuals interact with neighbouring environments
]
duration = 2910     # Duration

sol73 = run_sir73(u0, p, duration; saveat = 4) # saveat = 4 to give approximately 
    # 30 seconds of video with duration = 2910 and framerate = 24

# This function will save a video in the folder "outputvideos/" as "video73.mp4"
video_sir73(
    sol73; 
    # should the colour scale be constant throughout the video (vs each frame having a separate scale):
    fixmax = true,
    # attempt to find colormap with good differentiation between small values (especially with fixmax = true):
    colormap = :seaborn_colorblind#:gist_stern
)

### Repeat with defined starting points (one infectious individual in the middle and one in a corner)

u0_2 = sir73_u0(    # Initial conditions for the model
    101,                # size of grid 
    999,                # value of x0 in each cell 
    [5101, 9901],       # vector of cells with non-zero Y0 
    1,                  # value of Y0 in cells with non-zero Y0
    1000                # population size in each cell
)
p_2 = [             # Model parameters
    .4,                 # beta = infectiousness parameter 
    .2,                 # gamma = recovery rate 
    4e-5,               # mu = mortality rate
    .1                  # rho = rate at which individuals interact with neighbouring environments
]
duration_2 = 500    # Duration

sol73_2 = run_sir73(u0_2, p_2, duration_2; saveat = 1) 
video_sir73(sol73_2; filename = "video73_2.mp4", fixmax = false)

### A custom addition with a "firebreak" area with resistant population 

u0_3 = sir73_u0(    # Initial conditions for the model
    101,                # size of grid 
    1000,               # value of x0 in each cell 
    0,                  # number of cells with non-zero Y0 
    0,                  # value of Y0 in cells with non-zero Y0
    1000                # population size in each cell
)
for i ∈ axes(u0_3, 1), j ∈ axes(u0_3, 2)
    if j ∈ [50, 51, 52] && i <= 51 
        u0_3[i, j, 1] = 0 
        u0_3[i, j, 3] = 1000 
    elseif j == 53 && i == 1 
        u0_3[i, j, 1] = 999
        u0_3[i, j, 2] = 1
    end 
end 
p_3 = [             # Model parameters
    .4,                 # beta = infectiousness parameter 
    .2,                 # gamma = recovery rate 
    4e-5,               # mu = mortality rate
    .1                  # rho = rate at which individuals interact with neighbouring environments
]
duration_3 = 750    # Duration

sol73_3 = run_sir73(u0_3, p_3, duration_3; saveat = 1) 
video_sir73(sol73_3; filename = "video73_3.mp4", fixmax = true, colormap = :gist_stern)


## Programme 7.4

include("src/chapter7/p7.4.jl")
using .MID_74

u0 = u0_sir74(50)   # Initial conditions for the model (all susceptible)
p = [               # Model parameters
    1.,                 # tau = transmission rate between neighbours
    .1,                 # gamma = recovery rate
    .01,                # nu = birth or waning immunity rate
    .001                # epsilon = import rate
]
duration = 100      # Duration (is 1000 in other example code but 100 allows the 
                    # video to run slower and finish after 100 seconds)

tv, uv = run_sir74(u0, p, duration; seed = 74)
video_sir74(uv, tv)

### Re-run with recovered population 

u0_2 = u0_sir74(50, 0, 2500)  # Initial conditions for the model (all recovered)
p_2 = [             # Model parameters
    1.,                 # tau = transmission rate between neighbours
    .1,                 # gamma = recovery rate
    .01,                # nu = birth or waning immunity rate
    .001                # epsilon = import rate
]
duration_2 = 180    # Duration (longer video as slower to start)

tv_2, uv_2 = run_sir74(u0_2, p_2, duration_2; seed = 742)
video_sir74(uv_2, tv_2; filename = "video74_2.mp4")


## Programme 7.5

include("src/chapter7/p7.5.jl")
using .MID_75

u0 = u0_sis75(      # Initial conditions for the model
    1000,               # n = number of individuals in model  
    4,                  # Y0 = number of initially infectious individuals  
    10;                 # size = size of the grid that the individuals are in 
    seed = 75
) 
p = [               # Model parameters
    3.,                 # alpha = power law decay for the transmission kernal
    .01,                # beta = transmission parameter 
    .5                  # gamma = recovery rate
]
duration = 20       # Duration of the model (example code gives duration of 100)
tstep = .01         # Discrete time intervals used by model 

result75 = run_sis75(u0, p, duration; tstep, seed = 750)

video_sis75(result75; step = 1/48)


## Programme 7.6

include("src/chapter7/p7.6.jl")
using .MID_76

p = Parameters76(   # Model parameters
    [1., 10.5],         # susceptibility parameters for sheep and cows 
    [5.1e-7, 7.7e-7],   # transmissibility parameters for sheep and cows 
    .0                  # diameter of ring culling
)
u0 = u0_seirc76(    # Initial conditions for the model
    4000,               # number of farms
    1,                  # Y0 = initial number of farms with infections
    20,                 # size of grid 
    p;                  # the parameters defined above
    seed = 76
) 
duration = 400      # Duration

result76 = run_seirc76(u0, p, duration; seed = 760)
df76 = dataframe_seirc76(result76)
            
plot_seirc76(df76)
video_seirc76(result76, df76)

### Repeat with ring cull diameter > 0

p_2 = Parameters76( # Model parameters
    [1., 10.5],         # susceptibility parameters for sheep and cows 
    [5.1e-7, 7.7e-7],   # transmissibility parameters for sheep and cows 
    1.                  # diameter of ring culling
)
u0_2 = u0_seirc76(  # Initial conditions for the model
    4000,               # number of farms
    1,                  # Y0 = initial number of farms with infections
    20,                 # size of grid 
    p_2;                # the parameters defined above
    seed = 762
) 
duration_2 = 400    # Duration

result76_2 = run_seirc76(u0_2, p_2, duration_2; seed = 7620)
df76_2 = dataframe_seirc76(result76_2)
            
plot_seirc76(df76_2)
video_seirc76(result76_2, df76_2; filename = "video76_2.mp4")

### Aggressive ring culling (diameter = 2 km) and 4 initially infectious

p_agg = Parameters76(  # Model parameters
    [1., 10.5],         # susceptibility parameters for sheep and cows 
    [5.1e-7, 7.7e-7],   # transmissibility parameters for sheep and cows 
    2.                  # diameter of ring culling
)
u0_agg = u0_seirc76(   # Initial conditions for the model
    4000,               # number of farms
    4,                  # Y0 = initial number of farms with infections
    20,                 # size of grid 
    p_2;                # the parameters defined above
    seed = 763
) 
duration_agg = 400     # Duration

result76_agg = run_seirc76(u0_agg, p_agg, duration_agg; seed = 7630)
df76_agg = dataframe_seirc76(result76_agg)
            
plot_seirc76(df76_agg)
video_seirc76(result76_agg, df76_agg; filename = "video76_agg.mp4")


## Programme 7.7

include("src/chapter7/p7.7.jl")
using .MID_77

N = 100                 # Number of individuals in the model 
connections = 4         # Average number of connections per individual
y0 = 1                  # Initial number of infectious individuals 
p = [               # Model parameters 
    .1,                 # gamma = recovery rate  
    1.                  # tau = transmission rates to contacts 
]
duration = 50       # Duration

### Random network

u0_rand = u0_sis77(N, connections,  y0, :random; seed = 77)

u_rand, times_rand = run_sis77(u0_rand, p, duration; seed = 770)
result77_rand = dataframe_sis77(u_rand, times_rand, N)

### Lattice network

u0_lat = u0_sis77(N, connections, y0, :lattice; seed = 77)

u_lat, times_lat = run_sis77(u0_lat, p, duration; seed = 770)
result77_lat = dataframe_sis77(u_lat, times_lat, N)

### Small world network

u0_sw = u0_sis77(N, connections, y0, :smallworld; seed = 77)

u_sw, times_sw = run_sis77(u0_sw, p, duration; seed = 770)
result77_sw = dataframe_sis77(u_sw, times_sw, N)

### Figure with all three sets of results 
titles = ["Random", "Lattice", "Small world"]

fig = Figure() 
gl = GridLayout(fig[1, 1])
axs = [ Axis(gl[i, 1]) for i ∈ 1:3 ]
for (i, res) ∈ enumerate([result77_rand, result77_lat, result77_sw])
    plot_sis77!(axs[i], res)
    axs[i].title = titles[i]
end 
leg = Legend(gl[:, 2], axs[1])
fig 
