
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set up 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Enter the path to this file here 
loc = "C:\\Users\\yourname\\Documents\\GitHub\\ModelingInfectiousDiseases.jl"
cd(loc)

using Pkg 
Pkg.activate(loc)
Pkg.instantiate()

include("src/ModelingInfectiousDiseases.jl")

using BenchmarkTools


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Chapter 2 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Programme 2.1

using .MID_21

beta = 520 / 365        # infectiousness parameter
gamma = 1 / 7           # recovery rate
S0 = 1 - 1e-6           # initial proportion susceptible 
I0 = 1e-6               # initial proportion resistant
duration = 70           # duration of model

sol21 = run_sir21(; beta, gamma, S0, I0, duration, saveat = .125) # more frequent saveat to give smooth plot
plot_sir21(sol21)


## Programme 2.2

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

using .MID_25

beta = 520 / 365        # infectiousness parameter 
gamma = 1 / 7           # recovery rate 
I0 = 1e-6               # initial proportion infectious
duration = 70           # duration of model

sol25 = run_sis25(; beta, gamma, I0, duration, saveat = .125) # more frequent saveat to give a smooth plot
plot_sis25(sol25)


## Programme 2.6 

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

using .MID_31

beta = [10 .1; .1 1]    # matrix of infectiousness parameters [β_hh, β_hl; β_lh β_ll]
gamma = 1               # recovery rate 
nh = .2                 # proportion of high-risk individuals 
Ih = 1e-5               # initial proportion of infectious high-risk individuals
Il = 1e-3               # initial proportion of infectious low-risk individuals
duration = 15           # duration of model

sol31 = run_sir31(; beta, gamma, nh, Ih, Il, duration, saveat = .025) # more frequent saveat to give a smooth plot
plot_sir31(sol31)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Chapter 6 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Programme 6.1 

using .MID_61

u0 = [
    1e5,                # X0
    500,                # Y0
    1e6 - (1e5 + 500)   # Z0
]

# Run with no noise
p_nonoise = Parameters61(
    1.,                 # beta 
    .1,                 # gamma 
    1 / (50 * 365),     # mu 
    1 / (50 * 365),     # nu 
    0.                  # xi
)

results61_nonoise = run_sir61(u0, p_nonoise, 5 * 365)
plot_sir61(results61_nonoise)

# Run with noise parameter
p = Parameters61(
    1.,                 # beta 
    .1,                 # gamma 
    1 / (50 * 365),     # mu 
    1 / (50 * 365),     # nu 
    10.                 # xi
)

results61 = run_sir61(u0, p, 5 * 365; seed = 61)
plot_sir61(results61)


## Programme 6.2 

using .MID_62

u0 = [
    1e5,                # X0
    500,                # Y0
    1e6 - (1e5 + 500)   # Z0
]

# Run with no noise
p_nonoise = Parameters62(
    1.,                 # beta 
    .1,                 # gamma 
    1 / (50 * 365),     # mu 
    1 / (50 * 365),     # nu 
    0.                  # xi
)

results62_nonoise = run_sir62(u0, p_nonoise, 5 * 365)
plot_sir62(results62_nonoise)

# Run with noise parameter
p = Parameters62(
    1.,                 # beta 
    .1,                 # gamma 
    1 / (50 * 365),     # mu 
    1 / (50 * 365),     # nu 
    1.                  # xi
)

results62 = run_sir62(u0, p, 5 * 365; seed = 62)
plot_sir62(results62)


## Programme 6.3 

using .MID_63

u0 = [
    30,                 # X0 
    70                  # Y0
]
p = [
    .03,                # beta 
    .01                 # gamma 
]

results63 = run_sis63(u0, p, 10 * 365; seed = 63)
plot_sis63(results63)
