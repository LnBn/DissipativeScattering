# DissipativeScattering
This repository contains the Julia code used to generate the results in the paper "Probabilistic description of dissipative chaotic scattering".

## Packages
The simulations were performed mainly using the  [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) and [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystems.jl/dev/) libraries. Other third-party libraries used were:

1. [CSV.jl](https://csv.juliadata.org/stable/) and [DataFrames.jl](https://dataframes.juliadata.org/stable/) (for saving/loading/manipulating data) 
2. [LsqFit.jl](https://github.com/JuliaNLSolvers/LsqFit.jl) (for curve-fitting, where necessary)
3. [Plots.jl](https://docs.juliaplots.org/stable/) (for visualisations).

**Familiarity with the syntax and operation of these tools is necessary.**

## Setup
In this paper we considered the Hénon-Heiles system with an added friction $\gamma$ : 

$$
\begin{equation}
\begin{aligned}
    &\dot{x} = p_x,\\
    &\dot{y} = p_y,\\
    &\dot{p_x} = -\omega^2x - 2\lambda xy - \gamma p_x,\\
    &\dot{p_y} = -\omega^2y - \lambda(x^2 - y^2) - \gamma p_y,
\end{aligned}
\end{equation}
$$

We were interested in statistical properties of ensembles of trajectories. All of the source codes are included in `dissipative_scattering_functions.jl`, which can loaded with `include("/path/to/dissipative_scattering_functions.jl")`. 


## Example Usage
### Survival Probability

To simulate an ensemble of $npts$ initial conditions chosen uniformly from the configuration space $(x,y)$ with initial energy $E$ and dissipation strength $\gamma$,
you can call `EnsembleSurvivalTime_D(npts,E,γ)`. This returns an `EnsembleSolution` object. The information returned by the integrator is user specified (see the documentation for `DifferentialEquations.jl`),
but by default `EnsembleSurvivalTime_D` supplies only the escape/settling time of each trajectory (for memory reasons). When the escape time is presented as a `DataFrame` under a column `t`, the function `SurvivalProbabilityCurve`
will generate another `DataFrame` with two columns `t` and `p`, information which represents the survival probability. 

To turn this information into a plot of the survival probability, one can do the following (here we choose (`npts`,E,γ) = (10^5,0.3,0.01)) :

```julia
sim = EnsembleSurvivalTime_D(10^5,0.3,0.01)

df = DataFrame(t = sim.u)
dfs = SurvivalProbabilityCurve(df)
plot(dfs.t,dfs.p,yscale=:log10,lab=false,xlab="t",ylab="P(t)",title="Survival Probability")

```
which produces 

![Survival Probability](Images/survprob.png)

