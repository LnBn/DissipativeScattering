##############################################################
## Section 0: The required packages
##############################################################

# saving, loading, organizing data 
using DataFrames
using CSV

# integrating ODEs
using DifferentialEquations
using DynamicalSystems

# curve fitting 
using LsqFit

# some basic functions 
using Statistics
using LinearAlgebra
using StaticArrays
using Base.Threads

# plotting 
using Plots

#############################################################
## Section 1: Dynamical System setup
#############################################################

## System constants, parameters, functions:

# potential energy - accepts individual configuration space coordinates or a position vector q = [x,y]
V(x,y) = 1//2 * (x^2 + y^2 + 2x^2*y - 2//3 * y^3);
V(q) = 1//2 * (q[1]^2 + q[2]^2 + 2q[1]^2 * q[2]- 2//3 * q[2]^3);

# positions of saddle points in the potential 
saddles = [[0.0,1.0],[-√3/2,-1/2],[√3/2,-1/2]] #saddle points


#kinetic energy - accepts momentum vector p = [̇x,̇y]
T(p) = 1//2 * norm(p)^2;

#total energy/hamiltonian
H(x,y,dx,dy) = V(x,y) + 1//2 * (dx^2 + dy^2);
H(p,q, params) = T(p) + V(q);

const r_escape = 1.03 #escape radius from potential 
const r_attractor = 1e-3 #collision radius for attractor (dissipative case only)

# maximum integration timespan (this is more than enough)
tspan = (0.0,5000.0)

########
######## Hamiltonian system is setup as a SecondOrderODEProblem() from DifferentialEquations.jl
######## Dissipative system is setup as a ContinuousDynamicalSystem from DynamicalSystems.jl
########

## Equations of motion (Hamiltonian): 

function HH_acceleration!(dv,v,u,p,t)
    x,y  = u
    dx,dy = dv
    dv[1] = -x - 2x*y
    dv[2] = y^2 - y -x^2
end

## Equations of motion (Dissipative - α here is γ in the paper):

@inline @inbounds function loop(u, p, t)
    α = p[1]; λ = p[2]
    dx = u[3]
    dy = u[4]
    dẋ = -u[1] - 2*u[1]*u[2] - α*u[3]
    dẏ = -u[1]^2 - u[2] + u[2]^2 - α*u[4]
    dE = -α*(u[3]^2 + u[4]^2)
    return SVector{5}(dx,dy,dẋ,dẏ,dE)
end

# negated equations of motion (t → -t). I used these as an independent way of finding the unstable manifold. Otherwise not used.


function HH_acceleration_neg!(dv,v,u,p,t) # t → -t
    x,y  = u
    dx,dy = dv
    dv[1] = x + 2x*y
    dv[2] = -y^2 + y + x^2
end



@inline @inbounds function loop_neg(u, p, t)
    α = p[1]; λ = p[2]
    dx = -u[3]
    dy = -u[4]
    dẋ = u[1] + 2*u[1]*u[2] + α*u[3]
    dẏ = u[1]^2 + u[2] - u[2]^2 + α*u[4]
    dE = α*(u[3]^2 + u[4]^2)
    return SVector{5}(dx,dy,dẋ,dẏ,dE)
end



########################################
## Section 1.2: Initial Conditions
########################################


## Supply an initial position (x,y) and energy E. 
# These functions will throw an error if V(x,y) > E 
# I also include an automated way to generate a large number of initial conditions (functions called "goodICs_...")
# but it is slow (generates a random point, throws it out if outside potential)
# DynamicalSystems.jl has a better version (Systems.henonheiles_ics()) but it assumes a different 
# Poincare section


## Hamiltonian:

function IC_H(x,y,E)
    q = x,y
    q0 = [x,y]
    dx = -y/norm(q) * sqrt(2(E-V(q)))
    dy = x/norm(q) * sqrt(2(E-V(q)))
    p0 = [dx,dy]
    return p0,q0
end


## Dissipative:

# function for cartesian coordinates:

function IC_D(x,y,E);
    q = x,y
    dx = -y/norm(q) * sqrt(2(E-V(q)))
    dy = x/norm(q) * sqrt(2(E-V(q)))
    return @SVector [x,y,dx,dy,H(x,y,dx,dy)]
end

#function for polar coordinates:

x(r,θ) = r*cos(θ); y(r,θ) = r*sin(θ);

function IC_D_disk(r,θ,E);
    xs = x(r,θ); ys = y(r,θ)
    q = xs,ys
    dx = -ys/norm(q) * sqrt(2(E-V(q)))
    dy = xs/norm(q) * sqrt(2(E-V(q)))
    return @SVector [xs,ys,dx,dy,H(xs,ys,dx,dy)]
end



########################################
## Section 1.3: Event Handling functions
########################################

# The events to watch for during integration are: escape (hamiltonian case), the escape channels closing (dissipative case),
# and/or the trajectory coming too close to attractor (dissipative case).


## Hamiltonian:

function escape_condition(u,t,integrator);
    norm(u[3:4]) > r_escape
end

function escape_affect!(integrator);
    terminate!(integrator)
end

cb_H = DiscreteCallback(escape_condition,escape_affect!);

## Dissipative:

function closure_condition(u,t,integrator);
    u[5] < 1/6
end

function closure_affect!(integrator);
    terminate!(integrator)
end

function attractor_condition(u,t,integrator)
    norm(u[1:2]) < r_attractor
end

function escape_diss_condition(u,t,integrator)
    (u[5] < 1/6) && norm(u[1:2]) > r_escape
end

function escape_diss2(u,t,integrator)
    norm(u[1:2]) > r_escape
end


function psos_condition_time(u,t,integrator)
    (t_crossing_i < t < t_crossing_f) && (abs(u[1]*u[3]+u[2]*u[4]) < 1e-3) && ((-u[3]/u[2] > 0) && (u[4]/u[1]) > 0) 
end

function psos_condition_energy(u,t,integrator)
    (Emin < u[5] < Emax) && (abs(u[1]*u[3]+u[2]*u[4]) < 1e-3) && ((-u[3]/u[2] > 0) && (u[4]/u[1]) > 0) 
end




###############################
## Callbacks and Solver Options
###############################
# Not all of these are used - I was experimenting with different choices 

tspan_comp = (0.0,80.0) #integration window for cmeasure computations (here is approx 0 to τ₂ for E = 0.35)   
t_crossing_i = 17.0 #min time to start looking for intersections
t_crossing_f = 30.0 #max time to look for intersections
Emin = 0.295 #min energy to look for intersections
Emax = 0.305 #max energy to look for intersections

cb_D = DiscreteCallback(closure_condition,closure_affect!);
cb_attractor = DiscreteCallback(attractor_condition,closure_affect!)
cb_escape_diss = DiscreteCallback(escape_diss_condition,closure_affect!)
cb_escape_diss2 = DiscreteCallback(escape_diss2,escape_affect!)
# cb_D2 = CallbackSet(cb_attractor,cb_escape_diss)
cb_D2 = CallbackSet(cb_attractor,cb_escape_diss2)
cb_D_plt = CallbackSet(cb_attractor,cb_escape_diss2)
cb_D3 = CallbackSet(cb_escape_diss2,cb_D)

cb_psos_time = DiscreteCallback(psos_condition_time,escape_affect!)
cb_psos_energy = DiscreteCallback(psos_condition_energy,escape_affect!)
cb_D4 = CallbackSet(cb_escape_diss2,cb_psos,cb_D) #for comparing H to diss case 
diffeq_D8 = (alg = Vern9(),saveat=5.0,callback=cb_D4,abstol=1e-10,reltol=1e-10,maxiters=1e10);




diffeq_D = (alg = Vern9(),saveat=0.1,callback=cb_D,abstol=1e-12,reltol=1e-12,maxiters=1e10);



## attractor stop condition
diffeq_D2 = (alg = Vern9(),saveat=0.1,callback=cb_D2,abstol=1e-12,reltol=1e-12,maxiters=1e10);
# diffeq_D2_eps = (alg = Vern9(),save_everystep=false,callback=cb_D3,abstol=1e-14,reltol=1e-14,maxiters=1e10);
diffeq_D2_eps = (alg = Vern9(),save_everystep=false,callback=cb_D3,abstol=1e-10,reltol=1e-10,maxiters=1e10);
diffeq_D2_eps_cl = (alg = Vern9(),save_everystep=false,callback=cb_D,abstol=1e-14,reltol=1e-14,maxiters=1e10); #save endpoints, closure condition
diffeq_D2_plt = (alg = Vern9(),saveat=0.05,callback=cb_D_plt,abstol=1e-14,reltol=1e-14,maxiters=1e10);

diffeq_D3 = (alg = Vern9(),saveat=0.1,callback=cb_D3,abstol=1e-14,reltol=1e-14,maxiters=1e10);


#############################################################
## Section 2: Ensembling
#############################################################

## In this section, the functions for computing ensembles of trajectories. 
## many functions have similar names. The difference is usually in the 
## solver options/what information is being saved.  

l = Threads.ReentrantLock() #threadlock to prevent race conditions when saving data multihtreaded


function create_problem(Energy::Real)
    r = sqrt(0.25*rand()); θ = 2π*rand()
    du0, u0 = IC_H(r*cos(θ),r*sin(θ),Energy)
    prob = SecondOrderODEProblem(
        HH_acceleration!,
        du0,
        u0,
        tspan,
        callback = cb_H,
        )
    return prob
end

function create_problem_diss(Energy::Real,γ::Real)
    u0 = IC_D_disk(sqrt(0.25*rand()),2π*rand(),Energy)
    ds = ContinuousDynamicalSystem(loop,u0,[γ,1.0])
    prob = ODEProblem(ds,tspan)
    return prob
end

function create_problem_set(Energy::Real,n::Integer)
    probset = [create_problem(Energy) for i ∈ 1:n]
    return probset
end

function create_problem_set_diss(Energy::Real,γ::Real,n::Integer)
    probset = [create_problem_diss(Energy,γ) for i ∈ 1:n]
    return probset
end




function EnsembleSurvivalTime_D(npts::Integer,Energy::Real,γ::Float64)
    
    probset = create_problem_set_diss(Energy,γ,npts)

    function prob_funcD(prob,i,repeat)
        prob = probset[i]
    end
        
    ## save full solution or just stopping time?
    # output_funcD(sol,i) = (sol,false)
    output_funcD(sol,i) = (sol.t[end],false)

    prob = probset[1]

    EnsProb = EnsembleProblem(prob,prob_func = prob_funcD,output_func = output_funcD)
    sim = solve(EnsProb,Vern9(),EnsembleThreads(),trajectories=npts;diffeq_D2_eps...) #only save end points

    return sim

end


function EnsembleSurvivalTime_H(npts::Integer,Energy::Real)

        probset = create_problem_set(Energy,npts) #create array of ode problems

        function prob_funcH(prob,i,repeat)
            prob = probset[i]
        end
        
        output_funcH(sol,i) = (sol,false)

        du0b, u0b = IC_H(0.2,0.2,Energy)

        prob = SecondOrderODEProblem(
                HH_acceleration!,
                du0b,
                u0b,
                tspan,
                callback = cb_H,
                )



        EnsProb = EnsembleProblem(prob,prob_func = prob_funcH,output_func = output_funcH)

        # sim = solve(EnsProb,KahanLi8(),EnsembleThreads(),dt=0.01, saveat=0.1,maxiters=1e10, trajectories = npts) #save every 0.1s
        sim = solve(EnsProb,KahanLi8(),EnsembleThreads(),dt=0.01, save_everystep=false,maxiters=1e10, trajectories = npts) #save first and last point

        return sim
end


function EnsembleSurvivalTime_H_timeonly(npts::Integer,Energy::Real)

    probset = create_problem_set(Energy,npts) #create array of ode problems

    function prob_funcH(prob,i,repeat)
        prob = probset[i]
    end
    
    output_funcH(sol,i) = (sol.t[end],false)

    du0b, u0b = IC_H(0.2,0.2,Energy)

    prob = SecondOrderODEProblem(
            HH_acceleration!,
            du0b,
            u0b,
            tspan,
            callback = cb_H,
            )



    EnsProb = EnsembleProblem(prob,prob_func = prob_funcH,output_func = output_funcH)

    # sim = solve(EnsProb,KahanLi8(),EnsembleThreads(),dt=0.01, saveat=0.1,maxiters=1e10, trajectories = npts) #save every 0.1s
    sim = solve(EnsProb,KahanLi8(),EnsembleThreads(),dt=0.01, save_everystep=false,maxiters=1e10, trajectories = npts) #save first and last point

    return sim
end


function EnsembleSurvivalTime_POS(npts::Integer,Energy::Real)

    probset = create_problem_set(Energy,npts) #create array of ode problems
    
    function prob_func(prob,i,repeat)
        prob = probset[i]
    end

    output_func(sol,i) = (sol,false)

    du0b, u0b = IC_H(0.2,0.2,Energy)

    prob = SecondOrderODEProblem(
            HH_acceleration!,
            du0b,
            u0b,
            tspan,
            callback = cb_H,
            )



    EnsProb = EnsembleProblem(prob,prob_func = prob_func,output_func = output_func)

    # sim = solve(EnsProb,KahanLi8(),EnsembleThreads(),dt=0.01, saveat=0.1,maxiters=1e10, trajectories = npts) #save every 0.1s
    sim = solve(EnsProb,KahanLi8(),EnsembleThreads(),dt=0.01, save_everystep=false,maxiters=1e10, trajectories = npts) #save first and last point

    return sim

end




function goodICs(E,n::Int64)
    ICs = []
    count = 0
    while length(ICs) < n
        @threads for i in 1:length(nthreads())
            count += 1
            x = -1 + 2*rand()
            y = -1 + 2*rand()
            # sqrt(x^2 + y^2) ≥ r_escape && continue
            V(x,y) ≥ E && continue
            q = [x,y]
            dx = -y/norm(q) * sqrt(2(E-V(q)))
            dy = x/norm(q) * sqrt(2(E-V(q)))
            ic = @SVector [x,y,dx,dy,H(x,y,dx,dy)]
            lock(l)
                push!(ICs,ic)
            unlock(l)
        end
    end
    println("$count points were required.")
    return ICs
end



function goodICs_H(E,n::Int64)
    ICs = []
    while length(ICs) < n
        for i in 1:nthreads()
            x = -1 + 2*rand()
            y = -1 + 2*rand()
            # sqrt(x^2 + y^2) ≥ r_escape && continue
            V(x,y) ≥ E && continue
            q = [x,y]
            dx = -y/norm(q) * sqrt(2(E-V(q)))
            dy = x/norm(q) * sqrt(2(E-V(q)))
            p = [dx,dy]
            ic = [p,q]
            lock(l)
            push!(ICs,ic)
            unlock(l)
        end
    end
    return ICs
end


function create_problem_set_diss2(E::Real,γ::Real,n::Integer)
    probset = []
    ICs = goodICs(E,n)
    for u0 in ICs
        ds = ContinuousDynamicalSystem(loop,u0,[γ,1.0])
        prob = ODEProblem(ds,tspan)
        push!(probset,prob)
    end
    return probset
end

function create_problem_set_diss2_tspan(E::Real,γ::Real,n::Integer,tspan)
    probset = []
    ICs = goodICs(E,n)
    for u0 in ICs
        ds = ContinuousDynamicalSystem(loop,u0,[γ,1.0])
        prob = ODEProblem(ds,tspan)
        push!(probset,prob)
    end
    return probset
end


function create_problem_set_H2(E::Real,n::Integer) #create problems over whole potential, not just in ball 
    probset = []
    ICs = goodICs_H(E,n)
    for ic in ICs
        prob = SecondOrderODEProblem(
        HH_acceleration!,
        ic[1],
        ic[2],
        tspan,
        callback = cb_H,
        )
        push!(probset,prob)
    end
    return probset
end


function EnsembleSurvivalTime_D_goodICs(npts::Integer,Energy::Real,γ::Float64)
    
    probset = create_problem_set_diss2(Energy,γ,npts)

    function prob_funcD(prob,i,repeat)
        prob = probset[i]
    end
        
    ## save full solution or just stopping time?
    # output_funcD(sol,i) = (sol,false)
    # output_funcD(sol,i) = (sol.t[end],false)
    output_funcD(sol,i) = ([sol[1],sol[end],sol.t[end]],false)

    prob = probset[1]

    EnsProb = EnsembleProblem(prob,prob_func = prob_funcD,output_func = output_funcD)
    sim = solve(EnsProb,Vern9(),EnsembleThreads(),trajectories=npts;diffeq_D2_eps...) #only save end points, closure condition

    return sim

end

function EnsembleSurvivalTime_D_goodICs_save01(npts::Integer,Energy::Real,γ::Float64)
    
    probset = create_problem_set_diss2(Energy,γ,npts)

    function prob_funcD(prob,i,repeat)
        prob = probset[i]
    end
        
    ## save full solution or just stopping time?
    output_funcD(sol,i) = (sol,false)
    # output_funcD(sol,i) = (sol.t[end],false)
    # output_funcD(sol,i) = ([sol[1],sol[end],sol.t[end]],false)

    prob = probset[1]

    EnsProb = EnsembleProblem(prob,prob_func = prob_funcD,output_func = output_funcD)
    sim = solve(EnsProb,Vern9(),EnsembleThreads(),trajectories=npts;diffeq_D2...)

    return sim

end


function EnsembleSurvivalTime_D_START_AND_END_POS(npts::Integer,Energy::Real,γ::Float64)
    
    probset = create_problem_set_diss2(Energy,γ,npts)

    function prob_funcD(prob,i,repeat)
        prob = probset[i]
    end
        
    ## save full solution or just stopping time?
    # output_funcD(sol,i) = (sol,false)
    output_funcD(sol,i) = ([sol[1],sol[end],sol.t[end]],false)

    prob = probset[1]

    EnsProb = EnsembleProblem(prob,prob_func = prob_funcD,output_func = output_funcD)
    sim = solve(EnsProb,Vern9(),EnsembleThreads(),trajectories=npts;diffeq_D2_eps_cl...) #only save end points, closure condition

    return sim

end

function EnsembleSurvivalTime_D_START_AND_END_POS_ZOOM(npts::Integer,Energy::Real,γ::Float64)
    
    probset = create_problem_set_diss2_zoom(Energy,γ,npts)

    function prob_funcD(prob,i,repeat)
        prob = probset[i]
    end
        
    ## save full solution or just stopping time?
    # output_funcD(sol,i) = (sol,false)
    output_funcD(sol,i) = ([sol[1],sol[end],sol.t[end]],false)

    prob = probset[1]

    EnsProb = EnsembleProblem(prob,prob_func = prob_funcD,output_func = output_funcD)
    sim = solve(EnsProb,Vern9(),EnsembleThreads(),trajectories=npts;diffeq_D2_eps_cl...) #only save end points, closure condition

    return sim

end



function EnsembleSurvivalTime_H2(npts::Integer,Energy::Real)

    probset = create_problem_set_H2(Energy,npts) #create array of ode problems

    function prob_funcH(prob,i,repeat)
        prob = probset[i]
    end
    
    # output_funcH(sol,i) = ([sol[1],sol[end],sol.t[end]],false)
    output_funcH(sol,i) = (sol.t[end],false)

    du0b, u0b = IC_H(0.2,0.2,Energy)

    prob = SecondOrderODEProblem(
            HH_acceleration!,
            du0b,
            u0b,
            tspan,
            callback = cb_H,
            )



    EnsProb = EnsembleProblem(prob,prob_func = prob_funcH,output_func = output_funcH)

    # sim = solve(EnsProb,KahanLi8(),EnsembleThreads(),dt=0.01, saveat=0.1,maxiters=1e10, trajectories = npts) #save every 0.1s
    sim = solve(EnsProb,KahanLi8(),EnsembleThreads(),dt=0.01, save_everystep=false,maxiters=1e10, trajectories = npts) #save first and last point

    return sim
end



function EnsembleSurvivalTime_D4(npts::Integer,Energy::Real,γ::Float64)
    
    probset = create_problem_set_diss2_tspan(Energy,γ,npts,tspan_comp)

    function prob_funcD(prob,i,repeat)
        prob = probset[i]
    end
        
    ## save full solution or just stopping time?
    # output_funcD(sol,i) = (sol,false)
    # output_funcD(sol,i) = (sol.t[end],false)
    output_funcD(sol,i) = (sol,false)

    prob = probset[1]

    EnsProb = EnsembleProblem(prob,prob_func = prob_funcD,output_func = output_funcD)

    sim = solve(EnsProb,Vern9(),EnsembleThreads(),trajectories=npts;diffeq_D8...) #includes psos crossing condition t >15
    return sim

end


#####################################################################
#### Survival probability
#####################################################################

function SurvivalProbabilityCurve(df::DataFrame)
    # df = df[df.t .< quantile(df.t,0.99),:]
    if maximum(df.t) < 50.0
        tvec = 0.0:0.01:maximum(df.t)
    else
        tvec = 0.0:0.1:maximum(df.t)
    end 
    pvec = Real[]

    for t in tvec
        dft = df[df.t .< t,:]
        push!(pvec,1.0 .- length(dft.t)/length(df.t))
    end
    dfs = DataFrame(t = tvec,p = pvec)
    return dfs
end

## the above function can take a long time if tspan is long 
# here is a version with a maximum time 
function SurvivalProbabilityCurve_shorttime(df::DataFrame,tmax::Real)
    # df = df[df.t .< quantile(df.t,0.99),:]
    tvec = 0.0:0.1:tmax
    pvec = Real[]

    for t in tvec
            dft = df[df.t .< t,:]
            push!(pvec,1.0 .- length(dft.t)/length(df.t))
    end
    dfs = DataFrame(t = tvec,p = pvec)
    return dfs
end


######################################################################
#### Conditionally Invariant Measure
######################################################################



