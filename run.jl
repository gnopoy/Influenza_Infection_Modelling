include("SolveBalances.jl")

using ODE
using PyPlot

#Define the timescale
tStart = 0.0 ;
tStop = 14.0 ;
tStep = 0.1  ;

#Load Parameters values
Param = Float64[]

#SolverResults
(t,x) = SolveBalances(tStart,tStop,tStep,Param)

#Define species from SolverResults
T = x[:,1]
V = x[:,2]
R = x[:,3]
F = x[:,4]
I = x[:,5]

using PyPlot
figure(1)
plot(t,T,color="k",label="Target Cells",linewidth=2)
plot(t,I,color="r",label="Infected Cells",linewidth=2)
plot(t,R,color="b",label="Refractory Cells",linewidth=2)
legend(fontsize=18)
xlabel("Time in days",fontsize=20)
ylabel("Cells",fontsize=20)
