include("Data_Cybernetic.jl")
include("Balances.jl")

using ODE
using Sundials

function SolveBalances(tStart,tStop,tStep,Param)

  # Define the timescale
  t = collect(tStart:tStep:tStop)

  # Load the initial conditions
  Data_dict = Data_Cybernetic(tStart,tStop,tStep)
  IC = Data_dict["InitialConditions"]

  if isempty(Param)
    DF = Data_dict
  else
    rho     = Param[1]
    phi     = Param[2]
    kappa   = Param[3]
    p       = Param[4]
    q       = Param[5]
    c       = Param[6]
    d       = Param[7]
    mu      = Param[8]
    delta_I = Param[9]
    omega   = Param[10]

    # Rewrite Data Dictionary
    Data_dict["Clearance"]       = c
    Data_dict["IFNProdRate"]     = q
    Data_dict["IFNDeathRate"]    = d
    Data_dict["ViralProdRate"]   = p
    Data_dict["KillingRate"]     = kappa
    Data_dict["RefractoryRate"]  = rho
    Data_dict["AntiviralEff"]    = phi
    Data_dict["InfectedDeath"]   = delta_I
    Data_dict["InfectionRate"]   = beta
    Data_dict["AdaptiveDelay"]   = mu
    Data_dict["DeathFactor"]     = omega
    DF                           = Data_dict
  end

  # Run Solver
  f(t,x,dxdt_vector) = Balances(t,x,dxdt_vector,DF)
  x = Sundials.cvode(f,IC,t,reltol=10.0^(-3),abstol=10.0^(-6))
  #t,x = ODE.ode23s(Balances,IC,t)

  return (t,x)
end
