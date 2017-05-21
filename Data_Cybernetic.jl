function Data_Cybernetic(tStart,tStop,tStep)

  #=============================#
  # Values for the entry kinetics of the virus
  Kd = 6.023*10.0^(13)  ;   # virions/ml
  ke = 72.0             ;   # 1/(day.virions)
  kb = 1.04065*10.0^(29);   # virions/(ml.day)
  kf = 10454.4          ;   # 1/day

  #=============================#
  # infection rate, beta
  beta = (kf*ke/(Kd+kf*ke/kb));

  #=============================#
  # infected cells death rate, delta_I
  delta_I = 2.0;    # cells/day

  #=============================#
  # IFN induced antiviral efficacy
  phi = 6.9*10.0^(-2) ; # 1/(IFN fold change. day)

  #============================#
  # Reversion rate from the refractory state
  rho = 1.0*10.0^(-2); # 1/day

  #============================#
  # Killing rate of infected cells by NK cells
  kappa = 1.6; # 1/(IFN fold change. day)

  #============================#
  # Viral production rate
  p = 7.7*10.0^(-5); # RNA copies/(ml NS).day.cell

  #============================#
  # Clearance rate of free virions
  c = 13.71; # 1/day

  #============================#
  # Production rate of IFN
  q = 6.1*10.0^(-10); # IFN fold change/(day.cell)

  #============================#
  # Decay rate of IFN
  d = 0.85 ; # 1/day


  #============================#
  # Rate of infected cell death increase
  omega = 1.2 ; # No units

  #============================#
  # Adaptive immune system kickstarted
  mu = 7.0 ;  # day

  #============================#
  # initial conditions
  x0 = [
    3.5*10.0^(11) ; # T initial target cells
    100.0         ; # V initial load of virus
    0.0           ; # R Refractory cells
    1.0           ; # F IFN starting fold change
    0.0           ; # I Infected cells
  ];

  #===========================#
  # Parameters from Data_dict
  Data_dict = Dict()
  Data_dict["InitialConditions"] = x0
  Data_dict["FusionRate"] = kf
  Data_dict["EndocytosisRate"] = ke
  Data_dict["ViralParticleDissociation"] = Kd
  Data_dict["BindingConstant"] = kb
  Data_dict["Clearance"] = c
  Data_dict["IFNProdRate"] = q
  Data_dict["IFNDeathRate"] = d
  Data_dict["ViralProdRate"] = p
  Data_dict["KillingRate"] = kappa
  Data_dict["RefractoryRate"] = rho
  Data_dict["AntiviralEff"] = phi
  Data_dict["InfectedDeath"] = delta_I
  Data_dict["InfectionRate"] = beta
  Data_dict["AdaptiveDelay"] = mu
  Data_dict["DeathFactor"] = omega

  return Data_dict
end
