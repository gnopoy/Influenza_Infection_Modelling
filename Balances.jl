using Sundials


function Balances(t,x,dxdt_vector,DF)

  idx = find(x<10.0^(-9))
  x[idx] = 10.0^(-9)
  #dxdt_vector = similar(x)
  # Alias the species vector
  T = x[1];
  V = x[2];
  R = x[3];
  F = x[4];
  I = x[5];

  # Get Parameters from DF

  c       = DF["Clearance"]
  q       = DF["IFNProdRate"]
  d       = DF["IFNDeathRate"]
  p       = DF["ViralProdRate"]
  kappa   = DF["KillingRate"]
  rho     = DF["RefractoryRate"]
  phi     = DF["AntiviralEff"]
  delta_I = DF["InfectedDeath"]
  beta    = DF["InfectionRate"]
  mu      = DF["AdaptiveDelay"]
  omega   = DF["DeathFactor"]

  # Define the death rate of infected cells due to AdaptiveDelay
  delta_A = delta_I*exp(omega*(t-mu))

  # ODEs for the model
  if (t<mu)
    dxdt_vector[1]= rho*x[3]-beta*x[4]*x[1]-phi*x[5]*x[1]
    dxdt_vector[2]= beta*x[4]*x[1]-delta_I*x[2]-kappa*x[2]*x[5]
    dxdt_vector[3]= phi*x[5]*x[1]-rho*x[3]
    dxdt_vector[4]= p*x[2]-c*x[4]
    dxdt_vector[5]= q*x[2]-d*x[5]

  else
    dxdt_vector[1]= rho*x[3]-beta*x[4]*x[1]-phi*x[5]*x[1]
    dxdt_vector[2]= beta*x[4]*x[1]-delta_A*x[2]-kappa*x[2]*x[5]
    dxdt_vector[3]= phi*x[5]*x[1]-rho*x[3]
    dxdt_vector[4]= p*x[2]-c*x[4]
    dxdt_vector[5]= q*x[2]-d*x[5]
  end
  #dxdt_vector
end
