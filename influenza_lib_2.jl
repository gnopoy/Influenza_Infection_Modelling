include("Data_Cybernetic.jl")
include("SolveBalances.jl")

using PyCall
@pyimport numpy as np

#some global parameters
BIG = 1e10
SMALL = 1e-16

# Globally load the data so we don't have to load it every iteration
MEASURED_ARRAY_1 = readdlm("./data/Experimental_data_Pony_1_Viremia.dat")
MEASURED_ARRAY_2 = readdlm("./data/Experimental_data_Pony_1_IFN.dat")
# MEASURED_ARRAY_3 = readdlm("./data/Experimental_data_Pony_2_Viremia.dat")
# MEASURED_ARRAY_4 = readdlm("./data/Experimental_data_Pony_2_IFN.dat")
# MEASURED_ARRAY_5 = readdlm("./data/Experimental_data_Pony_3_Viremia.dat")
# MEASURED_ARRAY_6 = readdlm("./data/Experimental_data_Pony_3_IFN.dat")
# MEASURED_ARRAY_7 = readdlm("./data/Experimental_data_Pony_4_Viremia.dat")
# MEASURED_ARRAY_8 = readdlm("./data/Experimental_data_Pony_4_IFN.dat")
# MEASURED_ARRAY_9 = readdlm("./data/Experimental_data_Pony_5_Viremia.dat")
# MEASURED_ARRAY_10 = readdlm("./data/Experimental_data_Pony_5_IFN.dat")
# MEASURED_ARRAY_11 = readdlm("./data/Experimental_data_Pony_6_Viremia.dat")
# MEASURED_ARRAY_12 = readdlm("./data/Experimental_data_Pony_6_IFN.dat")


function local_refinement_step(parameter_array)

  SIGMA = 0.05

  # initalize -
  number_of_parameters = length(parameter_array)

  # Setup the bounds constraints
  LOWER_BOUND = SMALL
  UPPER_BOUND = [
  30.00 ; # 1  rho
  30.00 ; # 2  phi
  30.00 ; # 3  kappa
  1.00 ; # 4  p
  1.00 ; # 5  q
  30.00 ; # 6  c
  30.00 ; # 7  d
  30.00 ; # 8  mu
  6.00 ; # 9  delta_I
  1.00 ; # 10 omega
  ]
  # calculate the starting error
  parameter_array_best = parameter_array
  error_array = BIG*ones(2)
  erro_array[1] = sum(objective_function(parameter_array_best))

  # Main refinement loop
  iteration_counter = 1
  iteration_max = 1000
  while (iteration_counter<iteration_max)

    #take a step up -
    parameter_up = parameter_array_best.*(1+SIGMA*rand(number_of_parameters))
    parameter_up = parameter_bounds_function(parameter_up,LOWER_BOUND*ones(number_of_parameters),UPPER_BOUND)

    #take a step down -
    parameter_down = parameter_array_best.*(1-SIGMA*rand(number_of_parameters))
    parameter_down = parameter_bounds_function(parameter_down,LOWER_BOUND*ones(number_of_parameters),UPPER_BOUND)

    #Evaluate the objective_function
    error_array[2] = sum(objective_function(parameter_up))
    error_array[3] = sum(objective_function(parameter_down))

    # Calculate a correction factor
    a = error_array[2]+error_array[3]-2.0*error_array[1]
    parameter_corrected = parameter_array_best
    if (a>0.0)
      amda = -0.5*(error_array[3]-error_array[2])/a
      parameter_corrected = parameter_array_best+amda*rand(number_of_parameters)
      parameter_corrected = parameter_bounds_function(parameter_corrected,LOWER_BOUND*ones(number_of_parameters),UPPER_BOUND)
      error_array[4] = sum(objective_function(parameter_corrected))
    end

    #which step has the minimum error?
    min_index =indmin(error_array)
    if (min_index == 1)
      parameter_array_best = parameter_array_best
    elseif (min_index == 2)
      parameter_array_best = parameter_up
    elseif (min_index == 3)
      parameter_array_best = parameter_down
    elseif (min_index == 4)
      parameter_array_best = parameter_corrected
    end

    #update the local error
    error_array[1] = error_array[min_index]

    @show iteration_counter,error_array[min_index]

    #update local counter -
    iteration_counter =iteration_counter + 1
  end

  return parameter_array_best
end


# write the calculate_error_m_ functions for each Pony 1-6
function calculate_error_m_1(t,x)
  tStart = 0.0 ;
  tStop  = 14.0 ;
  tStep  = 0.1 ;

  obj_array = BIG*ones(2)

  # Need to interpolate the simulations onto  the experimental time scale -
  # for the viremia
  number_of_measurements  = 10
  time_experimental = linspace(tStart,tStop,number_of_measurements)
  AI = np.interp(time_experimental[:],t,x[:,4])
  BI = np.interp(time_experimental[:],t,x[:,5])


  # Interpolate the experimental data onto the same timescale -
  AMI = np.interp(time_experimental[:],MEASURED_ARRAY_1[:,1],MEASURED_ARRAY_1[:,2])
  BMI = np.interp(time_experimental[:],MEASURED_ARRAY_2[:,1],MEASURED_ARRAY_2[:,2])

  # Compute the error values -
  error_A = sum((AMI-AI).^2)
  error_B = sum((BMI-BI).^2)

  # Compute the objective array -
  obj_array[1] = error_A
  obj_array[2] = error_B

  return (sum(obj_array))
end
#
# function calculate_error_m_2(t,x)
#   tStart = 0.0 ;
#   tStop  = 14.0 ;
#   tStep  = 0.1 ;
#
#   obj_array = BIG*ones(2)
#
#   # Need to interpolate the simulations onto  the experimental time scale -
#   # for the viremia
#   number_of_measurements  = 10
#   time_experimental = linspace(tStart,tStop,number_of_measurements)
#   AI = np.interp(time_experimental[:],t,x[:,4])
#   BI = np.interp(time_experimental[:],t,x[:,5])
#
#
#   # Interpolate the experimental data onto the same timescale -
#   AMI = np.interp(time_experimental[:],MEASURED_ARRAY_3[:,1],MEASURED_ARRAY_3[:,2])
#   BMI = np.interp(time_experimental[:],MEASURED_ARRAY_4[:,1],MEASURED_ARRAY_4[:,2])
#
#   # Compute the error values -
#   error_A = sum((AMI-AI).^2)
#   error_B = sum((BMI-BI).^2)
#
#   # Compute the objective array -
#   obj_array[1] = error_A
#   obj_array[2] = error_B
#
#   return (sum(obj_array))
# end
#
# function calculate_error_m_3(t,x)
#   tStart = 0.0 ;
#   tStop  = 14.0 ;
#   tStep  = 0.1 ;
#
#   obj_array = BIG*ones(2)
#
#   # Need to interpolate the simulations onto  the experimental time scale -
#   # for the viremia
#   number_of_measurements  = 10
#   time_experimental = linspace(tStart,tStop,number_of_measurements)
#   AI = np.interp(time_experimental[:],t,x[:,4])
#   BI = np.interp(time_experimental[:],t,x[:,5])
#
#
#   # Interpolate the experimental data onto the same timescale -
#   AMI = np.interp(time_experimental[:],MEASURED_ARRAY_5[:,1],MEASURED_ARRAY_5[:,2])
#   BMI = np.interp(time_experimental[:],MEASURED_ARRAY_6[:,1],MEASURED_ARRAY_6[:,2])
#
#   # Compute the error values -
#   error_A = sum((AMI-AI).^2)
#   error_B = sum((BMI-BI).^2)
#
#   # Compute the objective array -
#   obj_array[1] = error_A
#   obj_array[2] = error_B
#
#   return (sum(obj_array))
# end
#
# function calculate_error_m_4(t,x)
#   tStart = 0.0 ;
#   tStop  = 14.0 ;
#   tStep  = 0.1 ;
#
#   obj_array = BIG*ones(2)
#
#   # Need to interpolate the simulations onto  the experimental time scale -
#   # for the viremia
#   number_of_measurements  = 10
#   time_experimental = linspace(tStart,tStop,number_of_measurements)
#   AI = np.interp(time_experimental[:],t,x[:,4])
#   BI = np.interp(time_experimental[:],t,x[:,5])
#
#
#   # Interpolate the experimental data onto the same timescale -
#   AMI = np.interp(time_experimental[:],MEASURED_ARRAY_7[:,1],MEASURED_ARRAY_7[:,2])
#   BMI = np.interp(time_experimental[:],MEASURED_ARRAY_8[:,1],MEASURED_ARRAY_8[:,2])
#
#   # Compute the error values -
#   error_A = sum((AMI-AI).^2)
#   error_B = sum((BMI-BI).^2)
#
#   # Compute the objective array -
#   obj_array[1] = error_A
#   obj_array[2] = error_B
#
#   return (sum(obj_array))
# end
#
# function calculate_error_m_5(t,x)
#   tStart = 0.0 ;
#   tStop  = 14.0 ;
#   tStep  = 0.1 ;
#
#   obj_array = BIG*ones(2)
#
#   # Need to interpolate the simulations onto  the experimental time scale -
#   # for the viremia
#   number_of_measurements  = 10
#   time_experimental = linspace(tStart,tStop,number_of_measurements)
#   AI = np.interp(time_experimental[:],t,x[:,4])
#   BI = np.interp(time_experimental[:],t,x[:,5])
#
#
#   # Interpolate the experimental data onto the same timescale -
#   AMI = np.interp(time_experimental[:],MEASURED_ARRAY_9[:,1],MEASURED_ARRAY_9[:,2])
#   BMI = np.interp(time_experimental[:],MEASURED_ARRAY_10[:,1],MEASURED_ARRAY_10[:,2])
#
#   # Compute the error values -
#   error_A = sum((AMI-AI).^2)
#   error_B = sum((BMI-BI).^2)
#
#   # Compute the objective array -
#   obj_array[1] = error_A
#   obj_array[2] = error_B
#
#   return (sum(obj_array))
# end
#
# function calculate_error_m_6(t,x)
#   tStart = 0.0 ;
#   tStop  = 14.0 ;
#   tStep  = 0.1 ;
#
#   obj_array = BIG*ones(2)
#
#   # Need to interpolate the simulations onto  the experimental time scale -
#   # for the viremia
#   number_of_measurements  = 10
#   time_experimental = linspace(tStart,tStop,number_of_measurements)
#   AI = np.interp(time_experimental[:],t,x[:,4])
#   BI = np.interp(time_experimental[:],t,x[:,5])
#
#
#   # Interpolate the experimental data onto the same timescale -
#   AMI = np.interp(time_experimental[:],MEASURED_ARRAY_11[:,1],MEASURED_ARRAY_11[:,2])
#   BMI = np.interp(time_experimental[:],MEASURED_ARRAY_12[:,1],MEASURED_ARRAY_12[:,2])
#
#   # Compute the error values -
#   error_A = sum((AMI-AI).^2)
#   error_B = sum((BMI-BI).^2)
#
#   # Compute the objective array -
#   obj_array[1] = error_A
#   obj_array[2] = error_B
#
#   return (sum(obj_array))
# end

function objective_function(parameter_array)
  # Calculate objective function array -
  obj_array = BIG*ones(1,1)

  # Solve the model with the update parameters -
  tStart =  0.0  ;
  tStop =   14.0 ;
  tStep =   0.1  ;
  (t,x) = SolveBalances(tStart,tStop,tStep,parameter_array)

  # Call the experiment functions
  obj_array[1] = calculate_error_m_1(t,x)
  # obj_array[2] = calculate_error_m_2(t,x)
  # obj_array[3] = calculate_error_m_3(t,x)
  # obj_array[4] = calculate_error_m_4(t,x)
  # obj_array[5] = calculate_error_m_5(t,x)
  # obj_array[6] = calculate_error_m_6(t,x)

  # return
  return obj_array
end

# Generates new parameter array, given current array -
function neighbor_function(parameter_array)

  SIGMA = 0.05
  number_of_parameters = length(parameter_array)

  # Calculate new parameters -
  new_parameter_array = parameter_array.*(1+SIGMA*randn(number_of_parameters))

  # Check the bound constraints -
  LOWER_BOUND = SMALL
  UPPER_BOUND = [
  30.00 ; # 1  rho
  30.00 ; # 2  phi
  30.00 ; # 3  kappa
  1.00 ; # 4  p
  1.00 ; # 5  q
  30.00 ; # 6  c
  30.00 ; # 7  d
  30.00 ; # 8  mu
  6.00 ; # 9  delta_I
  1.00 ; # 10 omega
  ]

  # return the corrected parameter arrays -
  return parameter_bounds_function(new_parameter_array, LOWER_BOUND*ones(number_of_parameters), UPPER_BOUND)
end

function acceptance_probability_function(rank_array, temperature)
  return (exp(-rank_array[end]/temperature))
end

function cooling_function(temperature)
  # define my new temperature -
  alpha = 0.9
  return alpha*temperature
end


# Helper functions

function parameter_bounds_function(parameter_array, lower_bound_array, upper_bound_array)
  # reflecttion_factor -
  epsilon = 0.01

  # iterate through and fix the parameters -
  new_parameter_array = copy(parameter_array)
  for (index,value) in enumerate(parameter_array)

    lower_bound = lower_bound_array[index]
    upper_bound = upper_bound_array[index]

    if (value<lower_bound)
      new_parameter_array[index] = lower_bound
    elseif (value>upper_bound)
      new_parameter_array[index] = upper_bound
    end
  end

  return new_parameter_array
end
