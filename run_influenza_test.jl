# Estimate the model parameters for the Within-Host Dynamics of Influenza Infection with Entry kinetics
using POETs
#include("influenza_lib.jl")
include("influenza_lib_2.jl")

function test_influenza_model()

  number_of_subdivisions = 10
  number_of_parameters   = 10
  number_of_objectives   = 6
  number_of_objectives = 1

  # generate parameter guess -
  initial_parameter_array     = zeros(10)
  initial_parameter_array[1]  = 1.0*10.0^(-2)
  initial_parameter_array[2]  = 6.9*10.0^(-2)
  initial_parameter_array[3]  = 1.6
  initial_parameter_array[4]  = 7.7*10.0^(-5)
  initial_parameter_array[5]  = 6.1*10.0^(-10)
  initial_parameter_array[6]  = 13.71
  initial_parameter_array[7]  = 0.85
  initial_parameter_array[8]  = 7.0
  initial_parameter_array[9]  = 2.0
  initial_parameter_array[10] = 0.95

  ec_array = zeros(number_of_objectives)
  pc_array = zeros(number_of_parameters)

  for index in collect(1:number_of_subdivisions)
    # Run JuPOETs -
    (EC,PC,RA) = estimate_ensemble(objective_function,neighbor_function,acceptance_probability_function,cooling_function,initial_parameter_array;rank_cutoff=4,maximum_number_of_iterations=40,show_trace=true)

    # Package -
    ec_array = [ec_array EC]
    pc_array = [pc_array PC]

    # Take the *best* value from the current ensemble, refine it and go around again -
    total_error = sum(ec_array[:,2:end],1)

    # Which col is the min error?
    min_index = indmin(total_error)
    @show index,total_error[min_index]

    # Refine the best solution -
    initial_parameter_array = pc_array[:,min_index].*(1+0.15*randn(number_of_parameters))
    initial_parameter_array = local_refinement_step(initial_parameter_array)
  end

  return (ec_array,pc_array)
end
