include("influenza_lib_2.jl")
using POETs
using Sundials

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


(EC,PC,RA) = estimate_ensemble(objective_function,neighbor_function,acceptance_probability_function,cooling_function,initial_parameter_array;rank_cutoff=4,maximum_number_of_iterations=40,show_trace=true)
