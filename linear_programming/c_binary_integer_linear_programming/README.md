# c_binary_integer_linear_programming

This project is a rational binary integer linear programming solver written in C using the Balas Additive algorithm that can be used to build custom binary integer linear programming optimization models for many applications.


## Inputs

Currently there are thirteen input files to the program in the form of tab delimited text files that the solver will use to find the optimal values for the decision variables.

### number_of_constraints.txt

This input file contains just one integer value which is the number of constraints which is chosen from the set of non-negative integers, i.e. 0 &le; n &lt; &infin;. This will be based on your individual problem.

### number_of_variables.txt

This input file contains just one integer value which is the number of decision variables which is chosen from the set of non-negative integers, i.e. 0 &le; n &lt; &infin;. This will be based on your individual problem.

### objective_function_maximization.txt

This input file containts just one integer value 0 or 1 indicating if the optimization problem type is maximization or not(minimization).

### objective_function_initial_constant_numerator.txt & objective_function_initial_constant_denominator.txt

These input files, representing the numerator and denominator, respectively, each contain one integer value representing the objective function's initial constant.

### objective_function_coefficient_vector_numerator.txt & objective_function_coefficient_vector_denominator.txt

These input files, representing the numerator and denominator, respectively, each contain a new line for each decision variable's objective function coefficient.

### constraint_inequality_direction_vector.txt

This input file contains integers in the range [-1, 1] with a new line for each of the constraint inequality directions.

### constraint_constant_vector_numerator.txt & constraint_constant_vector_denominator.txt

These input files, representing the numerator and denominator, respectively, each contain a new line for each constraint's inequality constant.

### constraint_coefficient_matrix_numerator.txt & constraint_coefficient_matrix_denominator.txt

These input files, representing the numerator and denominator, respectively, each contain a matrix of number of constraints by number of variables that represent for each new line that constraint's decision variable coefficients.


## Outputs

The outputs save the optimal decision variable values and the associated optimal value.

### final_BILP_optimal_objective_function_value.txt
This output file contains the final binary integer linear programming objective function's optimal value as a float.

### final_BILP_optimal_variable_values.txt
This output file contains the final binary integer linear programming optimal decision variable values as a float on a new line for each variable.
