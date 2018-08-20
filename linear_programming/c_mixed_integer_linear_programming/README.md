# c_mixed_integer_linear_programming

This project is a rational mixed integer linear programming solver written in C that can be used to build custom mixed integer linear programming optimization models for many applications.


## Inputs

Currently there are thirteen input files to the program in the form of tab delimited text files that the solver will use to find the optimal values for the decision variables.

### number_of_constraints.txt

This input file contains just one integer value which is the number of constraints which is chosen from the set of non-negative integers, i.e. 0 &le; n &lt; &infin;. This will be based on your individual problem.

### number_of_variables.txt

This input file contains just one integer value which is the number of decision variables which is chosen from the set of non-negative integers, i.e. 0 &le; n &lt; &infin;. This will be based on your individual problem.

### variable_special_requirements.txt

This input file contains integers with a new line for each of the decision variables indiciating any special requirements for each of the decision variables such as integer, binary, etc.

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

### initial_LP_objective_function_optimal_value.txt
This output file contains the initial linear programming (LP relaxation) objective function's optimal value as a float.

### initial_LP_basic_feasible_solution.txt
This output file contains the initial linear programming (LP relaxation) basic feasible solution of the decision variables as a float on a new line for each variable.

### initial_LP_tableau_matrix.txt
This output file contains the initial linear programming (LP relaxation) tableau matrix.

### final_LP_objective_function_optimal_value.txt
This output file contains the final linear programming (LP relaxation) objective function's optimal value as a float.

### final_LP_basic_feasible_solution.txt
This output file contains the final linear programming (LP relaxation) basic feasible solution of the decision variables as a float on a new line for each variable.

### final_LP_tableau_matrix.txt
This output file contains the final linear programming (LP relaxation) tableau matrix.

### final_MILP_optimal_objective_function_value.txt
This output file contains the final mixed integer linear programming objective function's optimal value as a float.

### final_MILP_optimal_variable_values.txt
This output file contains the final mixed integer linear programming optimal decision variable values as a float on a new line for each variable.
