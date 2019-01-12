#include <stdio.h> /* printf, scanf, puts */
#include <stdlib.h> /* realloc, free, exit, NULL */
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>

unsigned int max_branch_and_bound_recursion_depth = 6;

/*********************************************************************************************************/
/********************************************** STRUCTURES ***********************************************/
/*********************************************************************************************************/

struct BranchAndBoundState
{
	unsigned int tableau_current_size[2];
	unsigned int number_of_constraints;
	unsigned int number_of_slack_surplus_variables;
	unsigned int number_of_artificial_variables;

	unsigned int *basic_variables;
	long long **basic_feasible_solution;
	long long ***tableau_matrix;
};

/*********************************************************************************/
/********************************** PROTOTYPES ***********************************/
/*********************************************************************************/

/*********************************************************************************/
/********************************* READ INPUTS ***********************************/
/*********************************************************************************/

/* This function reads in if the objective function is going to be a maximization or minimization problem */
void ReadObjectiveFunctionMaximizationOrMinimization(int *maximization_problem);

/* This function reads the number of constraints and variables */
void ReadNumberOfConstraintsAndVariables(unsigned int *initial_number_of_constraints, unsigned int *number_of_constraints, unsigned int *number_of_variables);

/* This function reads and counts variable special requirements like needing to be an integer, etc. */
void ReadAndCountVariableSpecialRequirements(unsigned int number_of_variables, int **variable_special_requirements, unsigned int *number_of_variables_required_to_be_standard, unsigned int *number_of_variables_required_to_be_integer, unsigned int *number_of_variables_required_to_be_binary, unsigned int *number_of_variables_required_to_be_unrestricted, unsigned int *number_of_variables_not_required_to_be_binary);

/* This function reads in the objective function's initial constant */
void ReadObjectiveFunctionInitialConstant(long long *objective_function_initial_constant);

/* This function reads the objective functions variable coefficients */
void ReadObjectiveFunctionCoefficientVector(unsigned int number_of_variables, long long ***objective_function_coefficient_vector);

/* This function reads the constraint inequality directions */
void ReadConstraintInequalityDirectionVector(unsigned int number_of_constraints, int **constraint_inequality_direction_vector);

/* This function reads the constraint constants */
void ReadConstraintConstantVector(unsigned int number_of_constraints, long long ***constraint_constant_vector);

/* This function reads in the constraint coefficient matrix */
void ReadConstraintCoefficientMatrix(unsigned int number_of_constraints, unsigned int number_of_variables, long long ****constraint_coefficient_matrix);

/*********************************************************************************/
/******************************** OPTIMAL VALUES *********************************/
/*********************************************************************************/

/* This function creates the arrays to store the optimal objective function value and the corresponding variable values */
void CreateOptimalObjectiveFunctionAndVariableValues(int maximization_problem, unsigned int number_of_variables, long long *optimal_objective_function_value, long long ***optimal_variable_values);

/*********************************************************************************/
/*********************** MIXED INTEGER LINEAR PROGRAMMING ************************/
/*********************************************************************************/

/* This function performs mixed integer linear programming */
int MixedIntegerLinearProgramming(int maximization_problem, unsigned int *number_of_constraints, unsigned int number_of_variables, unsigned int number_of_variables_required_to_be_integer, unsigned int number_of_variables_required_to_be_binary, unsigned int number_of_variables_not_required_to_be_binary, int *variable_special_requirements, long long *objective_function_initial_constant, long long **objective_function_coefficient_vector, int *constraint_inequality_direction_vector, long long **constraint_constant_vector, long long ***constraint_coefficient_matrix, long long *optimal_objective_function_value, long long **optimal_variable_values);

/* This function modifies the constraint inequality directions and then counts instances of each direction type */
void ModifyConstraintInequalityDirectionsAndCountDirections(unsigned int number_of_constraints, unsigned int number_of_variables_required_to_be_binary, long long **constraint_constant_vector, int *constraint_inequality_direction_vector, int **modified_constraint_inequality_direction_vector, unsigned int *number_of_less_than_or_equal_to_constraints, unsigned int *number_of_equal_to_constraints, unsigned int *number_of_greater_than_or_equal_to_constraints);

/* This function initializes the current and max tableau row and column counts */
void InitializeCurrentAndMaxTableauRowsAndColumnsCounts(unsigned int number_of_constraints, unsigned int number_of_variables, unsigned int number_of_variables_required_to_be_binary, unsigned int number_of_slack_surplus_variables, unsigned int number_of_artificial_variables, unsigned int *tableau_current_size, unsigned int *tableau_max_size);

/* This function creates the simplex tableau matrix */
void CreateTableauMatrix(int maximization_problem, unsigned int *number_of_constraints, unsigned int number_of_variables, unsigned int number_of_variables_required_to_be_binary, unsigned int number_of_slack_surplus_variables, unsigned int number_of_artificial_variables, unsigned int *tableau_max_size, int *variable_special_requirements, long long *objective_function_initial_constant, long long **objective_function_coefficient_vector, int *constraint_inequality_direction_vector, long long **constraint_constant_vector, long long ***constraint_coefficient_matrix, long long ****tableau_matrix);

/* This function creates the basic variables which tell which basis we currently are in */
void CreateBasicVariables(unsigned int number_of_constraints, unsigned int number_of_variables, unsigned int number_of_slack_surplus_variables, unsigned int number_of_artificial_variables, long long ***tableau_matrix, unsigned int **basic_variables);

/* This function creates the basic feasible solution to keep track of the best variable values in the current basis */
void CreateBasicFeasibleSolution(unsigned int number_of_constraints, unsigned int number_of_variables, unsigned int number_of_slack_surplus_variables, unsigned int number_of_artificial_variables, long long ***tableau_matrix, unsigned int *basic_variables, long long ***basic_feasible_solution);

/* This function writes to disk the initial MILP array values */
void WriteInitialLPValues(unsigned int number_of_constraints, unsigned int number_of_variables, unsigned int number_of_slack_surplus_variables, unsigned int number_of_artificial_variables, unsigned int *tableau_current_size, unsigned int *basic_variables, long long **basic_feasible_solution, long long ***tableau_matrix);

/* This function prints the initial constraint, variable, etc. counts */
void PrintInitialCounts(unsigned int number_of_constraints, unsigned int number_of_variables, unsigned int number_of_less_than_or_equal_to_constraints, unsigned int number_of_equal_to_constraints, unsigned int number_of_greater_than_or_equal_to_constraints, unsigned int number_of_slack_surplus_variables, unsigned int number_of_artificial_variables, unsigned int *tableau_current_size);

/* This function writes to disk the final LP array values */
void WriteFinalLPValues(unsigned int number_of_constraints, unsigned int number_of_variables, unsigned int number_of_slack_surplus_variables, unsigned int number_of_artificial_variables, unsigned int *tableau_current_size, unsigned int *basic_variables, long long **basic_feasible_solution, long long ***tableau_matrix);

/* This function writes to disk the final MILP array values */
void WriteFinalMILPValues(unsigned int number_of_variables, long long *optimal_objective_function_value, long long **optimal_variable_values);

/* This function updates the optimal objective function and variable solution */
void UpdateOptimalObjectiveFunctionAndVariableSolution(unsigned int number_of_variables, long long **basic_feasible_solution, long long *objective_function_initial_constant, long long **objective_function_coefficient_vector, long long *optimal_objective_function_value, long long **optimal_variable_values);

/*********************************************************************************/
/*********************************** SIMPLEX *************************************/
/*********************************************************************************/

/* This function finds the optimal solution for the given variables and constraints */
int SimplexAlgorithm(unsigned int number_of_constraints, unsigned int number_of_variables, unsigned int number_of_slack_surplus_variables, unsigned int *number_of_artificial_variables, unsigned int *basic_variables, long long **basic_feasible_solution, long long ***tableau_matrix, unsigned int *tableau_current_size);

/* This function transforms the tableau by removing artificial variables to obtain a basic feasible solution */
int SimplexPhase1(unsigned int number_of_constraints, unsigned int number_of_variables, unsigned int number_of_slack_surplus_variables, unsigned int *number_of_artificial_variables, unsigned int *basic_variables, long long ***tableau_matrix, unsigned int *tableau_current_size);

/* This function starts from an basic feasible solution and iterates toward the optimal solution */
int SimplexPhase2(unsigned int number_of_constraints, unsigned int number_of_variables, unsigned int number_of_slack_surplus_variables, unsigned int *basic_variables, long long ***tableau_matrix, unsigned int *tableau_current_size);

/* This function performs Gauss-Jordan Elimination on the pivot column */
void PivotColumnGaussJordanElimnation(unsigned int number_of_rows, unsigned int number_of_columns, unsigned int pivot_row_index, unsigned int pivot_col_index, long long *pivot_value, long long ***tableau_matrix);

/* This function updates the basic feasible solution */
void UpdateBasicFeasibleSolution(unsigned int number_of_total_variables, unsigned int number_of_constraints, unsigned int *basic_variables, long long **basic_feasible_solution, long long ***tableau_matrix);

/*********************************************************************************/
/******************************* BRANCH AND BOUND ********************************/
/*********************************************************************************/

/* This function initiates the branch and bound */
int BranchAndBoundMILP(int maximization_problem, unsigned int *number_of_constraints, unsigned int number_of_variables, unsigned int *number_of_slack_surplus_variables, unsigned int *number_of_artificial_variables, unsigned int **basic_variables, long long ***basic_feasible_solution, long long ****tableau_matrix, unsigned int *tableau_current_size, unsigned int *tableau_max_size, unsigned int number_of_variables_required_to_be_integer, unsigned int number_of_variables_required_to_be_binary, int *variable_special_requirements, long long *optimal_objective_function_value, long long **optimal_variable_values, long long **objective_function_coefficient_vector, int lp_error_code);

/* This function counts the number of variables that need to be integer or binary AND already are */
void CountNumberOfVariablesNeedingToBeIntegerOrBinaryThatAlreadyAre(unsigned int number_of_variables, unsigned int number_of_variables_required_to_be_integer, unsigned int number_of_variables_required_to_be_binary, int *variable_special_requirements, long long **basic_feasible_solution, unsigned int *number_of_variables_needing_to_be_integer_or_binary_currently_integer_or_binary, unsigned int *last_variable_that_still_needs_to_become_integer_or_binary_index);

/* This function recursively applies branch and bound */
int BranchAndBoundMILPRecursive(int maximization_problem, unsigned int *number_of_constraints, unsigned int number_of_variables, unsigned int *number_of_slack_surplus_variables, unsigned int *number_of_artificial_variables, unsigned int **basic_variables, long long ***basic_feasible_solution, long long ****tableau_matrix, unsigned int *tableau_current_size, unsigned int *tableau_max_size, unsigned int number_of_variables_required_to_be_integer, unsigned int number_of_variables_required_to_be_binary, int *variable_special_requirements, unsigned int last_variable_that_still_needs_to_become_integer_or_binary_index, double *best_mixed_integer_optimal_value_double, long long *best_mixed_integer_optimal_value, long long **best_mixed_integer_variable_values, int *best_milp_error_code, unsigned int recursion_level);

/* This function handles the less than or greater than branches of the recursive tree of branch and bound */
int BranchAndBoundMILPRecursiveLessOrGreaterThan(int maximization_problem, unsigned int *number_of_constraints, unsigned int number_of_variables, unsigned int *number_of_slack_surplus_variables, unsigned int *number_of_artificial_variables, unsigned int **basic_variables, long long ***basic_feasible_solution, long long ****tableau_matrix, unsigned int *tableau_current_size, unsigned int *tableau_max_size, unsigned int number_of_variables_required_to_be_integer, unsigned int number_of_variables_required_to_be_binary, int *variable_special_requirements, unsigned int last_variable_that_still_needs_to_become_integer_or_binary_index, double last_variable_that_still_needs_to_become_integer_value, double *best_mixed_integer_optimal_value_double, long long *best_mixed_integer_optimal_value, long long **best_mixed_integer_variable_values, int *best_milp_error_code, struct BranchAndBoundState branch_and_bound_state_old, int new_constraint_inequality_direction, unsigned int recursion_level);

/* This function creates cache counts and arrays to reset the active arrays at the end of each recursion level */
void CreateCacheCountsAndArraysForBranchAndBoundMILPRecursion(unsigned int tableau_current_size[2], unsigned int number_of_constraints, unsigned int number_of_variables, unsigned int number_of_slack_surplus_variables, unsigned int number_of_artificial_variables, unsigned int *basic_variables, long long **basic_feasible_solution, long long ***tableau_matrix, struct BranchAndBoundState *branch_and_bound_state_old);

/* This function checks if more recursion of MILP branch and bound is necessary */
int CheckIfMoreBranchAndBoundMILPRecursionIsNecessary(int maximization_problem, unsigned int *number_of_constraints, unsigned int number_of_variables, unsigned int *number_of_slack_surplus_variables, unsigned int *number_of_artificial_variables, unsigned int **basic_variables, long long ***basic_feasible_solution, long long ****tableau_matrix, unsigned int *tableau_current_size, unsigned int *tableau_max_size, unsigned int number_of_variables_required_to_be_integer, unsigned int number_of_variables_required_to_be_binary, int *variable_special_requirements, unsigned int last_variable_that_still_needs_to_become_integer_or_binary_index_recursive, double *best_mixed_integer_optimal_value_double, long long *best_mixed_integer_optimal_value, long long **best_mixed_integer_variable_values, int *best_milp_error_code, unsigned int recursion_level);

/* This function saves the best mixed integer optimal variables and values */
void SaveBestMixedIntegerOptimalVariablesAndValues(int maximization_problem, unsigned int number_of_constraints, unsigned int number_of_variables, long long **basic_feasible_solution, long long ***tableau_matrix, double *best_mixed_integer_optimal_value_double, long long *best_mixed_integer_optimal_value, long long **best_mixed_integer_variable_values, int *best_milp_error_code);

/* This function resets the counts and arrays during recursive branch and bound MILP */
void ResetBranchAndBoundMILPRecursiveCountsAndArrays(unsigned int number_of_variables, struct BranchAndBoundState branch_and_bound_state_old, unsigned int *tableau_current_size, unsigned int *number_of_constraints, unsigned int *number_of_slack_surplus_variables, unsigned int *number_of_artificial_variables, unsigned int **basic_variables, long long ***basic_feasible_solution, long long ****tableau_matrix);

/* This function adds a constraint to the previous optimal solution (warm-start) */
int AddConstraint(unsigned int *number_of_constraints, unsigned int number_of_variables, unsigned int *number_of_slack_surplus_variables, unsigned int *number_of_artificial_variables, unsigned int **basic_variables, long long ***basic_feasible_solution, long long ****tableau_matrix, unsigned int *tableau_current_size, unsigned int *tableau_max_size, unsigned int last_variable_that_still_needs_to_become_integer_or_binary_index, int new_constraint_inequality_direction, long long new_constraint_constant, unsigned int recursion_level);

/*********************************************************************************/
/******************************* HELPER FUNCTIONS ********************************/
/*********************************************************************************/

/* This function prints the optimal objective function value and the variable values */
void PrintOptimalResults(unsigned int number_of_variables, long long optimal_objective_function_value_numerator, long long optimal_objective_function_value_denominator, long long **optimal_variable_values);

/* This function reallocates more memory to the passed 1d array */
void Realloc1DUnsignedInt(unsigned int **array1d, unsigned int oldsize1d, unsigned int newsize1d, unsigned int initializer);

/* This function reallocates more memory to the passed 2d array */
void Realloc2DLongLong(long long ***array2d, unsigned int oldsize1d, unsigned int newsize1d);

/* This function reallocates more memory to the passed 3d array */
void Realloc3DLongLong(long long ****array3d, unsigned int oldsize1d, unsigned int newsize1d, unsigned int oldsize2d, unsigned int newsize2d);

/* This function adds long long rationals and keeps their numerators and denominators each as small as possible */
void LongLongRationalAddition(long long A_numerator, long long A_denominator, long long B_numerator, long long B_denominator, long long *C_numerator, long long *C_denominator);

/* This function divides long long rationals and keeps their numerators and denominators each as small as possible */
void LongLongRationalDivision(long long A_numerator, long long A_denominator, long long B_numerator, long long B_denominator, long long *C_numerator, long long *C_denominator);

/* This function finds the greatest common denominator between a and b */
long long GreastestCommonDenominator(long long a, long long b);

/*********************************************************************************/
/********************************** PROTOTYPES ***********************************/
/*********************************************************************************/

int main(int argc, char *argv[])
{
	clock_t start = clock(), diff;

	unsigned int i, k;
	int error_code = 0;

	/*********************************************************************************/
	/********************************* READ INPUTS ***********************************/
	/*********************************************************************************/

	/* Optimization problem type */
	int maximization_problem = 1;
	
	ReadObjectiveFunctionMaximizationOrMinimization(&maximization_problem);

	/* Sizes */
	unsigned int initial_number_of_constraints = 0, number_of_constraints = 0, number_of_variables = 0;

	ReadNumberOfConstraintsAndVariables(&initial_number_of_constraints, &number_of_constraints, &number_of_variables);

	/* Variable special requirements */
	unsigned int number_of_variables_required_to_be_standard = 0, number_of_variables_required_to_be_integer = 0, number_of_variables_required_to_be_binary = 0, number_of_variables_required_to_be_unrestricted = 0, number_of_variables_not_required_to_be_binary = 0;
	int *variable_special_requirements;

	ReadAndCountVariableSpecialRequirements(number_of_variables, &variable_special_requirements, &number_of_variables_required_to_be_standard, &number_of_variables_required_to_be_integer, &number_of_variables_required_to_be_binary, &number_of_variables_required_to_be_unrestricted, &number_of_variables_not_required_to_be_binary);

	/* Objective function */
	long long objective_function_initial_constant[2];

	ReadObjectiveFunctionInitialConstant(objective_function_initial_constant);

	long long **objective_function_coefficient_vector;

	ReadObjectiveFunctionCoefficientVector(number_of_variables, &objective_function_coefficient_vector);

	/* Constraints */
	int *constraint_inequality_direction_vector;

	ReadConstraintInequalityDirectionVector(number_of_constraints, &constraint_inequality_direction_vector);

	long long **constraint_constant_vector;

	ReadConstraintConstantVector(number_of_constraints, &constraint_constant_vector);

	long long ***constraint_coefficient_matrix;

	ReadConstraintCoefficientMatrix(number_of_constraints, number_of_variables, &constraint_coefficient_matrix);


	/*********************************************************************************/
	/******************************** OPTIMAL VALUES *********************************/
	/*********************************************************************************/

	/* Optimums */
	long long optimal_objective_function_value[2];

	long long **optimal_variable_values;

	CreateOptimalObjectiveFunctionAndVariableValues(maximization_problem, number_of_variables, optimal_objective_function_value, &optimal_variable_values);


	/*********************************************************************************/
	/*********************** MIXED INTEGER LINEAR PROGRAMMING ************************/
	/*********************************************************************************/

	error_code = MixedIntegerLinearProgramming(maximization_problem, &number_of_constraints, number_of_variables, number_of_variables_required_to_be_integer, number_of_variables_required_to_be_binary, number_of_variables_not_required_to_be_binary, variable_special_requirements, objective_function_initial_constant, objective_function_coefficient_vector, constraint_inequality_direction_vector, constraint_constant_vector, constraint_coefficient_matrix, optimal_objective_function_value, optimal_variable_values);

	if (error_code == 0)
	{
		printf("\n\n************************************************************************************************************\n");
		printf("************************************************************************************************************\n");
		printf("*********************************************** FINAL SOLUTION ******************************************\n");
		printf("************************************************************************************************************\n");
		printf("************************************************************************************************************\n\n\n");

		PrintOptimalResults(number_of_variables, optimal_objective_function_value[0], optimal_objective_function_value[1], optimal_variable_values);
	}
	else if (error_code == 1)
	{
		printf("Problem is infeasible!\n");
	}
	else if (error_code == 2)
	{
		printf("Problem is unbounded!\n");
	}
	else
	{
		printf("How did I get an error_code > 2?\n");
	}

	/*********************************************************************************/
	/****************************** FREE DYNAMIC MEMORY ******************************/
	/*********************************************************************************/

	for (k = 0; k < 2; k++)
	{
		for (i = 0; i < initial_number_of_constraints; i++)
		{
			free(constraint_coefficient_matrix[k][i]);
		} // end of i loop
		free(optimal_variable_values[k]);
		free(constraint_coefficient_matrix[k]);
		free(constraint_constant_vector[k]);
		free(objective_function_coefficient_vector[k]);
	} // end of k loop
	free(optimal_variable_values);
	free(constraint_coefficient_matrix);
	free(constraint_constant_vector);
	free(constraint_inequality_direction_vector);
	free(objective_function_coefficient_vector);
	free(variable_special_requirements);

	/* Calculate execution time */
	diff = clock() - start;
	unsigned long long microsec = diff * 1000000 / CLOCKS_PER_SEC;
	printf("\nmain: Time taken %llu microseconds\n", microsec);

	return 0;
} // end of main

/*********************************************************************************/
/*********************************** FUNCTIONS ***********************************/
/*********************************************************************************/

/*********************************************************************************/
/********************************* READ INPUTS ***********************************/
/*********************************************************************************/

/* This function reads in if the objective function is going to be a maximization or minimization problem */
void ReadObjectiveFunctionMaximizationOrMinimization(int *maximization_problem)
{
	int systemreturn;

	FILE *infile_objective_function_maximization = fopen("inputs/objective_function_maximization.txt", "r"); // read only
	if (infile_objective_function_maximization == NULL)
	{
		printf("ReadObjectiveFunctionMaximizationOrMinimization: Unable to open inputs/objective_function_maximization.txt\n");
	}
	else
	{
		systemreturn = fscanf(infile_objective_function_maximization, "%d", &(*maximization_problem));
		if (systemreturn == -1)
		{
			printf("ReadObjectiveFunctionMaximizationOrMinimization: Failed reading inputs/objective_function_maximization.txt\n");
		}
		fclose(infile_objective_function_maximization);
	}
	
	return;
} // end of ReadObjectiveFunctionMaximizationOrMinimization function

/* This function reads the number of constraints and variables */
void ReadNumberOfConstraintsAndVariables(unsigned int *initial_number_of_constraints, unsigned int *number_of_constraints, unsigned int *number_of_variables)
{
	int systemreturn;

	FILE *infile_number_of_constraints = fopen("inputs/number_of_constraints.txt", "r"); // read only
	if (infile_number_of_constraints == NULL)
	{
		printf("ReadNumberOfConstraintsAndVariables: Unable to open inputs/number_of_constraints.txt\n");
	}
	else
	{
		systemreturn = fscanf(infile_number_of_constraints, "%u", &(*number_of_constraints));
		if (systemreturn == -1)
		{
			printf("ReadNumberOfConstraintsAndVariables: Failed reading inputs/number_of_constraints.txt\n");
		}
		(*initial_number_of_constraints) = (*number_of_constraints);
		fclose(infile_number_of_constraints);
	}

	FILE *infile_number_of_variables = fopen("inputs/number_of_variables.txt", "r"); // read only
	if (infile_number_of_variables == NULL)
	{
		printf("ReadNumberOfConstraintsAndVariables: Unable to open inputs/number_of_variables.txt\n");
	}
	else
	{
		systemreturn = fscanf(infile_number_of_variables, "%u", &(*number_of_variables));
		if (systemreturn == -1)
		{
			printf("ReadNumberOfConstraintsAndVariables: Failed reading inputs/number_of_variables.txt\n");
		}
		fclose(infile_number_of_variables);
	}
	
	return;
} // end of ReadNumberOfConstraintsAndVariables function

/* This function reads and counts variable special requirements like needing to be an integer, etc. */
void ReadAndCountVariableSpecialRequirements(unsigned int number_of_variables, int **variable_special_requirements, unsigned int *number_of_variables_required_to_be_standard, unsigned int *number_of_variables_required_to_be_integer, unsigned int *number_of_variables_required_to_be_binary, unsigned int *number_of_variables_required_to_be_unrestricted, unsigned int *number_of_variables_not_required_to_be_binary)
{
	unsigned int i;
	int systemreturn;

	(*variable_special_requirements) = malloc(sizeof(int) * number_of_variables);

	for (i = 0; i < number_of_variables; i++)
	{
		(*variable_special_requirements)[i] = 0;
	} // end of i loop

	FILE *infile_variable_special_requirements = fopen("inputs/variable_special_requirements.txt", "r"); // read only
	if (infile_variable_special_requirements == NULL)
	{
		printf("ReadAndCountVariableSpecialRequirements: Unable to open inputs/variable_special_requirements.txt\n");
	}
	else
	{
		for (i = 0; i < number_of_variables; i++)
		{
			systemreturn = fscanf(infile_variable_special_requirements, "%u\n", &(*variable_special_requirements)[i]);
			if (systemreturn == -1)
			{
				printf("ReadAndCountVariableSpecialRequirements: Failed reading inputs/variable_special_requirements.txt\n");
			}

			if ((*variable_special_requirements)[i] == 0)
			{
				(*number_of_variables_required_to_be_standard)++;
			}
			else if ((*variable_special_requirements)[i] == 1)
			{
				(*number_of_variables_required_to_be_integer)++;
			}
			else if ((*variable_special_requirements)[i] == 2)
			{
				(*number_of_variables_required_to_be_binary)++;
			}
			else if ((*variable_special_requirements)[i] == 3)
			{
				(*number_of_variables_required_to_be_unrestricted)++;
			}
			else
			{
				printf("ReadAndCountVariableSpecialRequirements: Variable %d has an unknown type of %u\n", i, (*variable_special_requirements)[i]);
			}
		} // end of i loop
		fclose(infile_variable_special_requirements);
	}

	(*number_of_variables_not_required_to_be_binary) = number_of_variables - (*number_of_variables_required_to_be_binary);
	
	return;
}// end of ReadAndCountVariableSpecialRequirements function

/* This function reads in the objective function's initial constant */
void ReadObjectiveFunctionInitialConstant(long long *objective_function_initial_constant)
{
	int systemreturn;

	objective_function_initial_constant[0] = 0;
	objective_function_initial_constant[1] = 1;

	FILE *infile_objective_function_initial_constant_numerator = fopen("inputs/objective_function_initial_constant_numerator.txt", "r"); // read only
	if (infile_objective_function_initial_constant_numerator == NULL)
	{
		printf("ReadObjectiveFunctionInitialConstant: Unable to open inputs/objective_function_initial_constant_numerator.txt\n");
	}
	else
	{
		systemreturn = fscanf(infile_objective_function_initial_constant_numerator, "%lld", &objective_function_initial_constant[0]);
		if (systemreturn == -1)
		{
			printf("ReadObjectiveFunctionInitialConstant: Failed reading inputs/objective_function_initial_constant_numerator.txt\n");
		}
		fclose(infile_objective_function_initial_constant_numerator);
	}

	FILE *infile_objective_function_initial_constant_denominator = fopen("inputs/objective_function_initial_constant_denominator.txt", "r"); // read only
	if (infile_objective_function_initial_constant_denominator == NULL)
	{
		printf("ReadObjectiveFunctionInitialConstant: Unable to open inputs/objective_function_initial_constant_denominator.txt\n");
	}
	else
	{
		systemreturn = fscanf(infile_objective_function_initial_constant_denominator, "%lld", &objective_function_initial_constant[1]);
		if (systemreturn == -1)
		{
			printf("ReadObjectiveFunctionInitialConstant: Failed reading inputs/objective_function_initial_constant_denominator.txt\n");
		}
		fclose(infile_objective_function_initial_constant_denominator);
	}
	
	return;
} // end of ReadObjectiveFunctionInitialConstant function

/* This function reads the objective functions variable coefficients */
void ReadObjectiveFunctionCoefficientVector(unsigned int number_of_variables, long long ***objective_function_coefficient_vector)
{
	unsigned int i, k;
	int systemreturn;

	(*objective_function_coefficient_vector) = malloc(sizeof(long long*) * 2);
	for (k = 0; k < 2; k++)
	{
		(*objective_function_coefficient_vector)[k] = malloc(sizeof(long long) * number_of_variables);
		for (i = 0; i < number_of_variables; i++)
		{
			(*objective_function_coefficient_vector)[k][i] = 0;
		} // end of i loop
	} // end of k loop

	FILE *infile_objective_function_coefficient_vector_numerator = fopen("inputs/objective_function_coefficient_vector_numerator.txt", "r"); // read only
	if (infile_objective_function_coefficient_vector_numerator == NULL)
	{
		printf("ReadObjectiveFunctionCoefficientVector: Unable to open inputs/objective_function_coefficient_vector_numerator.txt\n");
	}
	else
	{
		for (i = 0; i < number_of_variables; i++)
		{
			systemreturn = fscanf(infile_objective_function_coefficient_vector_numerator, "%lld\n", &(*objective_function_coefficient_vector)[0][i]);
			if (systemreturn == -1)
			{
				printf("ReadObjectiveFunctionCoefficientVector: Failed reading inputs/objective_function_coefficient_vector_numerator.txt\n");
			}
		} // end of i loop
		fclose(infile_objective_function_coefficient_vector_numerator);
	}

	FILE *infile_objective_function_coefficient_vector_denominator = fopen("inputs/objective_function_coefficient_vector_denominator.txt", "r"); // read only
	if (infile_objective_function_coefficient_vector_denominator == NULL)
	{
		printf("ReadObjectiveFunctionCoefficientVector: Unable to open inputs/objective_function_coefficient_vector_denominator.txt\n");
	}
	else
	{
		for (i = 0; i < number_of_variables; i++)
		{
			systemreturn = fscanf(infile_objective_function_coefficient_vector_denominator, "%lld\n", &(*objective_function_coefficient_vector)[1][i]);
			if (systemreturn == -1)
			{
				printf("ReadObjectiveFunctionCoefficientVector: Failed reading inputs/objective_function_coefficient_vector_denominator.txt\n");
			}
		} // end of i loop
		fclose(infile_objective_function_coefficient_vector_denominator);
	}
	
	return;
} // end of ReadObjectiveFunctionCoefficientVector function

/* This function reads the constraint inequality directions */
void ReadConstraintInequalityDirectionVector(unsigned int number_of_constraints, int **constraint_inequality_direction_vector)
{
	unsigned int i;
	int systemreturn;

	(*constraint_inequality_direction_vector) = malloc(sizeof(int) * number_of_constraints);
	for (i = 0; i < number_of_constraints; i++)
	{
		(*constraint_inequality_direction_vector)[i] = 0;
	} // end of i loop

	FILE *infile_constraint_inequality_direction_vector = fopen("inputs/constraint_inequality_direction_vector.txt", "r"); // read only
	if (infile_constraint_inequality_direction_vector == NULL)
	{
		printf("ReadConstraintInequalityDirectionVector: Unable to open inputs/constraint_inequality_direction_vector.txt\n");
	}
	else
	{
		for (i = 0; i < number_of_constraints; i++)
		{
			systemreturn = fscanf(infile_constraint_inequality_direction_vector, "%d\n", &(*constraint_inequality_direction_vector)[i]);
			if (systemreturn == -1)
			{
				printf("ReadConstraintInequalityDirectionVector: Failed reading inputs/constraint_inequality_direction_vector.txt\n");
			}
		} // end of i loop
		fclose(infile_constraint_inequality_direction_vector);
	}
	
	return;
} // end of ReadConstraintInequalityDirectionVector function

/* This function reads the constraint constants */
void ReadConstraintConstantVector(unsigned int number_of_constraints, long long ***constraint_constant_vector)
{
	unsigned int i, k;
	int systemreturn;

	(*constraint_constant_vector) = malloc(sizeof(long long*) * 2);
	for (k = 0; k < 2; k++)
	{
		(*constraint_constant_vector)[k] = malloc(sizeof(long long) * number_of_constraints);
		for (i = 0; i < number_of_constraints; i++)
		{
			(*constraint_constant_vector)[k][i] = 0;
		} // end of i loop
	} // end of k loop

	FILE *infile_constraint_constant_vector_numerator = fopen("inputs/constraint_constant_vector_numerator.txt", "r"); // read only
	if (infile_constraint_constant_vector_numerator == NULL)
	{
		printf("ReadConstraintConstantVector: Unable to open inputs/constraint_constant_vector_numerator.txt\n");
	}
	else
	{
		for (i = 0; i < number_of_constraints; i++)
		{
			systemreturn = fscanf(infile_constraint_constant_vector_numerator, "%lld\n", &(*constraint_constant_vector)[0][i]);
			if (systemreturn == -1)
			{
				printf("ReadConstraintConstantVector: Failed reading inputs/constraint_constant_vector_numerator.txt\n");
			}
		} // end of i loop
		fclose(infile_constraint_constant_vector_numerator);
	}

	FILE *infile_constraint_constant_vector_denominator = fopen("inputs/constraint_constant_vector_denominator.txt", "r"); // read only
	if (infile_constraint_constant_vector_denominator == NULL)
	{
		printf("ReadConstraintConstantVector: Unable to open inputs/constraint_constant_vector_denominator.txt\n");
	}
	else
	{
		for (i = 0; i < number_of_constraints; i++)
		{
			systemreturn = fscanf(infile_constraint_constant_vector_denominator, "%lld\n", &(*constraint_constant_vector)[1][i]);
			if (systemreturn == -1)
			{
				printf("ReadConstraintConstantVector: Failed reading inputs/constraint_constant_vector_denominator.txt\n");
			}
		} // end of i loop
		fclose(infile_constraint_constant_vector_denominator);
	}
	
	return;
} // end of ReadConstraintConstantVector function

/* This function reads in the constraint coefficient matrix */
void ReadConstraintCoefficientMatrix(unsigned int number_of_constraints, unsigned int number_of_variables, long long ****constraint_coefficient_matrix)
{
	unsigned int i, j, k;
	int systemreturn;

	(*constraint_coefficient_matrix) = malloc(sizeof(long long*) * 2);
	for (k = 0; k < 2; k++)
	{
		(*constraint_coefficient_matrix)[k] = malloc(sizeof(long long*) * number_of_constraints);
		for (i = 0; i < number_of_constraints; i++)
		{
			(*constraint_coefficient_matrix)[k][i] = malloc(sizeof(long long) * number_of_variables);
			for (j = 0; j < number_of_variables; j++)
			{
				(*constraint_coefficient_matrix)[k][i][j] = 0;
			} // end of j loop
		} // end of i loop
	} // end of k loop

	FILE *infile_constraint_coefficient_matrix_numerator = fopen("inputs/constraint_coefficient_matrix_numerator.txt", "r"); // read only
	if (infile_constraint_coefficient_matrix_numerator == NULL)
	{
		printf("ReadConstraintCoefficientMatrix: Unable to open inputs/constraint_coefficient_matrix_numerator.txt\n");
	}
	else
	{
		for (i = 0; i < number_of_constraints; i++)
		{
			for (j = 0; j < number_of_variables; j++)
			{
				systemreturn = fscanf(infile_constraint_coefficient_matrix_numerator, "%lld\n", &(*constraint_coefficient_matrix)[0][i][j]);
				if (systemreturn == -1)
				{
					printf("ReadConstraintCoefficientMatrix: Failed reading inputs/constraint_coefficient_matrix_numerator.txt\n");
				}
			} // end of j loop
		} // end of i loop
		fclose(infile_constraint_coefficient_matrix_numerator);
	}

	FILE *infile_constraint_coefficient_matrix_denominator = fopen("inputs/constraint_coefficient_matrix_denominator.txt", "r"); // read only
	if (infile_constraint_coefficient_matrix_denominator == NULL)
	{
		printf("ReadConstraintCoefficientMatrix: Unable to open inputs/constraint_coefficient_matrix_denominator.txt\n");
	}
	else
	{
		for (i = 0; i < number_of_constraints; i++)
		{
			for (j = 0; j < number_of_variables; j++)
			{
				systemreturn = fscanf(infile_constraint_coefficient_matrix_denominator, "%lld\n", &(*constraint_coefficient_matrix)[1][i][j]);
				if (systemreturn == -1)
				{
					printf("ReadConstraintCoefficientMatrix: Failed reading inputs/constraint_coefficient_matrix_denominator.txt\n");
				}
			} // end of j loop
		} // end of i loop
		fclose(infile_constraint_coefficient_matrix_denominator);
	}
	
	return;
} // end of ReadConstraintCoefficientMatrix function

/*********************************************************************************/
/******************************** OPTIMAL VALUES *********************************/
/*********************************************************************************/

/* This function creates the arrays to store the optimal objective function value and the corresponding variable values */
void CreateOptimalObjectiveFunctionAndVariableValues(int maximization_problem, unsigned int number_of_variables, long long *optimal_objective_function_value, long long ***optimal_variable_values)
{
	unsigned int i, k;

	if (maximization_problem == 1)
	{
		optimal_objective_function_value[0] = -LLONG_MAX;
	}
	else
	{
		optimal_objective_function_value[0] = LLONG_MAX;
	}
	optimal_objective_function_value[1] = 1;

	(*optimal_variable_values) = malloc(sizeof(long long*) * 2);
	for (k = 0; k < 2; k++)
	{
		(*optimal_variable_values)[k] = malloc(sizeof(long long) * number_of_variables);
	} // end of k loop

	for (i = 0; i < number_of_variables; i++)
	{
		(*optimal_variable_values)[0][i] = 0;
		(*optimal_variable_values)[1][i] = 1;
	} // end of i loop
	
	return;
} // end of CreateOptimalObjectiveFunctionAndVariableValues function

/*********************************************************************************/
/*********************** MIXED INTEGER LINEAR PROGRAMMING ************************/
/*********************************************************************************/

/* This function performs mixed integer linear programming */
int MixedIntegerLinearProgramming(int maximization_problem, unsigned int *number_of_constraints, unsigned int number_of_variables, unsigned int number_of_variables_required_to_be_integer, unsigned int number_of_variables_required_to_be_binary, unsigned int number_of_variables_not_required_to_be_binary, int *variable_special_requirements, long long *objective_function_initial_constant, long long **objective_function_coefficient_vector, int *constraint_inequality_direction_vector, long long **constraint_constant_vector, long long ***constraint_coefficient_matrix, long long *optimal_objective_function_value, long long **optimal_variable_values)
{
	int i, k, error_code = 0;

	printf("\n\n************************************************************************************************************\n");
	printf("************************************************************************************************************\n");
	printf("*************************************************** MILP ***************************************************\n");
	printf("************************************************************************************************************\n");
	printf("************************************************************************************************************\n\n\n");

	unsigned int number_of_less_than_or_equal_to_constraints = 0, number_of_equal_to_constraints = 0, number_of_greater_than_or_equal_to_constraints = 0;

	int *modified_constraint_inequality_direction_vector;

	ModifyConstraintInequalityDirectionsAndCountDirections((*number_of_constraints), number_of_variables_required_to_be_binary, constraint_constant_vector, constraint_inequality_direction_vector, &modified_constraint_inequality_direction_vector, &number_of_less_than_or_equal_to_constraints, &number_of_equal_to_constraints, &number_of_greater_than_or_equal_to_constraints);

	unsigned int number_of_slack_surplus_variables = number_of_less_than_or_equal_to_constraints + number_of_greater_than_or_equal_to_constraints, number_of_artificial_variables = number_of_equal_to_constraints + number_of_greater_than_or_equal_to_constraints;

	unsigned int tableau_current_size[2];
	tableau_current_size[0] = 0;
	tableau_current_size[1] = 0;
	
	unsigned int tableau_max_size[2];
	tableau_max_size[0] = 0;
	tableau_max_size[1] = 0;

	InitializeCurrentAndMaxTableauRowsAndColumnsCounts((*number_of_constraints), number_of_variables, number_of_variables_required_to_be_binary, number_of_slack_surplus_variables, number_of_artificial_variables, &tableau_current_size[0], &tableau_max_size[0]);

	long long ***tableau_matrix;
	CreateTableauMatrix(maximization_problem, number_of_constraints, number_of_variables, number_of_variables_required_to_be_binary, number_of_slack_surplus_variables, number_of_artificial_variables, tableau_max_size, variable_special_requirements, objective_function_initial_constant, objective_function_coefficient_vector, modified_constraint_inequality_direction_vector, constraint_constant_vector, constraint_coefficient_matrix, &tableau_matrix);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	unsigned int *basic_variables;

	CreateBasicVariables((*number_of_constraints), number_of_variables, number_of_slack_surplus_variables, number_of_artificial_variables, tableau_matrix, &basic_variables);

	long long **basic_feasible_solution;

	CreateBasicFeasibleSolution((*number_of_constraints), number_of_variables, number_of_slack_surplus_variables, number_of_artificial_variables, tableau_matrix, basic_variables, &basic_feasible_solution);

	WriteInitialLPValues((*number_of_constraints), number_of_variables, number_of_slack_surplus_variables, number_of_artificial_variables, tableau_current_size, basic_variables, basic_feasible_solution, tableau_matrix);

	PrintInitialCounts((*number_of_constraints), number_of_variables, number_of_less_than_or_equal_to_constraints, number_of_equal_to_constraints, number_of_greater_than_or_equal_to_constraints, number_of_slack_surplus_variables, number_of_artificial_variables, tableau_current_size);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/* Solve LP Relaxation */

	/* This function finds the optimal solution for the given variables and constraints */
	error_code = SimplexAlgorithm((*number_of_constraints), number_of_variables, number_of_slack_surplus_variables, &number_of_artificial_variables, basic_variables, basic_feasible_solution, tableau_matrix, &tableau_current_size[0]);

	WriteFinalLPValues((*number_of_constraints), number_of_variables, number_of_slack_surplus_variables, number_of_artificial_variables, tableau_current_size, basic_variables, basic_feasible_solution, tableau_matrix);

	if (error_code == 0) // if LP is optimal then MILP will be either infeasible or optimal
	{
		printf("MixedIntegerLinearProgramming: LP is optimal!\n");

		UpdateOptimalObjectiveFunctionAndVariableSolution(number_of_variables, basic_feasible_solution, objective_function_initial_constant, objective_function_coefficient_vector, optimal_objective_function_value, optimal_variable_values);

		PrintOptimalResults(number_of_variables, optimal_objective_function_value[0], optimal_objective_function_value[1], optimal_variable_values);

		if (number_of_variables_required_to_be_integer + number_of_variables_required_to_be_binary == 0) // if there are NO variables required to be integer or binary
		{
			/* Do nothing, you don't need to proceed any farther since you have no required integer or binary variables */
		} // end of if there are NO variables required to be integer
		else if (number_of_variables_required_to_be_integer + number_of_variables_required_to_be_binary > 0) // if there is at least one variable required to be integer or binary
		{
			printf("\n\n************************************************************************************************************\n");
			printf("************************************************************************************************************\n");
			printf("*********************************************  BRANCH AND BOUND ********************************************\n");
			printf("************************************************************************************************************\n");
			printf("************************************************************************************************************\n\n\n");

			error_code = BranchAndBoundMILP(maximization_problem, number_of_constraints, number_of_variables, &number_of_slack_surplus_variables, &number_of_artificial_variables, &basic_variables, &basic_feasible_solution, &tableau_matrix, &tableau_current_size[0], &tableau_max_size[0], number_of_variables_required_to_be_integer, number_of_variables_required_to_be_binary, variable_special_requirements, optimal_objective_function_value, optimal_variable_values, objective_function_coefficient_vector, error_code);

			if (error_code == 0)
			{
				printf("MixedIntegerLinearProgramming: LP was optimal, MILP is optimal!\n");

				WriteFinalMILPValues(number_of_variables, optimal_objective_function_value, optimal_variable_values);
			}
			else if (error_code == 1)
			{
				printf("MixedIntegerLinearProgramming: LP was optimal, MILP is infeasible!\n");
			}
			else if (error_code == 2)
			{
				printf("MixedIntegerLinearProgramming: LP was optimal, MILP is unbounded!  This is IMPOSSIBLE since LP was optimal!\n");
			}
			else
			{
				printf("MixedIntegerLinearProgramming: WEIRD ERROR_CODE = %d\n", error_code);
			}
		} // end of if there is at least one variable required to be integer
	}
	else if (error_code == 1) // if LP is infeasible then MILP will be infeasible since it is a subset of a NULL set
	{
		printf("MixedIntegerLinearProgramming: LP is infeasible!\n");

		if (number_of_variables_required_to_be_integer + number_of_variables_required_to_be_binary > 0) // if there is at least one variable required to be integer or binary
		{
			printf("MixedIntegerLinearProgramming: An infeasible LP ALWAYS leads to an infeasible MILP!\n");
		}
	}
	else if (error_code == 2) // if LP is unbounded and all LP coefficients are rational then the MILP will either be infeasible or unbounded, however if LP has some irrational coefficients then MILP might be optimal instead
	{
		printf("MixedIntegerLinearProgramming: LP is unbounded!\n");

		if (number_of_variables_required_to_be_integer + number_of_variables_required_to_be_binary > 0) // if there is at least one variable required to be integer or binary
		{
			printf("\n\n************************************************************************************************************\n");
			printf("************************************************************************************************************\n");
			printf("*********************************************  BRANCH AND BOUND ********************************************\n");
			printf("************************************************************************************************************\n");
			printf("************************************************************************************************************\n\n\n");

			error_code = BranchAndBoundMILP(maximization_problem, number_of_constraints, number_of_variables, &number_of_slack_surplus_variables, &number_of_artificial_variables, &basic_variables, &basic_feasible_solution, &tableau_matrix, &tableau_current_size[0], &tableau_max_size[0], number_of_variables_required_to_be_integer, number_of_variables_required_to_be_binary, variable_special_requirements, optimal_objective_function_value, optimal_variable_values, objective_function_coefficient_vector, error_code);

			if (error_code == 0)
			{
				printf("MixedIntegerLinearProgramming: LP was unbounded, MILP is optimal!  There must be some irrational coefficients and/or binary variables!\n");
			}
			else if (error_code == 1)
			{
				printf("MixedIntegerLinearProgramming: LP was unbounded, MILP is infeasible!\n");
			}
			else if (error_code == 2)
			{
				printf("MixedIntegerLinearProgramming: LP was unbounded, MILP is unbounded!\n");
			}
			else
			{
				printf("MixedIntegerLinearProgramming: WEIRD ERROR_CODE = %d\n", error_code);
			}
		} // end of if there is at least one variable required to be integer or binary
	}
	else
	{
		printf("MixedIntegerLinearProgramming: WEIRD ERROR_CODE = %d\n", error_code);
	}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/* Free dynamic memory */
	for (k = 0; k < 2; k++)
	{
		free(basic_feasible_solution[k]);

		for (i = 0; i < tableau_max_size[0]; i++)
		{
			free(tableau_matrix[k][i]);
		} // end of i loop
		free(tableau_matrix[k]);
	} // end of k loop
	free(basic_feasible_solution);
	free(basic_variables);

	free(tableau_matrix);

	free(modified_constraint_inequality_direction_vector);

	return error_code;
} // end of MixedIntegerLinearProgramming function

/* This function modifies the constraint inequality directions and then counts instances of each direction type */
void ModifyConstraintInequalityDirectionsAndCountDirections(unsigned int number_of_constraints, unsigned int number_of_variables_required_to_be_binary, long long **constraint_constant_vector, int *constraint_inequality_direction_vector, int **modified_constraint_inequality_direction_vector, unsigned int *number_of_less_than_or_equal_to_constraints, unsigned int *number_of_equal_to_constraints, unsigned int *number_of_greater_than_or_equal_to_constraints)
{
	unsigned int i;

	(*modified_constraint_inequality_direction_vector) = malloc(sizeof(int) * (number_of_constraints + number_of_variables_required_to_be_binary));
	for (i = 0; i < number_of_constraints; i++)
	{
		(*modified_constraint_inequality_direction_vector)[i] = constraint_inequality_direction_vector[i];
	} // end of i loop

	for (i = number_of_constraints; i < number_of_constraints + number_of_variables_required_to_be_binary; i++)
	{
		(*modified_constraint_inequality_direction_vector)[i] = 1;
	} // end of i loop
	(*number_of_less_than_or_equal_to_constraints) = number_of_variables_required_to_be_binary;

	for (i = 0; i < number_of_constraints; i++)
	{
		if (constraint_inequality_direction_vector[i] == -1) // if greater than or equal to
		{
			if (constraint_constant_vector[0][i] < 0)
			{
				(*number_of_less_than_or_equal_to_constraints)++; // reverse sign
				(*modified_constraint_inequality_direction_vector)[i] = 1;
			}
			else
			{
				(*number_of_greater_than_or_equal_to_constraints)++;
			}
		} // end if greater than or equal to
		else if (constraint_inequality_direction_vector[i] == 0) // if equal to
		{
			(*number_of_equal_to_constraints)++;
		} // end if equal to
		else if (constraint_inequality_direction_vector[i] == 1) // if less than or equal to
		{
			if (constraint_constant_vector[0][i] < 0)
			{
				(*number_of_greater_than_or_equal_to_constraints)++; // reverse sign
				(*modified_constraint_inequality_direction_vector)[i] = -1;
			}
			else
			{
				(*number_of_less_than_or_equal_to_constraints)++;
			}
		} // end if less than or equal to
	} // end of i loop
	
	return;
} // end of ModifyConstraintInequalityDirectionsAndCountDirections function

/* This function initializes the current and max tableau row and column counts */
void InitializeCurrentAndMaxTableauRowsAndColumnsCounts(unsigned int number_of_constraints, unsigned int number_of_variables, unsigned int number_of_variables_required_to_be_binary, unsigned int number_of_slack_surplus_variables, unsigned int number_of_artificial_variables, unsigned int *tableau_current_size, unsigned int *tableau_max_size)
{
	if (number_of_artificial_variables > 0)
	{
		tableau_max_size[0] = number_of_constraints + number_of_variables_required_to_be_binary + 2;
	}
	else
	{
		tableau_max_size[0] = number_of_constraints + number_of_variables_required_to_be_binary + 1;
	}
	tableau_max_size[1] = number_of_variables + number_of_slack_surplus_variables + number_of_artificial_variables + 1;
	
	tableau_current_size[0] = tableau_max_size[0];
	tableau_current_size[1] = tableau_max_size[1];
	
	return;
} // end of InitializeCurrentAndMaxTableauRowsAndColumnsCounts

/* This function creates the simplex tableau matrix */
void CreateTableauMatrix(int maximization_problem, unsigned int *number_of_constraints, unsigned int number_of_variables, unsigned int number_of_variables_required_to_be_binary, unsigned int number_of_slack_surplus_variables, unsigned int number_of_artificial_variables, unsigned int *tableau_max_size, int *variable_special_requirements, long long *objective_function_initial_constant, long long **objective_function_coefficient_vector, int *constraint_inequality_direction_vector, long long **constraint_constant_vector, long long ***constraint_coefficient_matrix, long long ****tableau_matrix)
{
	/*
		Tableau generic format is:

		con\var	|	b	x1	x2	...	xn	s1		s2		...	sq		A1	A2	...		Ar
		------------------------------------------------------------------------------------------------------------------------
		con1	|	b1	x11	x11	...	x1n	0/1/-1	0/1/-1	...	0/1/-1	0/1	0/1	...		0/1
		con2	|	b2	x21	x22	...	x2n	0/1/-1	0/1/-1	...	0/1/-1	0/1	0/1	...		0/1
		.		|	.	.	.	...	.	.		.		...	.		.	.	...		.
		.		|	.	.	.	...	.	.		.		...	.		.	.	...		.
		.		|	.	.	.	...	.	.		.		...	.		.	.	...		.
		conm	|	bm	xm1	xm2	...	xmn	0/1/-1	0/1/-1	...	0/1/-1	0/1	0/1	...		0/1
		z		|	z0	-c1	-c2	...	0	0		0		...	0		0	0	...		0
		zA		|	0	0	0	...	0	0		0		...	0		-1	-1	...		-1
	*/

	unsigned int i, j, k;

	(*tableau_matrix) = malloc(sizeof(long long*) * 2);
	for (k = 0; k < 2; k++)
	{
		(*tableau_matrix)[k] = malloc(sizeof(long long*) * tableau_max_size[0]);
		for (i = 0; i < tableau_max_size[0]; i++)
		{
			(*tableau_matrix)[k][i] = malloc(sizeof(long long) * tableau_max_size[1]);
		} // end of i loop
	} // end of k loop

	for (i = 0; i < tableau_max_size[0]; i++)
	{
		for (j = 0; j < tableau_max_size[1]; j++)
		{
			(*tableau_matrix)[0][i][j] = 0;
			(*tableau_matrix)[1][i][j] = 1;
		} // end of j loop
	} // end of i loop

	/* First add b constant vector */
	for (i = 0; i < (*number_of_constraints); i++)
	{
		if (constraint_constant_vector[0][i] < 0)
		{
			(*tableau_matrix)[0][i][0] = -constraint_constant_vector[0][i];
		}
		else
		{
			(*tableau_matrix)[0][i][0] = constraint_constant_vector[0][i];
		}
		(*tableau_matrix)[1][i][0] = constraint_constant_vector[1][i];
	} // end of i loop

	for (i = (*number_of_constraints); i < (*number_of_constraints) + number_of_variables_required_to_be_binary; i++)
	{
		for (k = 0; k < 2; k++)
		{
			(*tableau_matrix)[k][i][0] = 1;
		} // end of k loop
	} // end of i loop

	/* Next add A matrix */
	for (i = 0; i < (*number_of_constraints); i++)
	{
		if (constraint_constant_vector[0][i] < 0)
		{
			for (j = 0; j < number_of_variables; j++)
			{
				(*tableau_matrix)[0][i][j + 1] = -constraint_coefficient_matrix[0][i][j]; // multiply row by -1
				(*tableau_matrix)[1][i][j + 1] = constraint_coefficient_matrix[1][i][j];
			} // end of j loop
		}
		else
		{
			for (j = 0; j < number_of_variables; j++)
			{
				(*tableau_matrix)[0][i][j + 1] = constraint_coefficient_matrix[0][i][j];
				(*tableau_matrix)[1][i][j + 1] = constraint_coefficient_matrix[1][i][j];
			} // end of j loop
		}
	} // end of i loop

	i = (*number_of_constraints);
	for (j = 0; j < number_of_variables; j++)
	{
		if (variable_special_requirements[j] == 2) // if jth variable is required to be binary
		{
			for (k = 0; k < 2; k++)
			{
				(*tableau_matrix)[k][i][j + 1] = 1;
			} // end of k loop
			i++;
		}
	} // end of i loop

	/* Next add slack/surplus variable matrix */
	j = 0;
	for (i = 0; i < (*number_of_constraints) + number_of_variables_required_to_be_binary; i++)
	{
		if (constraint_inequality_direction_vector[i] == 1) // if less than or equal to
		{
			(*tableau_matrix)[0][i][number_of_variables + 1 + j] = 1; // slack
			(*tableau_matrix)[1][i][number_of_variables + 1 + j] = 1; // slack
			j++;
		} // end of if less than or equal to
		else if (constraint_inequality_direction_vector[i] == -1) // if greater than or equal to
		{
			(*tableau_matrix)[0][i][number_of_variables + 1 + j] = -1; // surplus
			(*tableau_matrix)[1][i][number_of_variables + 1 + j] = 1; // surplus
			j++;
		} // end of if greater than or equal to
	} // end of i loop

	/* Next add artificial variable matrix */
	j = 0;
	for (i = 0; i < (*number_of_constraints); i++)
	{
		if (constraint_inequality_direction_vector[i] != 1) // if NOT less than or equal to
		{
			for (k = 0; k < 2; k++)
			{
				(*tableau_matrix)[k][i][number_of_variables + number_of_slack_surplus_variables + 1 + j] = 1;
			} // end of k loop
			j++;
		} // end of if equal to OR greater than or equal to
	} // end of i loop

	/* Now we can update the number of constraints since we are done adding the extra binary variable information */
	(*number_of_constraints) += number_of_variables_required_to_be_binary;

	/* Next add original -cT vector */
	if (maximization_problem == 1)
	{
		for (j = 0; j < number_of_variables; j++)
		{
			(*tableau_matrix)[0][(*number_of_constraints)][j + 1] = -objective_function_coefficient_vector[0][j];
			(*tableau_matrix)[1][(*number_of_constraints)][j + 1] = objective_function_coefficient_vector[1][j];
		} // end of j loop
	}
	else
	{
		for (j = 0; j < number_of_variables; j++)
		{
			(*tableau_matrix)[0][(*number_of_constraints)][j + 1] = objective_function_coefficient_vector[0][j];
			(*tableau_matrix)[1][(*number_of_constraints)][j + 1] = objective_function_coefficient_vector[1][j];
		} // end of j loop
	}

	/* Next add artificial -cT vector */
	if (number_of_artificial_variables > 0)
	{
		for (j = 0; j < number_of_artificial_variables; j++)
		{
			(*tableau_matrix)[0][(*number_of_constraints) + 1][number_of_variables + number_of_slack_surplus_variables + 1 + j] = -1;
			(*tableau_matrix)[1][(*number_of_constraints) + 1][number_of_variables + number_of_slack_surplus_variables + 1 + j] = 1;
		} // end of j loop
	}

	/* Lastly add initial z */
	(*tableau_matrix)[0][(*number_of_constraints)][0] = objective_function_initial_constant[0];
	(*tableau_matrix)[1][(*number_of_constraints)][0] = objective_function_initial_constant[1];
	
	return;
} // end of CreateTableauMatrix function

/* This function creates the basic variables which tell which basis we currently are in */
void CreateBasicVariables(unsigned int number_of_constraints, unsigned int number_of_variables, unsigned int number_of_slack_surplus_variables, unsigned int number_of_artificial_variables, long long ***tableau_matrix, unsigned int **basic_variables)
{
	unsigned int i, j, number_of_positive_elements = 0, positive_element_row_index = 0;

	(*basic_variables) = malloc(sizeof(unsigned int) * number_of_constraints);
	for (i = 0; i < number_of_constraints; i++)
	{
		(*basic_variables)[i] = 0;
	} // end of i loop

	for (j = 0; j < (number_of_variables + number_of_slack_surplus_variables + number_of_artificial_variables); j++)
	{
		number_of_positive_elements = 0;
		for (i = 0; i < number_of_constraints; i++)
		{
			if (tableau_matrix[0][i][j + 1] > 0)
			{
				number_of_positive_elements++;
				positive_element_row_index = i;
			}
		} // end of i loop

		if (number_of_positive_elements == 1)
		{
			(*basic_variables)[positive_element_row_index] = j + 1;
		}
	} // end of j loop
	
	return;
} // end of CreateBasicVariables function

/* This function creates the basic feasible solution to keep track of the best variable values in the current basis */
void CreateBasicFeasibleSolution(unsigned int number_of_constraints, unsigned int number_of_variables, unsigned int number_of_slack_surplus_variables, unsigned int number_of_artificial_variables, long long ***tableau_matrix, unsigned int *basic_variables, long long ***basic_feasible_solution)
{
	unsigned int i, k;

	(*basic_feasible_solution) = malloc(sizeof(long long*) * 2);
	for (k = 0; k < 2; k++)
	{
		(*basic_feasible_solution)[k] = malloc(sizeof(long long) * (number_of_variables + number_of_slack_surplus_variables + number_of_artificial_variables));
	} // end of k loop

	for (i = 0; i < number_of_variables + number_of_slack_surplus_variables + number_of_artificial_variables; i++)
	{
		(*basic_feasible_solution)[0][i] = 0;
		(*basic_feasible_solution)[1][i] = 1;
	} // end of i loop

	for (i = 0; i < number_of_constraints; i++)
	{
		LongLongRationalDivision(tableau_matrix[0][i][0], tableau_matrix[1][i][0], tableau_matrix[0][i][basic_variables[i]], tableau_matrix[1][i][basic_variables[i]], &(*basic_feasible_solution)[0][basic_variables[i] - 1], &(*basic_feasible_solution)[1][basic_variables[i] - 1]);
	} // end of i loop
	
	return;
} // end of CreateBasicFeasibleSolution function

/* This function writes to disk the initial LP array values */
void WriteInitialLPValues(unsigned int number_of_constraints, unsigned int number_of_variables, unsigned int number_of_slack_surplus_variables, unsigned int number_of_artificial_variables, unsigned int *tableau_current_size, unsigned int *basic_variables, long long **basic_feasible_solution, long long ***tableau_matrix)
{
	unsigned int i, j;

	FILE *outfile_initial_LP_objective_function_optimal_value = fopen("outputs/initial_LP_objective_function_optimal_value.txt", "w"); // write only
	fprintf(outfile_initial_LP_objective_function_optimal_value, "%lf\n", (double)tableau_matrix[0][number_of_constraints][0] / tableau_matrix[1][number_of_constraints][0]);
	fclose(outfile_initial_LP_objective_function_optimal_value);

	FILE *outfile_initial_LP_basic_variables = fopen("outputs/initial_LP_basic_variables.txt", "w"); // write only
	for (i = 0; i < number_of_constraints; i++)
	{
		fprintf(outfile_initial_LP_basic_variables, "%u\n", basic_variables[i]);
	} // end of i loop
	fclose(outfile_initial_LP_basic_variables);

	FILE *outfile_initial_LP_basic_feasible_solution = fopen("outputs/initial_LP_basic_feasible_solution.txt", "w"); // write only
	for (j = 0; j < (number_of_variables + number_of_slack_surplus_variables + number_of_artificial_variables); j++)
	{
		fprintf(outfile_initial_LP_basic_feasible_solution, "%lf\n", (double)basic_feasible_solution[0][j] / basic_feasible_solution[1][j]);
	} // end of j loop
	fclose(outfile_initial_LP_basic_feasible_solution);

	FILE *outfile_initial_LP_tableau_matrix = fopen("outputs/initial_LP_tableau_matrix.txt", "w"); // write only
	for (i = 0; i < tableau_current_size[0]; i++)
	{
		for (j = 0; j < tableau_current_size[1]; j++)
		{
			fprintf(outfile_initial_LP_tableau_matrix, "%lf\t", (double)tableau_matrix[0][i][j] / tableau_matrix[1][i][j]);
		} // end of j loop
		fprintf(outfile_initial_LP_tableau_matrix, "\n");
	} // end of i loop
	fclose(outfile_initial_LP_tableau_matrix);
	
	return;
} // end of WriteInitialLPValues function

/* This function prints the initial constraint, variable, etc. counts */
void PrintInitialCounts(unsigned int number_of_constraints, unsigned int number_of_variables, unsigned int number_of_less_than_or_equal_to_constraints, unsigned int number_of_equal_to_constraints, unsigned int number_of_greater_than_or_equal_to_constraints, unsigned int number_of_slack_surplus_variables, unsigned int number_of_artificial_variables, unsigned int *tableau_current_size)
{
	printf("PrintInitialCounts: number_of_constraints = %u & number_of_variables = %u\n", number_of_constraints, number_of_variables);
	printf("PrintInitialCounts: number_of_less_than_or_equal_to_constraints = %u, number_of_equal_to_constraints = %u, number_of_greater_than_or_equal_to_constraints = %u\n", number_of_less_than_or_equal_to_constraints, number_of_equal_to_constraints, number_of_greater_than_or_equal_to_constraints);
	printf("PrintInitialCounts: number_of_slack_surplus_variables = %u & number_of_artificial_variables = %u\n", number_of_slack_surplus_variables, number_of_artificial_variables);
	printf("PrintInitialCounts: tableau_current_size[0] = %u & tableau_current_size[1] = %u\n", tableau_current_size[0], tableau_current_size[1]);
	
	return;
} // end of PrintInitialCounts function

/* This function writes to disk the final LP array values */
void WriteFinalLPValues(unsigned int number_of_constraints, unsigned int number_of_variables, unsigned int number_of_slack_surplus_variables, unsigned int number_of_artificial_variables, unsigned int *tableau_current_size, unsigned int *basic_variables, long long **basic_feasible_solution, long long ***tableau_matrix)
{
	unsigned int i, j;

	FILE *outfile_final_LP_objective_function_optimal_value = fopen("outputs/final_LP_objective_function_optimal_value.txt", "w"); // write only
	fprintf(outfile_final_LP_objective_function_optimal_value, "%lf\n", (double)tableau_matrix[0][number_of_constraints][0] / tableau_matrix[1][number_of_constraints][0]);
	fclose(outfile_final_LP_objective_function_optimal_value);

	FILE *outfile_final_LP_basic_variables = fopen("outputs/final_LP_basic_variables.txt", "w"); // write only
	for (i = 0; i < number_of_constraints; i++)
	{
		fprintf(outfile_final_LP_basic_variables, "%u\n", basic_variables[i]);
	} // end of i loop
	fclose(outfile_final_LP_basic_variables);

	FILE *outfile_final_LP_basic_feasible_solution = fopen("outputs/final_LP_basic_feasible_solution.txt", "w"); // write only
	for (j = 0; j < (number_of_variables + number_of_slack_surplus_variables + number_of_artificial_variables); j++)
	{
		fprintf(outfile_final_LP_basic_feasible_solution, "%lf\n", (double)basic_feasible_solution[0][j] / basic_feasible_solution[1][j]);
	} // end of j loop
	fclose(outfile_final_LP_basic_feasible_solution);

	FILE *outfile_final_LP_tableau_matrix = fopen("outputs/final_LP_tableau_matrix.txt", "w"); // write only
	for (i = 0; i < tableau_current_size[0]; i++)
	{
		for (j = 0; j < tableau_current_size[1]; j++)
		{
			fprintf(outfile_final_LP_tableau_matrix, "%lf\t", (double)tableau_matrix[0][i][j] / tableau_matrix[1][i][j]);
		} // end of j loop
		fprintf(outfile_final_LP_tableau_matrix, "\n");
	} // end of i loop
	fclose(outfile_final_LP_tableau_matrix);
	
	return;
} // end of WriteFinalLPValues function

/* This function writes to disk the final MILP array values */
void WriteFinalMILPValues(unsigned int number_of_variables, long long *optimal_objective_function_value, long long **optimal_variable_values)
{
	unsigned int i;

	FILE *outfile_final_MILP_optimal_objective_function_value = fopen("outputs/final_MILP_optimal_objective_function_value.txt", "w"); // write only
	fprintf(outfile_final_MILP_optimal_objective_function_value, "%lf\n", (double)optimal_objective_function_value[0] / optimal_objective_function_value[1]);
	fclose(outfile_final_MILP_optimal_objective_function_value);

	FILE *outfile_final_MILP_optimal_variable_values = fopen("outputs/final_MILP_optimal_variable_values.txt", "w"); // write only
	for (i = 0; i < number_of_variables; i++)
	{
		fprintf(outfile_final_MILP_optimal_variable_values, "%lf\n", (double)optimal_variable_values[0][i] / optimal_variable_values[1][i]);
	} // end of i loop
	fclose(outfile_final_MILP_optimal_variable_values);
	
	return;
} // end of WriteFinalMILPValues function

/* This function updates the optimal objective function and variable solution */
void UpdateOptimalObjectiveFunctionAndVariableSolution(unsigned int number_of_variables, long long **basic_feasible_solution, long long *objective_function_initial_constant, long long **objective_function_coefficient_vector, long long *optimal_objective_function_value, long long **optimal_variable_values)
{
	unsigned int i;

	optimal_objective_function_value[0] = objective_function_initial_constant[0];
	optimal_objective_function_value[1] = objective_function_initial_constant[1];

	for (i = 0; i < number_of_variables; i++)
	{
		LongLongRationalAddition(optimal_objective_function_value[0], optimal_objective_function_value[1], objective_function_coefficient_vector[0][i] * basic_feasible_solution[0][i], objective_function_coefficient_vector[1][i] * basic_feasible_solution[1][i], &optimal_objective_function_value[0], &optimal_objective_function_value[1]);

		optimal_variable_values[0][i] = basic_feasible_solution[0][i];
		optimal_variable_values[1][i] = basic_feasible_solution[1][i];
	} // end of i loop
	
	return;
} // end of UpdateOptimalObjectiveFunctionAndVariableSolution function

/*********************************************************************************/
/*********************************** SIMPLEX *************************************/
/*********************************************************************************/

/* This function finds the optimal solution for the given variables and constraints */
int SimplexAlgorithm(unsigned int number_of_constraints, unsigned int number_of_variables, unsigned int number_of_slack_surplus_variables, unsigned int *number_of_artificial_variables, unsigned int *basic_variables, long long **basic_feasible_solution, long long ***tableau_matrix, unsigned int *tableau_current_size)
{
	unsigned int i, k;
	int error_code = 0;

	/*********************************************************************************/
	/*********************************** PHASE 1 *************************************/
	/*********************************************************************************/

	if ((*number_of_artificial_variables) > 0)
	{
		/* This function transforms the tableau by removing artificial variables to obtain a basic feasible solution */
		error_code = SimplexPhase1(number_of_constraints, number_of_variables, number_of_slack_surplus_variables, number_of_artificial_variables, basic_variables, tableau_matrix, tableau_current_size);
		tableau_current_size[0]--; // decrement current rows since we've eliminated artificial objective function row
	}

	/*********************************************************************************/
	/*********************************** PHASE 2 *************************************/
	/*********************************************************************************/

	if (error_code == 0)
	{
		/* This function starts from an basic feasible solution and iterates toward the optimal solution */
		error_code = SimplexPhase2(number_of_constraints, number_of_variables, number_of_slack_surplus_variables, basic_variables, tableau_matrix, tableau_current_size);
	}

	/*********************************************************************************/
	/*********************************** SOLUTION *************************************/
	/*********************************************************************************/

	if (error_code == 0)
	{
		UpdateBasicFeasibleSolution(number_of_variables + number_of_slack_surplus_variables, number_of_constraints, basic_variables, basic_feasible_solution, tableau_matrix);
	}

	return error_code;
} // end of SimplexAlgorithm function

/* This function transforms the tableau by removing artificial variables to obtain a basic feasible solution */
int SimplexPhase1(unsigned int number_of_constraints, unsigned int number_of_variables, unsigned int number_of_slack_surplus_variables, unsigned int *number_of_artificial_variables, unsigned int *basic_variables, long long ***tableau_matrix, unsigned int *tableau_current_size)
{
	unsigned int i, j, k, l;
	int error_code = 0;

	/* PHASE I (first find the feasible region) */

	/* Remove artificial variables from objective function */
	unsigned int artificial_variable_col_index = 0;
	for (k = 0; k < (*number_of_artificial_variables); k++) // go through each artificial variable
	{
		artificial_variable_col_index = number_of_variables + number_of_slack_surplus_variables + 1 + k;

		for (i = 0; i < (tableau_current_size[0] - 2); i++) // try rows looking for kth artificial variable row
		{
			if (tableau_matrix[0][i][artificial_variable_col_index] == 1 && tableau_matrix[1][i][artificial_variable_col_index] == 1) // if this is the kth artificial variable's row
			{
				for (j = 0; j < tableau_current_size[1]; j++)
				{
					LongLongRationalAddition(tableau_matrix[0][number_of_constraints + 1][j], tableau_matrix[1][number_of_constraints + 1][j], tableau_matrix[0][i][j], tableau_matrix[1][i][j], &tableau_matrix[0][number_of_constraints + 1][j], &tableau_matrix[1][number_of_constraints + 1][j]);
				} // end of j loop

				break; // breaks i loop since we already found kth artificial variable's row
			} // end of if this is the kth artificial variable's row
		} // end of i loop
	} // end of k loop

	/* Count initial number of positive elements in the objective function row */
	unsigned int number_of_positive_objective_function_elements = 0;
	for (j = 0; j < (tableau_current_size[1] - 1); j++)
	{
		if (tableau_matrix[0][tableau_current_size[0] - 1][j + 1] > 0)
		{
			number_of_positive_objective_function_elements++;
		}
	} // end of j loop
	
	unsigned int most_positive_objective_function_element_index = 0, smallest_non_negative_ratio_index = 0, smallest_artificial_variable_non_negative_ratio_index = 0, pivot_col_index = 0, pivot_row_index = 0;
	double most_positive_objective_function_element_value = 0, b_a_ratio_double = 0;

	long long pivot_value[2];
	pivot_value[0] = 0;
	pivot_value[1] = 1;
	
	long long b_a_ratio[2];
	b_a_ratio[0] = 0;
	b_a_ratio[1] = 1;
	
	long long smallest_non_negative_ratio_value[2];
	smallest_non_negative_ratio_value[0] = LLONG_MAX;
	smallest_non_negative_ratio_value[1] = 1;
	
	long long smallest_artificial_variable_non_negative_ratio_value[2];
	smallest_artificial_variable_non_negative_ratio_value[0] = LLONG_MAX;
	smallest_artificial_variable_non_negative_ratio_value[1] = 1;

	unsigned int iteration = 0;
	while (number_of_positive_objective_function_elements > 0 && (*number_of_artificial_variables) > 0 && error_code == 0)
	{
		/* Find most positive objective function element amongst the variables which will become our entering variable */
		most_positive_objective_function_element_value = 0;
		most_positive_objective_function_element_index = 0;

		for (j = 0; j < (tableau_current_size[1] - 1); j++)
		{
			if ((double)tableau_matrix[0][tableau_current_size[0] - 1][j + 1] / tableau_matrix[1][tableau_current_size[0] - 1][j + 1] > most_positive_objective_function_element_value)
			{
				most_positive_objective_function_element_value = (double)tableau_matrix[0][tableau_current_size[0] - 1][j + 1] / tableau_matrix[1][tableau_current_size[0] - 1][j + 1];
				most_positive_objective_function_element_index = j + 1;
			}
		} // end of j loop

		if (most_positive_objective_function_element_value > 0) // if a pivot column was found
		{
			pivot_col_index = most_positive_objective_function_element_index;

			/* Search for smallest non negative ratio of bi / aij which will be departing variable */
			smallest_non_negative_ratio_index = 0;
			smallest_non_negative_ratio_value[0] = LLONG_MAX;
			smallest_non_negative_ratio_value[1] = 1;

			smallest_artificial_variable_non_negative_ratio_index = 0;
			smallest_artificial_variable_non_negative_ratio_value[0] = LLONG_MAX;
			smallest_artificial_variable_non_negative_ratio_value[1] = 1;

			for (i = 0; i < tableau_current_size[0] - 2; i++) // go through just constraint rows
			{
				if (tableau_matrix[0][i][pivot_col_index] > 0) // pivot element needs to be positive
				{
					LongLongRationalDivision(tableau_matrix[0][i][0], tableau_matrix[1][i][0], tableau_matrix[0][i][pivot_col_index], tableau_matrix[1][i][pivot_col_index], &b_a_ratio[0], &b_a_ratio[1]);
					
					b_a_ratio_double = (double)b_a_ratio[0] / b_a_ratio[1];
					
					if (b_a_ratio_double < (double)smallest_non_negative_ratio_value[0] / smallest_non_negative_ratio_value[1] && b_a_ratio_double >= 0) // this should already be using Bland's Rule since it is taking the lowest index in ties to avoid cycles
					{
						smallest_non_negative_ratio_value[0] = b_a_ratio[0];
						smallest_non_negative_ratio_value[1] = b_a_ratio[1];
						smallest_non_negative_ratio_index = i;
					}

					if (basic_variables[i] >= number_of_variables + number_of_slack_surplus_variables)
					{
						if (b_a_ratio_double < (double)smallest_artificial_variable_non_negative_ratio_value[0] / smallest_artificial_variable_non_negative_ratio_value[1] && b_a_ratio_double >= 0) // this should already be using Bland's Rule since it is taking the lowest index in ties to avoid cycles
						{
							smallest_artificial_variable_non_negative_ratio_value[0] = b_a_ratio[0];
							smallest_artificial_variable_non_negative_ratio_value[1] = b_a_ratio[1];
							smallest_artificial_variable_non_negative_ratio_index = i;
						}
					}
				}
			} // end of i loop
			
			if (!(smallest_non_negative_ratio_value[0] == LLONG_MAX && smallest_non_negative_ratio_value[1] == 1)) // if pivot row was found
			{
				/* IF LOWEST RATIO OCCURS BOTH IN AN ARTIFICIAL VARIABLE ROW AND A NON-NEGATIVE VARIABLE ROW THEN MUST CHOOSE PIVOT ROW AS NON-NEGATIVE VARIABLE ROW */
				if (smallest_artificial_variable_non_negative_ratio_index == 0) // if there are NO more artificial variables
				{
					pivot_row_index = smallest_non_negative_ratio_index;
				} // end of if there are NO more artificial variables
				else // if there are still artificial variables
				{
					if (smallest_non_negative_ratio_index != smallest_artificial_variable_non_negative_ratio_index)
					{
						if (smallest_non_negative_ratio_value[0] == smallest_artificial_variable_non_negative_ratio_value[0] && smallest_non_negative_ratio_value[1] == smallest_artificial_variable_non_negative_ratio_value[1])
						{
							pivot_row_index = smallest_artificial_variable_non_negative_ratio_index;
						}
						else
						{
							pivot_row_index = smallest_non_negative_ratio_index;
						}
					}
					else
					{
						pivot_row_index = smallest_non_negative_ratio_index;
					}
				} // end of if there are still artificial variables

				pivot_value[0] = tableau_matrix[0][pivot_row_index][pivot_col_index];
				pivot_value[1] = tableau_matrix[1][pivot_row_index][pivot_col_index];

				if (basic_variables[pivot_row_index] >= number_of_variables + number_of_slack_surplus_variables) // if basic variable is an artificial variable
				{
					(*number_of_artificial_variables)--;
					tableau_current_size[1]--;
				}

				/* Remove departing variable from and add entering variable to basic variables */
				basic_variables[pivot_row_index] = pivot_col_index;

				/* This function performs Gauss-Jordan Elimination on the pivot column */
				PivotColumnGaussJordanElimnation(tableau_current_size[0], tableau_current_size[1], pivot_row_index, pivot_col_index, pivot_value, tableau_matrix);

				/* Count again the number of variables that are positive still */
				number_of_positive_objective_function_elements = 0;
				for (j = 0; j < (tableau_current_size[1] - 1); j++)
				{
					if (tableau_matrix[0][tableau_current_size[0] - 1][j + 1] > 0)
					{
						number_of_positive_objective_function_elements++;
					}
				} // end of j loop

//				printf("SimplexPhase1: BREAKPOINT! iteration = %u, number_of_positive_objective_function_elements = %u & (*number_of_artificial_variables) = %u\n", iteration, number_of_positive_objective_function_elements, (*number_of_artificial_variables));
				iteration++;
			} // end of if pivot row was found
			else
			{
				error_code = 1; // infeasible
			}
		} // end of if a pivot column was found
		else
		{
			error_code = 1; // infeasible
		}
	} // end of while (number_of_positive_objective_function_elements > 0)

	if ((*number_of_artificial_variables) > 0)
	{
		printf("SimplexPhase1: ERROR: No basic feasible solution!\n");
		error_code = 1;
	}

	return error_code;
} // end of SimplexPhase1 function

/* This function starts from an basic feasible solution and iterates toward the optimal solution */
int SimplexPhase2(unsigned int number_of_constraints, unsigned int number_of_variables, unsigned int number_of_slack_surplus_variables, unsigned int *basic_variables, long long ***tableau_matrix, unsigned int *tableau_current_size)
{
	int i, j, k, l, error_code = 0;
	double old_optimum = -DBL_MAX;

	/* PHASE II (find the optimal solution in the feasible region we found above) */

	/* Count initial number of negative elements in the objective function row */
	unsigned int number_of_negative_objective_function_elements = 0;
	for (j = 0; j < (tableau_current_size[1] - 1); j++)
	{
		if (tableau_matrix[0][tableau_current_size[0] - 1][j + 1] < 0)
		{
			number_of_negative_objective_function_elements++;
		}
	} // end of j loop

	unsigned int most_negative_objective_function_element_index = 0, number_of_positive_elements = 0, smallest_non_negative_ratio_index = 0, pivot_col_index = 0, pivot_row_index = 0;
	double most_negative_objective_function_element_value, b_a_ratio = 0, smallest_non_negative_ratio_value = DBL_MAX;

	long long pivot_value[2];
	pivot_value[0] = 0;
	pivot_value[1] = 1;

	while (number_of_negative_objective_function_elements > 0 && error_code == 0) // while there is still a optimal solution
	{
		/* Search for most negative element in objective function row */
		most_negative_objective_function_element_index = 0;
		most_negative_objective_function_element_value = 0;

		for (j = 0; j < (tableau_current_size[1] - 1); j++)
		{
			if (tableau_matrix[0][tableau_current_size[0] - 1][j + 1] < 0) // if jth objective function coefficient is negative
			{
				if ((double)tableau_matrix[0][tableau_current_size[0] - 1][j + 1] / tableau_matrix[1][tableau_current_size[0] - 1][j + 1] < most_negative_objective_function_element_value) // this should already be using Bland's Rule since it is taking the lowest index in ties to avoid cycles
				{
					most_negative_objective_function_element_value = (double)tableau_matrix[0][tableau_current_size[0] - 1][j + 1] / tableau_matrix[1][tableau_current_size[0] - 1][j + 1];
					most_negative_objective_function_element_index = j + 1;
				}
			} // end of if jth objective function coefficient is negative
		} // end of j loop

		if (most_negative_objective_function_element_value < 0) // if pivot column was found
		{
			pivot_col_index = most_negative_objective_function_element_index;

			/* Check to see if the problem is unbounded */
			number_of_positive_elements = 0;
			for (i = 0; i < (tableau_current_size[0] - 1); i++)
			{
				if (tableau_matrix[0][i][pivot_col_index] > 0)
				{
					number_of_positive_elements++;
				}
			} // end of i loop

			if (number_of_positive_elements > 0) // if the problem is bounded still
			{
				/* Search for smallest non negative ratio of bi / aij */
				smallest_non_negative_ratio_index = 0;
				smallest_non_negative_ratio_value = DBL_MAX;

				for (i = 0; i < (tableau_current_size[0] - 1); i++)
				{
					if (tableau_matrix[0][i][pivot_col_index] > 0) // pivot element needs to be positive
					{
						b_a_ratio = ((double)tableau_matrix[0][i][0] / tableau_matrix[1][i][0]) / ((double)tableau_matrix[0][i][pivot_col_index] / tableau_matrix[1][i][pivot_col_index]);

						if (b_a_ratio < smallest_non_negative_ratio_value && b_a_ratio >= 0) // this should already be using Bland's Rule since it is taking the lowest index in ties to avoid cycles
						{
							smallest_non_negative_ratio_value = b_a_ratio;
							smallest_non_negative_ratio_index = i;
						}
					}
				} // end of i loop

				if (smallest_non_negative_ratio_value < DBL_MAX) // if pivot row was found
				{
					pivot_row_index = smallest_non_negative_ratio_index;

					pivot_value[0] = tableau_matrix[0][pivot_row_index][pivot_col_index];
					pivot_value[1] = tableau_matrix[1][pivot_row_index][pivot_col_index];

					/* Remove departing variable from and add entering variable to basic variables */
					basic_variables[pivot_row_index] = pivot_col_index;

					/* This function performs Gauss-Jordan Elimination on the pivot column */
					PivotColumnGaussJordanElimnation(tableau_current_size[0], tableau_current_size[1], pivot_row_index, pivot_col_index, pivot_value, tableau_matrix);

					/* Count new number of negative elements in the objective function row */
					number_of_negative_objective_function_elements = 0;
					for (j = 0; j < (tableau_current_size[1] - 1); j++)
					{
						if (tableau_matrix[0][tableau_current_size[0] - 1][j + 1] < 0)
						{
							number_of_negative_objective_function_elements++;
						}
					} // end of j loop

					if (old_optimum == (double)tableau_matrix[0][tableau_current_size[0] - 1][0] / tableau_matrix[1][tableau_current_size[0] - 1][0])
					{
						printf("SimplexPhase2: NOTE: Problem is Degenerate!\n");
					}
					else
					{
						old_optimum = (double)tableau_matrix[0][tableau_current_size[0] - 1][0] / tableau_matrix[1][tableau_current_size[0] - 1][0];
					}
				} // end of if pivot row was found
				else
				{
					error_code = 1; // infeasible
				}
			} // end of if the problem is bounded still
			else
			{
				printf("SimplexPhase2: ERROR: Problem is Unbounded!\n");
				error_code = 2; // unbounded
				break; // break while loop
			}
		} // end of if pivot column was found
		else
		{
			error_code = 1; // infeasible
		}
	} // end of while there is still a optimal solution

	return error_code;
} // end of SimplexPhase2 function

/* This function performs Gauss-Jordan Elimination on the pivot column */
void PivotColumnGaussJordanElimnation(unsigned int number_of_rows, unsigned int number_of_columns, unsigned int pivot_row_index, unsigned int pivot_col_index, long long *pivot_value, long long ***tableau_matrix)
{
	unsigned int i, j, k;

	/* Perform Gauss-Jordan Elimination on pivot column around pivot row */

	/* First take care of pivot row */
	if (tableau_matrix[0][pivot_row_index][pivot_col_index] != 1 || tableau_matrix[1][pivot_row_index][pivot_col_index] != 1)
	{
		for (j = 0; j < number_of_columns; j++)
		{
			if (j == pivot_col_index)
			{
				/* Set to 1 */
				tableau_matrix[0][pivot_row_index][j] = 1;
				tableau_matrix[1][pivot_row_index][j] = 1;
			}
			else
			{
				LongLongRationalDivision(tableau_matrix[0][pivot_row_index][j], tableau_matrix[1][pivot_row_index][j], pivot_value[0], pivot_value[1], &tableau_matrix[0][pivot_row_index][j], &tableau_matrix[1][pivot_row_index][j]);
			}
		} // end of j loop
	}

	/* Now take care of other rows */
	for (i = 0; i < number_of_rows; i++)
	{
		if (i != pivot_row_index) // if ith row is NOT pivot row, we already took care of the pivot row
		{
			if (tableau_matrix[0][i][pivot_col_index] != 0) // if the element in the ith row of the pivot column is NOT already 0
			{
				pivot_value[0] = tableau_matrix[0][i][pivot_col_index];
				pivot_value[1] = tableau_matrix[1][i][pivot_col_index];

				for (j = 0; j < number_of_columns; j++)
				{
					if (j == pivot_col_index)
					{
						/* Set to 0 */
						tableau_matrix[0][i][j] = 0;
						tableau_matrix[1][i][j] = 1;
					}
					else
					{
						LongLongRationalAddition(tableau_matrix[0][i][j], tableau_matrix[1][i][j], -pivot_value[0] * tableau_matrix[0][pivot_row_index][j], pivot_value[1] * tableau_matrix[1][pivot_row_index][j], &tableau_matrix[0][i][j], &tableau_matrix[1][i][j]);
					}
				} // end of j loop
			} // end of if the element in the ith row of the pivot column is NOT already 0
		} // end of if ith row is NOT pivot row, we already took care of the pivot row
	} // end of i loop
	
	return;
} // end of PivotColumnGaussJordanElimnation function

/* This function updates the basic feasible solution */
void UpdateBasicFeasibleSolution(unsigned int number_of_total_variables, unsigned int number_of_constraints, unsigned int *basic_variables, long long **basic_feasible_solution, long long ***tableau_matrix)
{
	unsigned int i;

	/* Update Basic Feasible Solution */
	for (i = 0; i < number_of_total_variables; i++)
	{
		basic_feasible_solution[0][i] = 0;
		basic_feasible_solution[1][i] = 1;
	} // end of i loop

	for (i = 0; i < number_of_constraints; i++)
	{
		LongLongRationalDivision(tableau_matrix[0][i][0], tableau_matrix[1][i][0], tableau_matrix[0][i][basic_variables[i]], tableau_matrix[1][i][basic_variables[i]], &basic_feasible_solution[0][basic_variables[i] - 1], &basic_feasible_solution[1][basic_variables[i] - 1]);
	} // end of i loop
	
	return;
} // end of updateBasicFeasibleSolution function

/*********************************************************************************/
/******************************* BRANCH AND BOUND ********************************/
/*********************************************************************************/

/* This function initiates the branch and bound */
int BranchAndBoundMILP(int maximization_problem, unsigned int *number_of_constraints, unsigned int number_of_variables, unsigned int *number_of_slack_surplus_variables, unsigned int *number_of_artificial_variables, unsigned int **basic_variables, long long ***basic_feasible_solution, long long ****tableau_matrix, unsigned int *tableau_current_size, unsigned int *tableau_max_size, unsigned int number_of_variables_required_to_be_integer, unsigned int number_of_variables_required_to_be_binary, int *variable_special_requirements, long long *optimal_objective_function_value, long long **optimal_variable_values, long long **objective_function_coefficient_vector, int lp_error_code)
{
	unsigned int i, j, k, number_of_variables_needing_to_be_integer_or_binary_currently_integer_or_binary = 0, last_variable_that_still_needs_to_become_integer_or_binary_index = 0;
	int error_code = 0;

	CountNumberOfVariablesNeedingToBeIntegerOrBinaryThatAlreadyAre(number_of_variables, number_of_variables_required_to_be_integer, number_of_variables_required_to_be_binary, variable_special_requirements, (*basic_feasible_solution), &number_of_variables_needing_to_be_integer_or_binary_currently_integer_or_binary, &last_variable_that_still_needs_to_become_integer_or_binary_index);

	if (number_of_variables_needing_to_be_integer_or_binary_currently_integer_or_binary < number_of_variables_required_to_be_integer + number_of_variables_required_to_be_binary)
	{
		/* Create variables to track the best optimal values and the corresponding variables */
		int best_milp_error_code = 1;

		double best_mixed_integer_optimal_value_double;

		long long best_mixed_integer_optimal_value[2];

		if (maximization_problem == 1)
		{
			best_mixed_integer_optimal_value_double = -DBL_MAX; // negative infinity
			best_mixed_integer_optimal_value[0] = -LLONG_MAX;
		}
		else
		{
			best_mixed_integer_optimal_value_double = DBL_MAX; // positive infinity
			best_mixed_integer_optimal_value[0] = LLONG_MAX;
		}
		best_mixed_integer_optimal_value[1] = 1;

		long long **best_mixed_integer_variable_values;
		best_mixed_integer_variable_values = malloc(sizeof(long long*) * 2);
		for (k = 0; k < 2; k++)
		{
			best_mixed_integer_variable_values[k] = malloc(sizeof(long long) * number_of_variables);
		} // end of k loop

		for (i = 0; i < number_of_variables; i++)
		{
			best_mixed_integer_variable_values[0][i] = 0;
			best_mixed_integer_variable_values[1][i] = 1;
		} // end of i loop

		/* Call recursive branch and bound function */
		error_code = BranchAndBoundMILPRecursive(maximization_problem, number_of_constraints, number_of_variables, number_of_slack_surplus_variables, number_of_artificial_variables, basic_variables, basic_feasible_solution, tableau_matrix, tableau_current_size, tableau_max_size, number_of_variables_required_to_be_integer, number_of_variables_required_to_be_binary, variable_special_requirements, last_variable_that_still_needs_to_become_integer_or_binary_index, &best_mixed_integer_optimal_value_double, best_mixed_integer_optimal_value, best_mixed_integer_variable_values, &best_milp_error_code, 0);

		error_code = best_milp_error_code;

		if (best_milp_error_code == 0) // feasible
		{
			printf("\nBranchAndBoundMILP: Optimal MILP!\n");
			for (i = 0; i < number_of_variables; i++)
			{
				(*basic_feasible_solution)[0][i] = best_mixed_integer_variable_values[0][i];
				(*basic_feasible_solution)[1][i] = best_mixed_integer_variable_values[1][i];

				/* Update overall optimal solution */
				optimal_variable_values[0][i] = (*basic_feasible_solution)[0][i];
				optimal_variable_values[1][i] = (*basic_feasible_solution)[1][i];
			} // end of i loop
			optimal_objective_function_value[0] = best_mixed_integer_optimal_value[0];
			optimal_objective_function_value[1] = best_mixed_integer_optimal_value[1];
		}
		else if (best_milp_error_code == 1) // infeasible
		{
			printf("\nBranchAndBoundMILP: Infeasible MILP!\n");
		}
		else if (best_milp_error_code == 2) // unbounded
		{
			printf("\nBranchAndBoundMILP: Unbounded MILP!\n");
		}

		/* Free dynamic memory */
		for (k = 0; k < 2; k++)
		{
			free(best_mixed_integer_variable_values[k]);
		} // end of k loop
		free(best_mixed_integer_variable_values);
	}
	else // if all variables that need to be integer or binary are already
	{
		if (lp_error_code == 0)
		{
			error_code = 0;

			printf("\nBranchAndBoundMILP: If there were any integer or binary constraints, they were already satisfied before MILP B&B by optimal LP!\n");
		}
		else if (lp_error_code == 2)
		{
			error_code = 2;

			printf("\nBranchAndBoundMILP: LP was unbounded and since computers can't represent irrational numbers, then MILP is unbounded!\n");
		}
		else
		{
			error_code = 999;
			printf("\nBranchAndBoundMILP: HOW DID I GET IN MILP B&B?! ERROR!\n");
		}
	}

	return error_code;
} // end of BranchAndBoundMILP function

/* This function counts the number of variables that need to be integer or binary AND already are */
void CountNumberOfVariablesNeedingToBeIntegerOrBinaryThatAlreadyAre(unsigned int number_of_variables, unsigned int number_of_variables_required_to_be_integer, unsigned int number_of_variables_required_to_be_binary, int *variable_special_requirements, long long **basic_feasible_solution, unsigned int *number_of_variables_needing_to_be_integer_or_binary_currently_integer_or_binary, unsigned int *last_variable_that_still_needs_to_become_integer_or_binary_index)
{
	unsigned int  i;

	for (i = 0; i < number_of_variables; i++)
	{
		if (number_of_variables_required_to_be_integer + number_of_variables_required_to_be_binary > 0)
		{
			if (variable_special_requirements[i] == 1) // if ith variable is required to be integer
			{
				if (basic_feasible_solution[1][i] == 1) // if denominator is 1, then already integer
				{
					(*number_of_variables_needing_to_be_integer_or_binary_currently_integer_or_binary)++;
				}
				else
				{
					(*last_variable_that_still_needs_to_become_integer_or_binary_index) = i;
				}
			}
			else if (variable_special_requirements[i] == 2) // if ith variable is required to be binary
			{
				if (basic_feasible_solution[0][i] == 0) // if numerator is 0, then already 0
				{
					(*number_of_variables_needing_to_be_integer_or_binary_currently_integer_or_binary)++;
				}
				else if (basic_feasible_solution[0][i] == basic_feasible_solution[1][i]) // if numerator = denominator, then already 1
				{
					(*number_of_variables_needing_to_be_integer_or_binary_currently_integer_or_binary)++;
				}
				else
				{
					if ((double)basic_feasible_solution[0][i] / basic_feasible_solution[1][i] >= 0 && (double)basic_feasible_solution[0][i] / basic_feasible_solution[1][i] <= 1)
					{
						(*last_variable_that_still_needs_to_become_integer_or_binary_index) = i;
					}
					else
					{
						printf("CountNumberOfVariablesNeedingToBeIntegerOrBinaryThatAlreadyAre: INFEASIBLE SOLUTION MADE IT THIS FAR INTO <= MILP B&B!  ERROR!\n");
						(*last_variable_that_still_needs_to_become_integer_or_binary_index) = i;
					}
				}
			}
		}
	} // end of i loop
	
	return;
} // end of CountNumberOfVariablesNeedingToBeIntegerOrBinaryThatAlreadyAre function

/* This function recursively applies branch and bound */
int BranchAndBoundMILPRecursive(int maximization_problem, unsigned int *number_of_constraints, unsigned int number_of_variables, unsigned int *number_of_slack_surplus_variables, unsigned int *number_of_artificial_variables, unsigned int **basic_variables, long long ***basic_feasible_solution, long long ****tableau_matrix, unsigned int *tableau_current_size, unsigned int *tableau_max_size, unsigned int number_of_variables_required_to_be_integer, unsigned int number_of_variables_required_to_be_binary, int *variable_special_requirements, unsigned int last_variable_that_still_needs_to_become_integer_or_binary_index, double *best_mixed_integer_optimal_value_double, long long *best_mixed_integer_optimal_value, long long **best_mixed_integer_variable_values, int *best_milp_error_code, unsigned int recursion_level)
{
	unsigned int i, j, k;
	int error_code = 0;

	struct BranchAndBoundState branch_and_bound_state_old;
	
	CreateCacheCountsAndArraysForBranchAndBoundMILPRecursion(tableau_current_size, (*number_of_constraints), number_of_variables, (*number_of_slack_surplus_variables), (*number_of_artificial_variables), (*basic_variables), (*basic_feasible_solution), (*tableau_matrix), &branch_and_bound_state_old);

	/* This will be the value that we try next for the corresponding variable */
	double last_variable_that_still_needs_to_become_integer_value = (double)(*basic_feasible_solution)[0][last_variable_that_still_needs_to_become_integer_or_binary_index] / (*basic_feasible_solution)[1][last_variable_that_still_needs_to_become_integer_or_binary_index];

	/*********************************************************************************/
	/*********************************** LESS THAN ***********************************/
	/*********************************************************************************/

	BranchAndBoundMILPRecursiveLessOrGreaterThan(maximization_problem, number_of_constraints, number_of_variables, number_of_slack_surplus_variables, number_of_artificial_variables, basic_variables, basic_feasible_solution, tableau_matrix, tableau_current_size, tableau_max_size, number_of_variables_required_to_be_integer, number_of_variables_required_to_be_binary, variable_special_requirements, last_variable_that_still_needs_to_become_integer_or_binary_index, last_variable_that_still_needs_to_become_integer_value, best_mixed_integer_optimal_value_double, best_mixed_integer_optimal_value, best_mixed_integer_variable_values, best_milp_error_code, branch_and_bound_state_old, 1, recursion_level);

	/*********************************************************************************/
	/********************************* GREATER THAN **********************************/
	/*********************************************************************************/

	BranchAndBoundMILPRecursiveLessOrGreaterThan(maximization_problem, number_of_constraints, number_of_variables, number_of_slack_surplus_variables, number_of_artificial_variables, basic_variables, basic_feasible_solution, tableau_matrix, tableau_current_size, tableau_max_size, number_of_variables_required_to_be_integer, number_of_variables_required_to_be_binary, variable_special_requirements, last_variable_that_still_needs_to_become_integer_or_binary_index, last_variable_that_still_needs_to_become_integer_value, best_mixed_integer_optimal_value_double, best_mixed_integer_optimal_value, best_mixed_integer_variable_values, best_milp_error_code, branch_and_bound_state_old, -1, recursion_level);

	/* Free dynamic memory */
	for (k = 0; k < 2; k++)
	{
		for (i = 0; i < branch_and_bound_state_old.tableau_current_size[0]; i++)
		{
			free(branch_and_bound_state_old.tableau_matrix[k][i]);
		} // end of i loop
		free(branch_and_bound_state_old.tableau_matrix[k]);
		free(branch_and_bound_state_old.basic_feasible_solution[k]);
	} // end of k loop
	free(branch_and_bound_state_old.tableau_matrix);
	free(branch_and_bound_state_old.basic_feasible_solution);
	free(branch_and_bound_state_old.basic_variables);
	
	return error_code;
} // end of BranchAndBoundMILPRecursive function

/* This function creates cache counts and arrays to reset the active arrays at the end of each recursion level */
void CreateCacheCountsAndArraysForBranchAndBoundMILPRecursion(unsigned int tableau_current_size[2], unsigned int number_of_constraints, unsigned int number_of_variables, unsigned int number_of_slack_surplus_variables, unsigned int number_of_artificial_variables, unsigned int *basic_variables, long long **basic_feasible_solution, long long ***tableau_matrix, struct BranchAndBoundState *branch_and_bound_state_old)
{
	unsigned int i, j, k;

	branch_and_bound_state_old->tableau_current_size[0] = tableau_current_size[0];
	branch_and_bound_state_old->tableau_current_size[1] = tableau_current_size[1];
	branch_and_bound_state_old->number_of_constraints = number_of_constraints;
	branch_and_bound_state_old->number_of_slack_surplus_variables = number_of_slack_surplus_variables;
	branch_and_bound_state_old->number_of_artificial_variables = number_of_artificial_variables;
	
	branch_and_bound_state_old->basic_variables = malloc(sizeof(unsigned int) * branch_and_bound_state_old->number_of_constraints);
	for (i = 0; i < branch_and_bound_state_old->number_of_constraints; i++)
	{
		branch_and_bound_state_old->basic_variables[i] = basic_variables[i];
	} // end of i loop

	branch_and_bound_state_old->basic_feasible_solution = malloc(sizeof(long long*) * 2);
	for (k = 0; k < 2; k++)
	{
		branch_and_bound_state_old->basic_feasible_solution[k] = malloc(sizeof(long long) * (number_of_variables + branch_and_bound_state_old->number_of_slack_surplus_variables + branch_and_bound_state_old->number_of_artificial_variables));
		for (i = 0; i < (number_of_variables + branch_and_bound_state_old->number_of_slack_surplus_variables + branch_and_bound_state_old->number_of_artificial_variables); i++)
		{
			branch_and_bound_state_old->basic_feasible_solution[k][i] = basic_feasible_solution[k][i];
		} // end of i loop
	} // end of k loop

	branch_and_bound_state_old->tableau_matrix = malloc(sizeof(long long*) * 2);
	for (k = 0; k < 2; k++)
	{
		branch_and_bound_state_old->tableau_matrix[k] = malloc(sizeof(long long*) * branch_and_bound_state_old->tableau_current_size[0]);
		for (i = 0; i < branch_and_bound_state_old->tableau_current_size[0]; i++)
		{
			branch_and_bound_state_old->tableau_matrix[k][i] = malloc(sizeof(long long) * branch_and_bound_state_old->tableau_current_size[1]);
			for (j = 0; j < branch_and_bound_state_old->tableau_current_size[1]; j++)
			{
				branch_and_bound_state_old->tableau_matrix[k][i][j] = tableau_matrix[k][i][j];
			} // end of j loop
		} // end of i loop
	} // end of k loop
	
	return;
} // end of CreateCacheCountsAndArraysForBranchAndBoundMILPRecursion function

/* This function handles the less than or greater than branches of the recursive tree of branch and bound */
int BranchAndBoundMILPRecursiveLessOrGreaterThan(int maximization_problem, unsigned int *number_of_constraints, unsigned int number_of_variables, unsigned int *number_of_slack_surplus_variables, unsigned int *number_of_artificial_variables, unsigned int **basic_variables, long long ***basic_feasible_solution, long long ****tableau_matrix, unsigned int *tableau_current_size, unsigned int *tableau_max_size, unsigned int number_of_variables_required_to_be_integer, unsigned int number_of_variables_required_to_be_binary, int *variable_special_requirements, unsigned int last_variable_that_still_needs_to_become_integer_or_binary_index, double last_variable_that_still_needs_to_become_integer_value, double *best_mixed_integer_optimal_value_double, long long *best_mixed_integer_optimal_value, long long **best_mixed_integer_variable_values, int *best_milp_error_code, struct BranchAndBoundState branch_and_bound_state_old, int new_constraint_inequality_direction, unsigned int recursion_level)
{
	int error_code;

	if (new_constraint_inequality_direction == 1) // less than
	{
		/* Last_variable_that_still_needs_to_become_integer_or_binary_index <= new_constraint_constant */
		error_code = AddConstraint(number_of_constraints, number_of_variables, number_of_slack_surplus_variables, number_of_artificial_variables, basic_variables, basic_feasible_solution, tableau_matrix, tableau_current_size, tableau_max_size, last_variable_that_still_needs_to_become_integer_or_binary_index, new_constraint_inequality_direction, (long long)floor(last_variable_that_still_needs_to_become_integer_value), recursion_level);
	}
	else // greater than
	{
		/* Last_variable_that_still_needs_to_become_integer_or_binary_index >= new_constraint_constant */
		error_code = AddConstraint(number_of_constraints, number_of_variables, number_of_slack_surplus_variables, number_of_artificial_variables, basic_variables, basic_feasible_solution, tableau_matrix, tableau_current_size, tableau_max_size, last_variable_that_still_needs_to_become_integer_or_binary_index, new_constraint_inequality_direction, (long long)ceil(last_variable_that_still_needs_to_become_integer_value), recursion_level);
	}

	if (error_code == 0)
	{
		unsigned int number_of_variables_needing_to_be_integer_or_binary_currently_integer_or_binary = 0, last_variable_that_still_needs_to_become_integer_or_binary_index_recursive = 0;

		CountNumberOfVariablesNeedingToBeIntegerOrBinaryThatAlreadyAre(number_of_variables, number_of_variables_required_to_be_integer, number_of_variables_required_to_be_binary, variable_special_requirements, (*basic_feasible_solution), &number_of_variables_needing_to_be_integer_or_binary_currently_integer_or_binary, &last_variable_that_still_needs_to_become_integer_or_binary_index_recursive);

		if (number_of_variables_needing_to_be_integer_or_binary_currently_integer_or_binary < number_of_variables_required_to_be_integer + number_of_variables_required_to_be_binary)
		{
			error_code = CheckIfMoreBranchAndBoundMILPRecursionIsNecessary(maximization_problem, number_of_constraints, number_of_variables, number_of_slack_surplus_variables, number_of_artificial_variables, basic_variables, basic_feasible_solution, tableau_matrix, tableau_current_size, tableau_max_size, number_of_variables_required_to_be_integer, number_of_variables_required_to_be_binary, variable_special_requirements, last_variable_that_still_needs_to_become_integer_or_binary_index_recursive, best_mixed_integer_optimal_value_double, best_mixed_integer_optimal_value, best_mixed_integer_variable_values, best_milp_error_code, recursion_level);
		}
		else
		{
			PrintOptimalResults(number_of_variables, (*tableau_matrix)[0][(*number_of_constraints)][0], (*tableau_matrix)[1][(*number_of_constraints)][0], (*basic_feasible_solution));
			SaveBestMixedIntegerOptimalVariablesAndValues(maximization_problem, (*number_of_constraints), number_of_variables, (*basic_feasible_solution), (*tableau_matrix), best_mixed_integer_optimal_value_double, best_mixed_integer_optimal_value, best_mixed_integer_variable_values, best_milp_error_code);
		}
	}
	else if (error_code == 1)
	{
		printf("BranchAndBoundMILPRecursiveLessOrGreaterThan: Infeasible!\n");
	}
	else if (error_code == 2)
	{
		if ((*best_milp_error_code) == 1)
		{
			(*best_milp_error_code) = 2;
		}
		printf("BranchAndBoundMILPRecursiveLessOrGreaterThan: Unbounded!\n");
	}

	ResetBranchAndBoundMILPRecursiveCountsAndArrays(number_of_variables, branch_and_bound_state_old, tableau_current_size, number_of_constraints, number_of_slack_surplus_variables, number_of_artificial_variables, basic_variables, basic_feasible_solution, tableau_matrix);

	return error_code;
} // end of BranchAndBoundMILPRecursiveLessOrGreaterThan function

/* This function checks if more recursion of MILP branch and bound is necessary */
int CheckIfMoreBranchAndBoundMILPRecursionIsNecessary(int maximization_problem, unsigned int *number_of_constraints, unsigned int number_of_variables, unsigned int *number_of_slack_surplus_variables, unsigned int *number_of_artificial_variables, unsigned int **basic_variables, long long ***basic_feasible_solution, long long ****tableau_matrix, unsigned int *tableau_current_size, unsigned int *tableau_max_size, unsigned int number_of_variables_required_to_be_integer, unsigned int number_of_variables_required_to_be_binary, int *variable_special_requirements, unsigned int last_variable_that_still_needs_to_become_integer_or_binary_index_recursive, double *best_mixed_integer_optimal_value_double, long long *best_mixed_integer_optimal_value, long long **best_mixed_integer_variable_values, int *best_milp_error_code, unsigned int recursion_level)
{
	int error_code = 0;
	
	if (recursion_level < max_branch_and_bound_recursion_depth)
	{
		if (maximization_problem == 1)
		{
			if ((double)(*tableau_matrix)[0][(*number_of_constraints)][0] / (*tableau_matrix)[1][(*number_of_constraints)][0] > (*best_mixed_integer_optimal_value_double))
			{
				/* Keep going down the rabbit hole */
				error_code = BranchAndBoundMILPRecursive(maximization_problem, number_of_constraints, number_of_variables, number_of_slack_surplus_variables, number_of_artificial_variables, basic_variables, basic_feasible_solution, tableau_matrix, tableau_current_size, tableau_max_size, number_of_variables_required_to_be_integer, number_of_variables_required_to_be_binary, variable_special_requirements, last_variable_that_still_needs_to_become_integer_or_binary_index_recursive, best_mixed_integer_optimal_value_double, best_mixed_integer_optimal_value, best_mixed_integer_variable_values, best_milp_error_code, recursion_level + 1);
			}
		}
		else
		{
			if (-(double)(*tableau_matrix)[0][(*number_of_constraints)][0] / (*tableau_matrix)[1][(*number_of_constraints)][0] < (*best_mixed_integer_optimal_value_double))
			{
				/* Keep going down the rabbit hole */
				error_code = BranchAndBoundMILPRecursive(maximization_problem, number_of_constraints, number_of_variables, number_of_slack_surplus_variables, number_of_artificial_variables, basic_variables, basic_feasible_solution, tableau_matrix, tableau_current_size, tableau_max_size, number_of_variables_required_to_be_integer, number_of_variables_required_to_be_binary, variable_special_requirements, last_variable_that_still_needs_to_become_integer_or_binary_index_recursive, best_mixed_integer_optimal_value_double, best_mixed_integer_optimal_value, best_mixed_integer_variable_values, best_milp_error_code, recursion_level + 1);
			}
		}
	}

	return error_code;
} // end of CheckIfMoreBranchAndBoundMILPRecursionIsNecessary function

/* This function saves the best mixed integer optimal variables and values */
void SaveBestMixedIntegerOptimalVariablesAndValues(int maximization_problem, unsigned int number_of_constraints, unsigned int number_of_variables, long long **basic_feasible_solution, long long ***tableau_matrix, double *best_mixed_integer_optimal_value_double, long long *best_mixed_integer_optimal_value, long long **best_mixed_integer_variable_values, int *best_milp_error_code)
{
	unsigned int i;

	if (maximization_problem == 1)
	{
		if ((double)tableau_matrix[0][number_of_constraints][0] / tableau_matrix[1][number_of_constraints][0] > (*best_mixed_integer_optimal_value_double))
		{
			(*best_milp_error_code) = 0;

			(*best_mixed_integer_optimal_value_double) = (double)tableau_matrix[0][number_of_constraints][0] / tableau_matrix[1][number_of_constraints][0];

			best_mixed_integer_optimal_value[0] = tableau_matrix[0][number_of_constraints][0];
			best_mixed_integer_optimal_value[1] = tableau_matrix[1][number_of_constraints][0];

			for (i = 0; i < number_of_variables; i++)
			{
				best_mixed_integer_variable_values[0][i] = basic_feasible_solution[0][i];
				best_mixed_integer_variable_values[1][i] = basic_feasible_solution[1][i];
			} // end of i loop
		}
	}
	else
	{
		if (-(double)tableau_matrix[0][number_of_constraints][0] / tableau_matrix[1][number_of_constraints][0] < (*best_mixed_integer_optimal_value_double))
		{
			(*best_milp_error_code) = 0;

			(*best_mixed_integer_optimal_value_double) = -(double)tableau_matrix[0][number_of_constraints][0] / tableau_matrix[1][number_of_constraints][0];

			best_mixed_integer_optimal_value[0] = -tableau_matrix[0][number_of_constraints][0];
			best_mixed_integer_optimal_value[1] = tableau_matrix[1][number_of_constraints][0];

			for (i = 0; i < number_of_variables; i++)
			{
				best_mixed_integer_variable_values[0][i] = basic_feasible_solution[0][i];
				best_mixed_integer_variable_values[1][i] = basic_feasible_solution[1][i];
			} // end of i loop
		}
	}
	
	return;
} // end of SaveBestMixedIntegerOptimalVariablesAndValues function

/* This function resets the counts and arrays during recursive branch and bound MILP */
void ResetBranchAndBoundMILPRecursiveCountsAndArrays(unsigned int number_of_variables, struct BranchAndBoundState branch_and_bound_state_old, unsigned int *tableau_current_size, unsigned int *number_of_constraints, unsigned int *number_of_slack_surplus_variables, unsigned int *number_of_artificial_variables, unsigned int **basic_variables, long long ***basic_feasible_solution, long long ****tableau_matrix)
{
	unsigned int i, j, k;

	/* Reset counts */
	tableau_current_size[0] = branch_and_bound_state_old.tableau_current_size[0];
	tableau_current_size[1] = branch_and_bound_state_old.tableau_current_size[1];
	(*number_of_constraints) = branch_and_bound_state_old.number_of_constraints;
	(*number_of_slack_surplus_variables) = branch_and_bound_state_old.number_of_slack_surplus_variables;
	(*number_of_artificial_variables) = branch_and_bound_state_old.number_of_artificial_variables;

	/* Reset arrays */
	for (i = 0; i < branch_and_bound_state_old.number_of_constraints; i++)
	{
		(*basic_variables)[i] = branch_and_bound_state_old.basic_variables[i];
	} // end of i loop
	
	for (k = 0; k < 2; k++)
	{
		for (i = 0; i < (number_of_variables + branch_and_bound_state_old.number_of_slack_surplus_variables + branch_and_bound_state_old.number_of_artificial_variables); i++)
		{
			(*basic_feasible_solution)[k][i] = branch_and_bound_state_old.basic_feasible_solution[k][i];
		} // end of i loop
	} // end of k loop

	for (k = 0; k < 2; k++)
	{
		for (i = 0; i < branch_and_bound_state_old.tableau_current_size[0]; i++)
		{
			for (j = 0; j < branch_and_bound_state_old.tableau_current_size[1]; j++)
			{
				(*tableau_matrix)[k][i][j] = branch_and_bound_state_old.tableau_matrix[k][i][j];
			} // end of j loop
		} // end of i loop
	} // end of k loop
	
	return;
} // end of ResetBranchAndBoundMILPRecursiveCountsAndArrays function

/* This function adds a constraint to the previous optimal solution (warm-start) */
int AddConstraint(unsigned int *number_of_constraints, unsigned int number_of_variables, unsigned int *number_of_slack_surplus_variables, unsigned int *number_of_artificial_variables, unsigned int **basic_variables, long long ***basic_feasible_solution, long long ****tableau_matrix, unsigned int *tableau_current_size, unsigned int *tableau_max_size, unsigned int last_variable_that_still_needs_to_become_integer_or_binary_index, int new_constraint_inequality_direction, long long new_constraint_constant, unsigned int recursion_level)
{
	unsigned int i, j;
	int error_code;

	if (new_constraint_inequality_direction == 1) // if less than or equal to
	{
		/* Create new >= constraint for first variable that still needs to become an interger */
		printf("\nAddConstraint: New constraint is x%u <= %lld -> -x%u >= %lld at recursion level %u\n", last_variable_that_still_needs_to_become_integer_or_binary_index + 1, new_constraint_constant, last_variable_that_still_needs_to_become_integer_or_binary_index + 1, -new_constraint_constant, recursion_level);
	}
	else // if greater than or equal to
	{
		/* Create new >= constraint for first variable that still needs to become an interger */
		printf("\nAddConstraint: New constraint is x%u >= %lld at recursion level %u\n", last_variable_that_still_needs_to_become_integer_or_binary_index + 1, new_constraint_constant, recursion_level);
	}

	printf("AddConstraint: Tableau: tableau_current_size = %u x %u & tableau_max_size = %u x %u\n", tableau_current_size[0], tableau_current_size[1], tableau_max_size[0], tableau_max_size[1]);

	/* Increase matrix array size if needed */
	if (tableau_current_size[0] + 2 > tableau_max_size[0]) // if need to realloc rows
	{
		if (tableau_current_size[1] + 2 > tableau_max_size[1]) // if need to realloc cols
		{
			Realloc3DLongLong(tableau_matrix, tableau_max_size[0], tableau_current_size[0] + 2, tableau_max_size[1], tableau_current_size[1] + 2);
			tableau_max_size[1] = tableau_current_size[1] + 2;

			Realloc2DLongLong(basic_feasible_solution, number_of_variables + (*number_of_slack_surplus_variables) + (*number_of_artificial_variables), number_of_variables + (*number_of_slack_surplus_variables) + (*number_of_artificial_variables) + 2);
		}
		else // if do NOT need to realloc cols
		{
			Realloc3DLongLong(tableau_matrix, tableau_max_size[0], tableau_current_size[0] + 2, tableau_max_size[1], tableau_max_size[1]);
		}
		tableau_max_size[0] = tableau_current_size[0] + 2;
	}
	else // if do NOT need to realloc rows
	{
		if (tableau_current_size[1] + 2 > tableau_max_size[1]) // if need to realloc cols
		{
			Realloc3DLongLong(tableau_matrix, tableau_max_size[0], tableau_max_size[0], tableau_max_size[1], tableau_current_size[1] + 2);
			tableau_max_size[1] = tableau_current_size[1] + 2;

			Realloc2DLongLong(basic_feasible_solution, number_of_variables + (*number_of_slack_surplus_variables) + (*number_of_artificial_variables), number_of_variables + (*number_of_slack_surplus_variables) + (*number_of_artificial_variables) + 2);
		}
	}

	tableau_current_size[0] += 2;
	tableau_current_size[1] += 2;

	/* Reset columns for new surplus and artificial variable */
	for (j = tableau_current_size[1] - 2; j < tableau_current_size[1]; j++)
	{
		for (i = 0; i < tableau_current_size[0]; i++)
		{
			(*tableau_matrix)[0][i][j] = 0;
			(*tableau_matrix)[1][i][j] = 1;
		} // end of i loop
	} // end of j loop

	/* Reset bottom row for new phase 1 objective function */
	for (j = 0; j < tableau_current_size[1] - 2; j++)
	{
		(*tableau_matrix)[0][tableau_current_size[0] - 1][j] = 0;
		(*tableau_matrix)[1][tableau_current_size[0] - 1][j] = 1;
	} // end of j loop

	/* Move phase 2 objective function down 1 row */
	for (j = 0; j < tableau_current_size[1] - 2; j++)
	{
		(*tableau_matrix)[0][tableau_current_size[0] - 2][j] = (*tableau_matrix)[0][tableau_current_size[0] - 3][j];
		(*tableau_matrix)[1][tableau_current_size[0] - 2][j] = (*tableau_matrix)[1][tableau_current_size[0] - 3][j];
		
		/* Reset new constraint row */
		(*tableau_matrix)[0][tableau_current_size[0] - 3][j] = 0;
		(*tableau_matrix)[1][tableau_current_size[0] - 3][j] = 1;
	} // end of j loop

	if (new_constraint_inequality_direction == 1) // if less than or equal to
	{
		/* Add new constraint constant */
		(*tableau_matrix)[0][tableau_current_size[0] - 3][0] = -new_constraint_constant;
		(*tableau_matrix)[1][tableau_current_size[0] - 3][0] = 1;

		/* Add variable into its column in the new constraint row */
		(*tableau_matrix)[0][tableau_current_size[0] - 3][last_variable_that_still_needs_to_become_integer_or_binary_index + 1] = -1;
		(*tableau_matrix)[1][tableau_current_size[0] - 3][last_variable_that_still_needs_to_become_integer_or_binary_index + 1] = 1;
	}
	else // if greater than or equal to
	{
		/* Add new constraint constant */
		(*tableau_matrix)[0][tableau_current_size[0] - 3][0] = new_constraint_constant;
		(*tableau_matrix)[1][tableau_current_size[0] - 3][0] = 1;

		/* Add variable into its column in the new constraint row */
		(*tableau_matrix)[0][tableau_current_size[0] - 3][last_variable_that_still_needs_to_become_integer_or_binary_index + 1] = 1;
		(*tableau_matrix)[1][tableau_current_size[0] - 3][last_variable_that_still_needs_to_become_integer_or_binary_index + 1] = 1;
	}

	/* Add surplus variable into second to last column */
	(*tableau_matrix)[0][tableau_current_size[0] - 3][tableau_current_size[1] - 2] = -1;
	(*tableau_matrix)[1][tableau_current_size[0] - 3][tableau_current_size[1] - 2] = 1;

	/* Add artificial variable into last column in new constraint row */
	(*tableau_matrix)[0][tableau_current_size[0] - 3][tableau_current_size[1] - 1] = 1;
	(*tableau_matrix)[1][tableau_current_size[0] - 3][tableau_current_size[1] - 1] = 1;

	/* Add artificial variable into last column in artificial objective row */
	(*tableau_matrix)[0][tableau_current_size[0] - 1][tableau_current_size[1] - 1] = -1;
	(*tableau_matrix)[1][tableau_current_size[0] - 1][tableau_current_size[1] - 1] = 1;

	if (new_constraint_inequality_direction == 1) // if less than or equal to
	{
		/* Add variable's basic row to new constraint row */
		for (i = 0; i < (*number_of_constraints); i++)
		{
			if ((*basic_variables)[i] == last_variable_that_still_needs_to_become_integer_or_binary_index + 1)
			{
				for (j = 0; j < tableau_current_size[1]; j++)
				{
					LongLongRationalAddition((*tableau_matrix)[0][tableau_current_size[0] - 3][j], (*tableau_matrix)[1][tableau_current_size[0] - 3][j], (*tableau_matrix)[0][i][j], (*tableau_matrix)[1][i][j], &(*tableau_matrix)[0][tableau_current_size[0] - 3][j], &(*tableau_matrix)[1][tableau_current_size[0] - 3][j]);
				} // end of j loop

				break; // break i loop, we found the variable already
			}
		} // end of i loop
	}
	else // if greater than or equal to
	{
		/* Subtract variable's basic row from new constraint row */
		for (i = 0; i < (*number_of_constraints); i++)
		{
			if ((*basic_variables)[i] == last_variable_that_still_needs_to_become_integer_or_binary_index + 1)
			{
				for (j = 0; j < tableau_current_size[1]; j++)
				{
					LongLongRationalAddition((*tableau_matrix)[0][tableau_current_size[0] - 3][j], (*tableau_matrix)[1][tableau_current_size[0] - 3][j], -(*tableau_matrix)[0][i][j], (*tableau_matrix)[1][i][j], &(*tableau_matrix)[0][tableau_current_size[0] - 3][j], &(*tableau_matrix)[1][tableau_current_size[0] - 3][j]);
				} // end of j loop

				break; // break i loop, we found the variable already
			}
		} // end of i loop
	}

	(*number_of_artificial_variables)++;

	/* This function reallocates more memory to the passed 1d array */
	Realloc1DUnsignedInt(basic_variables, (*number_of_constraints), (*number_of_constraints) + 1, tableau_current_size[1] - 1);

	(*number_of_constraints)++;
	(*number_of_slack_surplus_variables)++;

	/* This function updates the basic feasible solution */
	UpdateBasicFeasibleSolution(tableau_current_size[1] - 1, (*number_of_constraints), (*basic_variables), (*basic_feasible_solution), (*tableau_matrix));

	/* This function finds the optimal solution for the given variables and constraints */
	error_code = SimplexAlgorithm((*number_of_constraints), number_of_variables, (*number_of_slack_surplus_variables), number_of_artificial_variables, (*basic_variables), (*basic_feasible_solution), (*tableau_matrix), tableau_current_size);
	
	return error_code;
} // end of AddConstraint function

/*********************************************************************************/
/******************************* HELPER FUNCTIONS ********************************/
/*********************************************************************************/

/* This function prints the optimal objective function value and the variable values */
void PrintOptimalResults(unsigned int number_of_variables, long long optimal_objective_function_value_numerator, long long optimal_objective_function_value_denominator, long long **optimal_variable_values)
{
	unsigned int i;

	printf("PrintOptimalResults: optimal_objective_function_value = %lld\t%lld\n", optimal_objective_function_value_numerator, optimal_objective_function_value_denominator);
	printf("PrintOptimalResults: (double)optimal_objective_function_value = %.16f\n", (double)optimal_objective_function_value_numerator / optimal_objective_function_value_denominator);
	printf("PrintOptimalResults: optimal_variable_values:\n");
	printf("variable_index\tvariable_numerator\tvaraible_denominator\tvariable_value\n");
	for (i = 0; i < number_of_variables; i++)
	{
		printf("%u\t%lld\t%lld\t%.16f\n", i, optimal_variable_values[0][i], optimal_variable_values[1][i], (double)optimal_variable_values[0][i] / optimal_variable_values[1][i]);
	} // end of i loop
	
	return;
} // end of PrintOptimalResults function

/* This function reallocates more memory to the passed 1d array */
void Realloc1DUnsignedInt(unsigned int **array1d, unsigned int oldsize1d, unsigned int newsize1d, unsigned int initializer)
{
	unsigned int i;

	/* Create temp pointer */
	unsigned int *temp = NULL;
	temp = realloc(*array1d, sizeof(unsigned int) * newsize1d);
	if (temp == NULL)
	{
		printf("Realloc1DUnsignedInt: ERROR: out of memory\n");
	}
	else
	{
		*array1d = temp;
	}

	/* Initialize new elements */
	for (i = oldsize1d; i < newsize1d; i++)
	{
		(*array1d)[i] = initializer;
	} // end of i loop
	
	return;
} // end of Realloc1DUnsignedInt function

/* This function reallocates more memory to the passed 2d array */
void Realloc2DLongLong(long long ***array2d, unsigned int oldsize1d, unsigned int newsize1d)
{
	unsigned int i, k;

	/* Create temp pointer */
	long long **temp = NULL;
	temp = malloc(sizeof(long long*) * 2);
	for (k = 0; k < 2; k++)
	{
		temp[k] = malloc(sizeof(long long) * newsize1d);
	} // end of k loop

	for (k = 0; k < 2; k++)
	{
		for (i = 0; i < oldsize1d; i++)
		{
			temp[k][i] = (*array2d)[k][i];
		} // end of i loop
	} // end of k loop

	/* Free 1d array */
	for (k = 0; k < 2; k++)
	{
		free((*array2d)[k]);
	} // end of k loop
	free(*array2d);

	(*array2d) = malloc(sizeof(long long*) * 2);
	for (k = 0; k < 2; k++)
	{
		(*array2d)[k] = malloc(sizeof(long long) * newsize1d);
	} // end of k loop

	for (k = 0; k < 2; k++)
	{
		for (i = 0; i < oldsize1d; i++)
		{
			(*array2d)[k][i] = temp[k][i];
		} // end of i loop
	} // end of k loop

	for (k = 0; k < 2; k++)
	{
		free(temp[k]);
	} // end of k loop
	free(temp);

	/* Initialize new elements */
	for (i = oldsize1d; i < newsize1d; i++)
	{
		(*array2d)[0][i] = 0;
		(*array2d)[1][i] = 1;
	} // end of i loop
	
	return;
} // end of Realloc2DLongLong function

/* This function reallocates more memory to the passed 3d array */
void Realloc3DLongLong(long long ****array3d, unsigned int oldsize1d, unsigned int newsize1d, unsigned int oldsize2d, unsigned int newsize2d)
{
	unsigned int i, j, k;

	/* Create temp pointer */
	long long ***temp = NULL;

	if (newsize1d != oldsize1d || newsize2d != oldsize2d) // if there is any increase in size for either dimension
	{
		temp = malloc(sizeof(long long*) * 2);
		for (k = 0; k < 2; k++)
		{
			temp[k] = malloc(sizeof(long long*) * newsize1d);

			for (i = 0; i < newsize1d; i++)
			{
				temp[k][i] = malloc(sizeof(long long) * newsize2d);
			} // end of i loop
		} // end of k loop

		for (k = 0; k < 2; k++)
		{
			for (i = 0; i < oldsize1d; i++)
			{
				for (j = 0; j < oldsize2d; j++)
				{
					temp[k][i][j] = (*array3d)[k][i][j];
				} // end of j loop
			} // end of i loop
		} // end of k loop

		/* Free 2d array */
		for (k = 0; k < 2; k++)
		{
			for (i = 0; i < oldsize1d; i++)
			{
				free((*array3d)[k][i]);
			} // end of i loop
			free((*array3d)[k]);
		} // end of k loop
		free(*array3d);

		(*array3d) = malloc(sizeof(long long*) * 2);
		for (k = 0; k < 2; k++)
		{
			(*array3d)[k] = malloc(sizeof(long long*) * newsize1d);
			for (i = 0; i < newsize1d; i++)
			{
				(*array3d)[k][i] = malloc(sizeof(double) * newsize2d);
			} // end of i loop
		} // end of k loop

		for (k = 0; k < 2; k++)
		{
			for (i = 0; i < oldsize1d; i++)
			{
				for (j = 0; j < oldsize2d; j++)
				{
					(*array3d)[k][i][j] = temp[k][i][j];
				} // end of j loop
			} // end of i loop
		} // end of k loop

		for (k = 0; k < 2; k++)
		{
			for (i = 0; i < newsize1d; i++)
			{
				free(temp[k][i]);
			} // end of i loop
			free(temp[k]);
		} // end of k loop
		free(temp);
	} // end of if there is any increase in size for either dimension

	/* Initialize new elements */
	if ((newsize1d != oldsize1d && newsize2d == oldsize2d) || (newsize1d != oldsize1d && newsize2d != oldsize2d)) // if only 1st dim's size has changed or both
	{
		for (i = oldsize1d; i < newsize1d; i++) // just new rows
		{
			for (j = 0; j < newsize2d; j++) // fill in all columns
			{
				(*array3d)[0][i][j] = 0;
				(*array3d)[1][i][j] = 1;
			} // end of j loop
		} // end of i loop
	} // end of if only 1st dim's size has changed or both

	if ((newsize1d == oldsize1d && newsize2d != oldsize2d) || (newsize1d != oldsize1d && newsize2d != oldsize2d)) // if only 2nd dim's size has changed or both
	{
		for (i = 0; i < newsize1d; i++) // for all rows
		{
			for (j = oldsize2d; j < newsize2d; j++) // fill in just new columns
			{
				(*array3d)[0][i][j] = 0;
				(*array3d)[1][i][j] = 1;
			} // end of j loop
		} // end of i loop
	} // end of if only 2nd dim's size has changed or both
	
	return;
} // end of Realloc3DLongLong function

/* This function adds long long rationals and keeps their numerators and denominators each as small as possible */
void LongLongRationalAddition(long long A_numerator, long long A_denominator, long long B_numerator, long long B_denominator, long long *C_numerator, long long *C_denominator)
{
	long long temp_gcd;

	(*C_numerator) = A_numerator * B_denominator + A_denominator * B_numerator;
	(*C_denominator) = A_denominator * B_denominator;
	temp_gcd = GreastestCommonDenominator((*C_numerator), (*C_denominator));
	(*C_numerator) /= temp_gcd;
	(*C_denominator) /= temp_gcd;
	
	return;
} // end of LongLongRationalAddition function

/* This function divides long long rationals and keeps their numerators and denominators each as small as possible */
void LongLongRationalDivision(long long A_numerator, long long A_denominator, long long B_numerator, long long B_denominator, long long *C_numerator, long long *C_denominator)
{
	long long temp_gcd;

	(*C_numerator) = A_numerator * B_denominator;
	(*C_denominator) = A_denominator * B_numerator;
	temp_gcd = GreastestCommonDenominator((*C_numerator), (*C_denominator));
	(*C_numerator) /= temp_gcd;
	(*C_denominator) /= temp_gcd;

	if ((*C_denominator) < 0) // if denominator is negative
	{
		(*C_numerator) = -(*C_numerator);
		(*C_denominator) = -(*C_denominator);
	} // end of if denominator is negative
	
	return;
} // end of LongLongRationalDivision function

/* This function finds the greatest common denominator between a and b */
long long GreastestCommonDenominator(long long a, long long b)
{
	long long t;
	while (b != 0)
	{
		t = b;
		b = a % b;
		a = t;
	}
	return llabs(a);
} // end of GreastestCommonDenominator function