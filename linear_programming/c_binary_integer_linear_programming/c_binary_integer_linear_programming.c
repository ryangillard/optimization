#include <stdio.h>      /* printf, scanf, puts */
#include <stdlib.h>     /* realloc, free, exit, NULL */
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>

int maximization_problem = 1; // default is maximization
double epsilon = 0.000000001; // the maximum allowed calculation error when determining integers

/*********************************************************************************/
/********************************** PROTOTYPES ***********************************/
/*********************************************************************************/

/*********************************************************************************/
/********************************* READ INPUTS ***********************************/
/*********************************************************************************/

/* This function reads the number of constraints and variables */
void ReadNumberOfConstraintsAndVariables(unsigned int *initial_number_of_constraints, unsigned int *number_of_constraints, unsigned int *number_of_variables);

/* This function reads in if the objective function is going to be a maximization or minimization problem */
void ReadObjectiveFunctionMaximizationOrMinimization(int *maximization_problem);

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
void CreateOptimalObjectiveFunctionAndVariableValues(unsigned int number_of_variables, long long *optimal_objective_function_value, long long ***optimal_variable_values);

/*********************************************************************************/
/********************** BINARY INTEGER LINEAR PROGRAMMING ************************/
/*********************************************************************************/

/* This function performs binary integer programming */
int BinaryIntegerProgramming(unsigned int number_of_constraints, unsigned int number_of_variables, long long *objective_function_initial_constant, long long **objective_function_coefficient_vector, int *constraint_inequality_direction_vector, long long **constraint_constant_vector, long long ***constraint_coefficient_matrix, long long *optimal_objective_function_value, long long **optimal_variable_values);

/* This function created the sorted objective function arrays since that is needed for the Balas Additive Algorithm */
void CreateSortedObjectiveFunctionArrays(unsigned int number_of_variables, long long **objective_function_coefficient_vector, double **sorted_objective_function_coefficient_vector_value_double, unsigned int **sorted_objective_function_coefficient_vector_index, long long ***sorted_objective_function_coefficient_vector_value);

/* This function applies quicksort on doubles and keeps track of the change in indexing */
void QuickSortDouble(unsigned int n, double *a, unsigned int *index);

/* This function creates the binary constraint arrays */
void CreateBinaryConstraintArrays(unsigned int number_of_binary_constraints, unsigned int number_of_variables, long long ***binary_constraint_constant_vector, long long ****binary_constraint_coefficient_matrix);

/* This function initializes the binary constraint coefficient matrix from the input constraint coefficient matrix based on inequality direction and reversed sign objective function coefficients */
void InitializeBinaryConstraintCoefficientMatrix(unsigned int number_of_constraints, unsigned int number_of_variables, int *constraint_inequality_direction_vector, long long **constraint_constant_vector, long long ***constraint_coefficient_matrix, unsigned int *sorted_objective_function_coefficient_vector_index, int *binary_objective_function_coefficient_reversed_sign_vector, long long **binary_constraint_constant_vector, long long ***binary_constraint_coefficient_matrix);

/* This function initializes the binary constraint coefficient matrix from the input constraint coefficient matrix based on less than inequality direction and reversed sign objective function coefficients */
unsigned int InitializeBinaryConstraintCoefficientMatrixLessThan(unsigned int number_of_variables, long long **constraint_constant_vector, long long ***constraint_coefficient_matrix, unsigned int *sorted_objective_function_coefficient_vector_index, int *binary_objective_function_coefficient_reversed_sign_vector, long long **binary_constraint_constant_vector, long long ***binary_constraint_coefficient_matrix, unsigned int i, unsigned int k);

/* This function initializes the binary constraint coefficient matrix from the input constraint coefficient matrix based on greater than inequality direction and reversed sign objective function coefficients */
unsigned int InitializeBinaryConstraintCoefficientMatrixGreaterThan(unsigned int number_of_variables, long long **constraint_constant_vector, long long ***constraint_coefficient_matrix, unsigned int *sorted_objective_function_coefficient_vector_index, int *binary_objective_function_coefficient_reversed_sign_vector, long long **binary_constraint_constant_vector, long long ***binary_constraint_coefficient_matrix, unsigned int i, unsigned int k);

/* This function initiates branch and bound for binary variables */
void BalasAdditiveBranchAndBoundBinary(unsigned int number_of_binary_constraints, unsigned int number_of_variables, long long **binary_objective_function_coefficient_vector, long long **binary_constraint_constant_vector, long long ***binary_constraint_coefficient_matrix, int *current_binary_feasible_solution, int *best_binary_feasible_solution, double *best_binary_optimal_value_double, long long *best_binary_optimal_value);

/* This function evaluates branch and bound for binary variables */
void BalasAdditiveBranchAndBoundBinaryEvaluation(unsigned int number_of_binary_constraints, unsigned int number_of_variables, long long **binary_objective_function_coefficient_vector, long long **binary_constraint_constant_vector, long long ***binary_constraint_coefficient_matrix, int *current_binary_feasible_solution, int *best_binary_feasible_solution, double *best_binary_optimal_value_double, long long *best_binary_optimal_value, unsigned int variable_index);

/* This function performs branch and bound recursively for binary variables */
void BalasAdditiveBranchAndBoundBinaryRecursive(unsigned int number_of_binary_constraints, unsigned int number_of_variables, long long **binary_objective_function_coefficient_vector, long long **binary_constraint_constant_vector, long long ***binary_constraint_coefficient_matrix, int *current_binary_feasible_solution, int *best_binary_feasible_solution, double *best_binary_optimal_value_double, long long *best_binary_optimal_value, unsigned int variable_index);

/* This function writes to disk the final BILP array values */
void WriteFinalBILPValues(unsigned int number_of_variables, long long *optimal_objective_function_value, long long **optimal_variable_values);

/*********************************************************************************/
/******************************* HELPER FUNCTIONS ********************************/
/*********************************************************************************/

/* This function prints the optimal objective function value and the variable values */
void PrintOptimalResults(unsigned int number_of_variables, double optimal_objective_function_value, long long **optimal_variable_values);

// This function reallocates more memory to the passed 1d array
void Realloc1DUnsignedInt(unsigned int **array1d, unsigned int oldsize1d, unsigned int newsize1d, unsigned int initializer);

// This function reallocates more memory to the passed 1d array
void Realloc1DLongLong(long long ***array1d, unsigned int oldsize1d, unsigned int newsize1d);

// This function reallocates more memory to the passed 2d array
void Realloc2DLongLong(long long ****array2d, unsigned int oldsize1d, unsigned int newsize1d, unsigned int oldsize2d, unsigned int newsize2d);

/* This function adds long long rationals and keeps their numerators and denominators each as small as possible */
void LongLongRationalAddition(long long A_numerator, long long A_denominator, long long B_numerator, long long B_denominator, long long *C_numerator, long long *C_denominator);

/* This function divides long long rationals and keeps their numerators and denominators each as small as possible */
void LongLongRationalDivision(long long A_numerator, long long A_denominator, long long B_numerator, long long B_denominator, long long *C_numerator, long long *C_denominator);

// Find the greatest common denominator between a and b
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

	/* Sizes */
	unsigned int initial_number_of_constraints = 0, number_of_constraints = 0, number_of_variables = 0;

	ReadNumberOfConstraintsAndVariables(&initial_number_of_constraints, &number_of_constraints, &number_of_variables);

	/* Optimization problem type */
	ReadObjectiveFunctionMaximizationOrMinimization(&maximization_problem);

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

	CreateOptimalObjectiveFunctionAndVariableValues(number_of_variables, optimal_objective_function_value, &optimal_variable_values);

	/*********************************************************************************/
	/********************** BINARY INTEGER LINEAR PROGRAMMING ************************/
	/*********************************************************************************/

	error_code = BinaryIntegerProgramming(number_of_constraints, number_of_variables, objective_function_initial_constant, objective_function_coefficient_vector, constraint_inequality_direction_vector, constraint_constant_vector, constraint_coefficient_matrix, optimal_objective_function_value, optimal_variable_values);

	if (error_code == 0)
	{
		printf("\n\n************************************************************************************************************\n");
		printf("************************************************************************************************************\n");
		printf("*********************************************** FINAL SOLUTION ******************************************\n");
		printf("************************************************************************************************************\n");
		printf("************************************************************************************************************\n\n\n");

		PrintOptimalResults(number_of_variables, (double)optimal_objective_function_value[0] / optimal_objective_function_value[1], optimal_variable_values);
	}
	else
	{
		printf("main: Infeasible problem!\n");
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
} // end of ReadNumberOfConstraintsAndVariables function

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
} // end of ReadObjectiveFunctionMaximizationOrMinimization function

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
} // end of ReadConstraintCoefficientMatrix function

/*********************************************************************************/
/******************************** OPTIMAL VALUES *********************************/
/*********************************************************************************/

/* This function creates the arrays to store the optimal objective function value and the corresponding variable values */
void CreateOptimalObjectiveFunctionAndVariableValues(unsigned int number_of_variables, long long *optimal_objective_function_value, long long ***optimal_variable_values)
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
} // end of CreateOptimalObjectiveFunctionAndVariableValues function

/*********************************************************************************/
/********************** BINARY INTEGER LINEAR PROGRAMMING ************************/
/*********************************************************************************/

/* This function performs binary integer programming */
int BinaryIntegerProgramming(unsigned int number_of_constraints, unsigned int number_of_variables, long long *objective_function_initial_constant, long long **objective_function_coefficient_vector, int *constraint_inequality_direction_vector, long long **constraint_constant_vector, long long ***constraint_coefficient_matrix, long long *optimal_objective_function_value, long long **optimal_variable_values)
{
	unsigned int i, j, k;
	int error_code = 0;

	printf("\n\n************************************************************************************************************\n");
	printf("************************************************************************************************************\n");
	printf("****************************************** BINARY BRANCH AND BOUND *****************************************\n");
	printf("************************************************************************************************************\n");
	printf("************************************************************************************************************\n\n\n");

	/* Convert max to min by multiplying by -1, have all objective function coefficients be positive by xi -> (1 - yi), and have all constraints be >= by multiplying by -1 (negative constants allowed) */
	if (maximization_problem == 1)
	{
		objective_function_initial_constant[0] *= -1;
	}

	int *binary_objective_function_coefficient_reversed_sign_vector;
	binary_objective_function_coefficient_reversed_sign_vector = malloc(sizeof(int) * number_of_variables);
	for (i = 0; i < number_of_variables; i++)
	{
		binary_objective_function_coefficient_reversed_sign_vector[i] = 0;
	} // end of i loop

	printf("BinaryIntegerProgramming: binary_objective_function_coefficient_reversed_sign_vector:\n");
	for (i = 0; i < number_of_variables; i++)
	{
		if (maximization_problem == 1)
		{
			if (objective_function_coefficient_vector[0][i] > 0)
			{
				binary_objective_function_coefficient_reversed_sign_vector[i] = 1;

				LongLongRationalAddition(objective_function_initial_constant[0], objective_function_initial_constant[1], -objective_function_coefficient_vector[0][i], objective_function_coefficient_vector[1][i], &objective_function_initial_constant[0], &objective_function_initial_constant[1]);
			}
			else if (objective_function_coefficient_vector[0][i] < 0)
			{
				objective_function_coefficient_vector[0][i] *= -1;
			}
		}
		else
		{
			if (objective_function_coefficient_vector[0][i] < 0)
			{
				objective_function_coefficient_vector[0][i] *= -1;

				binary_objective_function_coefficient_reversed_sign_vector[i] = 1;

				LongLongRationalAddition(objective_function_initial_constant[0], objective_function_initial_constant[1], objective_function_coefficient_vector[0][i], objective_function_coefficient_vector[1][i], &objective_function_initial_constant[0], &objective_function_initial_constant[1]);
			}
		}

		printf("%d\n", binary_objective_function_coefficient_reversed_sign_vector[i]);
	} // end of i loop

	/* For the Balas Additive Algorithm, we want to sort the objective function coefficients in ascending order */
	double *sorted_objective_function_coefficient_vector_value_double;
	unsigned int *sorted_objective_function_coefficient_vector_index; // index hash map
	long long **sorted_objective_function_coefficient_vector_value;

	CreateSortedObjectiveFunctionArrays(number_of_variables, objective_function_coefficient_vector, &sorted_objective_function_coefficient_vector_value_double, &sorted_objective_function_coefficient_vector_index, &sorted_objective_function_coefficient_vector_value);

	/* Break equal to constraints into two constraints, so one becomes two */
	unsigned int number_of_binary_constraints = number_of_constraints;
	for (i = 0; i < number_of_constraints; i++)
	{
		if (constraint_inequality_direction_vector[i] == 0) // if equal to
		{
			number_of_binary_constraints++;
		}
	} // end of i loop

	long long **binary_constraint_constant_vector;
	long long ***binary_constraint_coefficient_matrix;

	CreateBinaryConstraintArrays(number_of_binary_constraints, number_of_variables, &binary_constraint_constant_vector, &binary_constraint_coefficient_matrix);

	InitializeBinaryConstraintCoefficientMatrix(number_of_constraints, number_of_variables, constraint_inequality_direction_vector, constraint_constant_vector, constraint_coefficient_matrix, sorted_objective_function_coefficient_vector_index, binary_objective_function_coefficient_reversed_sign_vector, binary_constraint_constant_vector, binary_constraint_coefficient_matrix);

	/* Create variables to keep track of the best optimal solutions */
	double best_binary_optimal_value_double = DBL_MAX;

	long long best_binary_optimal_value[2];
	best_binary_optimal_value[0] = 0;
	best_binary_optimal_value[1] = 1;

	int *current_binary_feasible_solution;
	current_binary_feasible_solution = malloc(sizeof(double) * number_of_variables);
	for (i = 0; i < number_of_variables; i++)
	{
		current_binary_feasible_solution[i] = 0;
	} // end of i loop

	int *best_binary_feasible_solution;
	best_binary_feasible_solution = malloc(sizeof(double) * number_of_variables);
	for (i = 0; i < number_of_variables; i++)
	{
		best_binary_feasible_solution[i] = 0;
	} // end of i loop

	/* Now that everything is setup, run the Balas Additive Branch and Bound algorithm */
	BalasAdditiveBranchAndBoundBinary(number_of_binary_constraints, number_of_variables, sorted_objective_function_coefficient_vector_value, binary_constraint_constant_vector, binary_constraint_coefficient_matrix, current_binary_feasible_solution, best_binary_feasible_solution, &best_binary_optimal_value_double, best_binary_optimal_value);

	if (best_binary_optimal_value_double != DBL_MAX) // if best binary optimal value was changed, thus a feasible solution was found
	{
		error_code = 0;

		int *best_unsorted_binary_feasible_solution;
		best_unsorted_binary_feasible_solution = malloc(sizeof(int) * number_of_variables);
		for (i = 0; i < number_of_variables; i++)
		{
			if (binary_objective_function_coefficient_reversed_sign_vector[sorted_objective_function_coefficient_vector_index[i]] == 0)
			{
				best_unsorted_binary_feasible_solution[sorted_objective_function_coefficient_vector_index[i]] = best_binary_feasible_solution[i];
			}
			else
			{
				best_unsorted_binary_feasible_solution[sorted_objective_function_coefficient_vector_index[i]] = 1 - best_binary_feasible_solution[i];
			}
		} // end of i loop

		if (maximization_problem == 1)
		{
			best_binary_optimal_value_double = -(((double)objective_function_initial_constant[0] / objective_function_initial_constant[1]) + best_binary_optimal_value_double);
			LongLongRationalAddition(-best_binary_optimal_value[0], best_binary_optimal_value[1], -objective_function_initial_constant[0], objective_function_initial_constant[1], &best_binary_optimal_value[0], &best_binary_optimal_value[1]);
		}
		else
		{
			best_binary_optimal_value_double = (((double)objective_function_initial_constant[0] / objective_function_initial_constant[1]) + best_binary_optimal_value_double);
			LongLongRationalAddition(best_binary_optimal_value[0], best_binary_optimal_value[1], objective_function_initial_constant[0], objective_function_initial_constant[1], &best_binary_optimal_value[0], &best_binary_optimal_value[1]);
		}

		objective_function_initial_constant[0] = best_binary_optimal_value[0];
		objective_function_initial_constant[1] = best_binary_optimal_value[1];

		/* Update overall optimal solution */
		optimal_objective_function_value[0] = best_binary_optimal_value[0];
		optimal_objective_function_value[1] = best_binary_optimal_value[1];

		for (i = 0; i < number_of_variables; i++)
		{
			optimal_variable_values[0][i] = best_unsorted_binary_feasible_solution[i];
		} // end of i loop

		WriteFinalBILPValues(number_of_variables, optimal_objective_function_value, optimal_variable_values);

		/* Free dynamic memory */
		free(best_unsorted_binary_feasible_solution);
	}
	else
	{
		error_code = 1; // infeasible
	}

	/* Free dynamic memory */
	free(best_binary_feasible_solution);
	free(current_binary_feasible_solution);

	for (k = 0; k < 2; k++)
	{
		for (i = 0; i < number_of_binary_constraints; i++)
		{
			free(binary_constraint_coefficient_matrix[k][i]);
		} // end of i loop
		free(binary_constraint_coefficient_matrix[k]);
		free(binary_constraint_constant_vector[k]);
	} // end of k loop
	free(binary_constraint_coefficient_matrix);
	free(binary_constraint_constant_vector);

	for (k = 0; k < 2; k++)
	{
		free(sorted_objective_function_coefficient_vector_value[k]);
	} // end of k loop
	free(sorted_objective_function_coefficient_vector_value);
	free(sorted_objective_function_coefficient_vector_index);
	free(sorted_objective_function_coefficient_vector_value_double);
	free(binary_objective_function_coefficient_reversed_sign_vector);

	return error_code;
} // end of BinaryIntegerProgramming function

/* This function created the sorted objective function arrays since that is needed for the Balas Additive Algorithm */
void CreateSortedObjectiveFunctionArrays(unsigned int number_of_variables, long long **objective_function_coefficient_vector, double **sorted_objective_function_coefficient_vector_value_double, unsigned int **sorted_objective_function_coefficient_vector_index, long long ***sorted_objective_function_coefficient_vector_value)
{
	unsigned int i, j, k;

	(*sorted_objective_function_coefficient_vector_value_double) = malloc(sizeof(double) * number_of_variables);

	(*sorted_objective_function_coefficient_vector_index) = malloc(sizeof(unsigned int) * number_of_variables);

	for (i = 0; i < number_of_variables; i++)
	{
		(*sorted_objective_function_coefficient_vector_value_double)[i] = (double)objective_function_coefficient_vector[0][i] / objective_function_coefficient_vector[1][i];
		(*sorted_objective_function_coefficient_vector_index)[i] = i;
	} // end of i loop

	/* Sort based on value and keep track of the indices */
	QuickSortDouble(number_of_variables, (*sorted_objective_function_coefficient_vector_value_double), (*sorted_objective_function_coefficient_vector_index));

	(*sorted_objective_function_coefficient_vector_value) = malloc(sizeof(long long*) * 2);
	for (k = 0; k < 2; k++)
	{
		(*sorted_objective_function_coefficient_vector_value)[k] = malloc(sizeof(long long) * number_of_variables);
	} // end of k loop

	for (k = 0; k < 2; k++)
	{
		for (j = 0; j < number_of_variables; j++)
		{
			(*sorted_objective_function_coefficient_vector_value)[k][j] = objective_function_coefficient_vector[k][(*sorted_objective_function_coefficient_vector_index)[j]];
		} // end of j loop
	} // end of k loop

	printf("CreateSortedObjectiveFunctionArrays: sorted_objective_function_coefficient_vector_value_double:\n");
	for (j = 0; j < number_of_variables; j++)
	{
		printf("%lf\t", (*sorted_objective_function_coefficient_vector_value_double)[j]);
	} // end of i loop
	printf("\n");

	printf("CreateSortedObjectiveFunctionArrays: sorted_objective_function_coefficient_vector_index:\n");
	for (j = 0; j < number_of_variables; j++)
	{
		printf("%u\t", (*sorted_objective_function_coefficient_vector_index)[j]);
	} // end of i loop
	printf("\n");
} // end of CreateSortedObjectiveFunctionArrays function

/* This function applies quicksort and keeps track of the change in indexing */
void QuickSortDouble(unsigned int n, double *a, unsigned int *index)
{
	unsigned int i, j;
	double p, t;

	if (n < 2)
	{
		return;
	}

 // Choose pivot as midpoint of array
	p = a[n / 2];

	for (i = 0, j = n - 1;; i++, j--)
	{
		while (a[i] < p)
		{
			i++;
		}

		while (p < a[j])
		{
			j--;
		}

		if (i >= j)
		{
			break;
		}

	 // Swap
		t = a[i];
		a[i] = a[j];
		a[j] = t;

		t = index[i];
		index[i] = index[j];
		index[j] = t;
	} // end of i, j loop

	QuickSortDouble(i, a, index);
	QuickSortDouble(n - i, a + i, index + i);
} // end of QuickSortDouble function

/* This function creates the binary constraint arrays */
void CreateBinaryConstraintArrays(unsigned int number_of_binary_constraints, unsigned int number_of_variables, long long ***binary_constraint_constant_vector, long long ****binary_constraint_coefficient_matrix)
{
	unsigned int i, j, k;

	(*binary_constraint_constant_vector) = malloc(sizeof(long long*) * 2);
	for (k = 0; k < 2; k++)
	{
		(*binary_constraint_constant_vector)[k] = malloc(sizeof(long long) * number_of_binary_constraints);
	} // end of k loop


	(*binary_constraint_coefficient_matrix) = malloc(sizeof(long long*) * 2);
	for (k = 0; k < 2; k++)
	{
		(*binary_constraint_coefficient_matrix)[k] = malloc(sizeof(long long*) * number_of_binary_constraints);
		for (i = 0; i < number_of_binary_constraints; i++)
		{
			(*binary_constraint_coefficient_matrix)[k][i] = malloc(sizeof(long long) * number_of_variables);
			for (j = 0; j < number_of_variables; j++)
			{
				(*binary_constraint_coefficient_matrix)[k][i][j] = 0;
			} // end of j loop
		} // end of i loop
	} // end of k loop
} // end of CreateBinaryConstraintArrays function

/* This function initializes the binary constraint coefficient matrix from the input constraint coefficient matrix based on inequality direction and reversed sign objective function coefficients */
void InitializeBinaryConstraintCoefficientMatrix(unsigned int number_of_constraints, unsigned int number_of_variables, int *constraint_inequality_direction_vector, long long **constraint_constant_vector, long long ***constraint_coefficient_matrix, unsigned int *sorted_objective_function_coefficient_vector_index, int *binary_objective_function_coefficient_reversed_sign_vector, long long **binary_constraint_constant_vector, long long ***binary_constraint_coefficient_matrix)
{
	unsigned int i, j, k;

	k = 0;
	for (i = 0; i < number_of_constraints; i++)
	{
		if (constraint_inequality_direction_vector[i] == 1) // if less than or equal to
		{
			k = InitializeBinaryConstraintCoefficientMatrixLessThan(number_of_variables, constraint_constant_vector, constraint_coefficient_matrix, sorted_objective_function_coefficient_vector_index, binary_objective_function_coefficient_reversed_sign_vector, binary_constraint_constant_vector, binary_constraint_coefficient_matrix, i, k);
		}
		else if (constraint_inequality_direction_vector[i] == 0) // if equal to
		{
			k = InitializeBinaryConstraintCoefficientMatrixLessThan(number_of_variables, constraint_constant_vector, constraint_coefficient_matrix, sorted_objective_function_coefficient_vector_index, binary_objective_function_coefficient_reversed_sign_vector, binary_constraint_constant_vector, binary_constraint_coefficient_matrix, i, k);

			k = InitializeBinaryConstraintCoefficientMatrixGreaterThan(number_of_variables, constraint_constant_vector, constraint_coefficient_matrix, sorted_objective_function_coefficient_vector_index, binary_objective_function_coefficient_reversed_sign_vector, binary_constraint_constant_vector, binary_constraint_coefficient_matrix, i, k);
		}
		else if (constraint_inequality_direction_vector[i] == -1) // if greater than or equal to
		{
			k = InitializeBinaryConstraintCoefficientMatrixGreaterThan(number_of_variables, constraint_constant_vector, constraint_coefficient_matrix, sorted_objective_function_coefficient_vector_index, binary_objective_function_coefficient_reversed_sign_vector, binary_constraint_constant_vector, binary_constraint_coefficient_matrix, i, k);
		}
	} // end of i loop
} // end of InitializeBinaryConstraintCoefficientMatrix function

/* This function initializes the binary constraint coefficient matrix from the input constraint coefficient matrix based on less than inequality direction and reversed sign objective function coefficients */
unsigned int InitializeBinaryConstraintCoefficientMatrixLessThan(unsigned int number_of_variables, long long **constraint_constant_vector, long long ***constraint_coefficient_matrix, unsigned int *sorted_objective_function_coefficient_vector_index, int *binary_objective_function_coefficient_reversed_sign_vector, long long **binary_constraint_constant_vector, long long ***binary_constraint_coefficient_matrix, unsigned int i, unsigned int k)
{
	unsigned int j;

	binary_constraint_constant_vector[0][k] = -constraint_constant_vector[0][i];
	binary_constraint_constant_vector[1][k] = constraint_constant_vector[1][i];

	for (j = 0; j < number_of_variables; j++)
	{
		if (binary_objective_function_coefficient_reversed_sign_vector[sorted_objective_function_coefficient_vector_index[j]] == 0)
		{
			binary_constraint_coefficient_matrix[0][k][j] = -constraint_coefficient_matrix[0][i][sorted_objective_function_coefficient_vector_index[j]];
			binary_constraint_coefficient_matrix[1][k][j] = constraint_coefficient_matrix[1][i][sorted_objective_function_coefficient_vector_index[j]];
		}
		else
		{
			binary_constraint_coefficient_matrix[0][k][j] = constraint_coefficient_matrix[0][i][sorted_objective_function_coefficient_vector_index[j]];
			binary_constraint_coefficient_matrix[1][k][j] = constraint_coefficient_matrix[1][i][sorted_objective_function_coefficient_vector_index[j]];

			LongLongRationalAddition(binary_constraint_constant_vector[0][k], binary_constraint_constant_vector[1][k], constraint_coefficient_matrix[0][i][sorted_objective_function_coefficient_vector_index[j]], constraint_coefficient_matrix[1][i][sorted_objective_function_coefficient_vector_index[j]], &binary_constraint_constant_vector[0][k], &binary_constraint_constant_vector[1][k]);
		}
	} // end of j loop

	return k + 1;
} // end of InitializeBinaryConstraintCoefficientMatrixLessThan function

/* This function initializes the binary constraint coefficient matrix from the input constraint coefficient matrix based on greater than inequality direction and reversed sign objective function coefficients */
unsigned int InitializeBinaryConstraintCoefficientMatrixGreaterThan(unsigned int number_of_variables, long long **constraint_constant_vector, long long ***constraint_coefficient_matrix, unsigned int *sorted_objective_function_coefficient_vector_index, int *binary_objective_function_coefficient_reversed_sign_vector, long long **binary_constraint_constant_vector, long long ***binary_constraint_coefficient_matrix, unsigned int i, unsigned int k)
{
	unsigned int j;

	binary_constraint_constant_vector[0][k] = constraint_constant_vector[0][i];
	binary_constraint_constant_vector[1][k] = constraint_constant_vector[1][i];

	for (j = 0; j < number_of_variables; j++)
	{
		if (binary_objective_function_coefficient_reversed_sign_vector[sorted_objective_function_coefficient_vector_index[j]] == 0)
		{
			binary_constraint_coefficient_matrix[0][k][j] = constraint_coefficient_matrix[0][i][sorted_objective_function_coefficient_vector_index[j]];
			binary_constraint_coefficient_matrix[1][k][j] = constraint_coefficient_matrix[1][i][sorted_objective_function_coefficient_vector_index[j]];
		}
		else
		{
			binary_constraint_coefficient_matrix[0][k][j] = -constraint_coefficient_matrix[0][i][sorted_objective_function_coefficient_vector_index[j]];
			binary_constraint_coefficient_matrix[1][k][j] = constraint_coefficient_matrix[1][i][sorted_objective_function_coefficient_vector_index[j]];

			LongLongRationalAddition(binary_constraint_constant_vector[0][k], binary_constraint_constant_vector[1][k], -constraint_coefficient_matrix[0][i][sorted_objective_function_coefficient_vector_index[j]], constraint_coefficient_matrix[1][i][sorted_objective_function_coefficient_vector_index[j]], &binary_constraint_constant_vector[0][k], &binary_constraint_constant_vector[1][k]);
		}
	} // end of j loop

	return k + 1;
} // end of InitializeBinaryConstraintCoefficientMatrixLessThan function

/* This function initiates branch and bound for binary variables */
void BalasAdditiveBranchAndBoundBinary(unsigned int number_of_binary_constraints, unsigned int number_of_variables, long long **binary_objective_function_coefficient_vector, long long **binary_constraint_constant_vector, long long ***binary_constraint_coefficient_matrix, int *current_binary_feasible_solution, int *best_binary_feasible_solution, double *best_binary_optimal_value_double, long long *best_binary_optimal_value)
{
	/* Start with first coefficient being 1*/
	current_binary_feasible_solution[0] = 1;

	BalasAdditiveBranchAndBoundBinaryEvaluation(number_of_binary_constraints, number_of_variables, binary_objective_function_coefficient_vector, binary_constraint_constant_vector, binary_constraint_coefficient_matrix, current_binary_feasible_solution, best_binary_feasible_solution, best_binary_optimal_value_double, best_binary_optimal_value, 0);

	/* Start with first coefficient being 0*/
	current_binary_feasible_solution[0] = 0;

	BalasAdditiveBranchAndBoundBinaryEvaluation(number_of_binary_constraints, number_of_variables, binary_objective_function_coefficient_vector, binary_constraint_constant_vector, binary_constraint_coefficient_matrix, current_binary_feasible_solution, best_binary_feasible_solution, best_binary_optimal_value_double, best_binary_optimal_value, 0);
} // end of BalasAdditiveBranchAndBoundBinary function

// This function evaluates branch and bound for binary variables
void BalasAdditiveBranchAndBoundBinaryEvaluation(unsigned int number_of_binary_constraints, unsigned int number_of_variables, long long **binary_objective_function_coefficient_vector, long long **binary_constraint_constant_vector, long long ***binary_constraint_coefficient_matrix, int *current_binary_feasible_solution, int *best_binary_feasible_solution, double *best_binary_optimal_value_double, long long *best_binary_optimal_value, unsigned int variable_index)
{
	unsigned int i, j;
	int possible = 1, feasible = 1;
	double z = 0, constraint_sum = 0, constraint_max = 0;

	/* Evaluate z */
	z = 0;
	for (i = 0; i < number_of_variables; i++)
	{
		z += ((double)binary_objective_function_coefficient_vector[0][i] / binary_objective_function_coefficient_vector[1][i] * current_binary_feasible_solution[i]);
	} // end of i loop

	/* Check if possible & feasible with constraints */
	possible = 1;
	feasible = 1;
	for (i = 0; i < number_of_binary_constraints; i++)
	{
		constraint_sum = 0;
		for (j = 0; j < variable_index + 1; j++)
		{
			constraint_sum += ((double)binary_constraint_coefficient_matrix[0][i][j] / binary_constraint_coefficient_matrix[1][i][j] * current_binary_feasible_solution[j]);
		} // end of j loop

		constraint_max = constraint_sum;
		for (j = variable_index + 1; j < number_of_variables; j++)
		{
			if (binary_constraint_coefficient_matrix[0][i][j] > 0)
			{
				constraint_max += ((double)binary_constraint_coefficient_matrix[0][i][j] / binary_constraint_coefficient_matrix[1][i][j]);
			}
		} // end of j loop

		if (constraint_max < (double)binary_constraint_constant_vector[0][i] / binary_constraint_constant_vector[1][i])
		{
			possible = 0;
			break; // break i loop since problem is already impossible
		}
		else
		{
			for (j = variable_index + 1; j < number_of_variables; j++)
			{
				constraint_sum += ((double)binary_constraint_coefficient_matrix[0][i][j] / binary_constraint_coefficient_matrix[1][i][j] * current_binary_feasible_solution[j]);
			} // end of j loop

			if (constraint_sum < (double)binary_constraint_constant_vector[0][i] / binary_constraint_constant_vector[1][i])
			{
				feasible = 0;
				break; // break i loop since problem is currently infeasible
			}
		}
	} // end of i loop

	if (possible == 1)
	{
		if (feasible == 1)
		{
			/* Keep best solution */
			if (z < (*best_binary_optimal_value_double))
			{
				(*best_binary_optimal_value_double) = z;

				best_binary_optimal_value[0] = 0;
				best_binary_optimal_value[1] = 1;

				for (i = 0; i < number_of_variables; i++)
				{
					LongLongRationalAddition(best_binary_optimal_value[0], best_binary_optimal_value[1], binary_objective_function_coefficient_vector[0][i] * current_binary_feasible_solution[i], binary_objective_function_coefficient_vector[1][i], &best_binary_optimal_value[0], &best_binary_optimal_value[1]);

					best_binary_feasible_solution[i] = current_binary_feasible_solution[i];
				} // end of i loop
			}
			/* else: Don't continue down rabbit hole, so no more recursion */
		}
		else // if NOT feasible
		{
			if (variable_index + 1 < number_of_variables)
			{
				/* Down and down we go! */
				BalasAdditiveBranchAndBoundBinaryRecursive(number_of_binary_constraints, number_of_variables, binary_objective_function_coefficient_vector, binary_constraint_constant_vector, binary_constraint_coefficient_matrix, current_binary_feasible_solution, best_binary_feasible_solution, best_binary_optimal_value_double, best_binary_optimal_value, variable_index + 1);
			}
			/* else: Don't continue down rabbit hole, so no more recursion */
		}
	}
	/* else: NOT possible, don't continue down rabbit hole, so no more recursion */
} // end of BalasAdditiveBranchAndBoundBinaryEvaluation

/* This function performs branch and bound recursively for binary variables*/
void BalasAdditiveBranchAndBoundBinaryRecursive(unsigned int number_of_binary_constraints, unsigned int number_of_variables, long long **binary_objective_function_coefficient_vector, long long **binary_constraint_constant_vector, long long ***binary_constraint_coefficient_matrix, int *current_binary_feasible_solution, int *best_binary_feasible_solution, double *best_binary_optimal_value_double, long long *best_binary_optimal_value, unsigned int variable_index)
{
	/* Start with variable coefficient being 1*/
	current_binary_feasible_solution[variable_index] = 1;

	BalasAdditiveBranchAndBoundBinaryEvaluation(number_of_binary_constraints, number_of_variables, binary_objective_function_coefficient_vector, binary_constraint_constant_vector, binary_constraint_coefficient_matrix, current_binary_feasible_solution, best_binary_feasible_solution, best_binary_optimal_value_double, best_binary_optimal_value, variable_index);

	/* Start with variable coefficient being 0*/
	current_binary_feasible_solution[variable_index] = 0;

	BalasAdditiveBranchAndBoundBinaryEvaluation(number_of_binary_constraints, number_of_variables, binary_objective_function_coefficient_vector, binary_constraint_constant_vector, binary_constraint_coefficient_matrix, current_binary_feasible_solution, best_binary_feasible_solution, best_binary_optimal_value_double, best_binary_optimal_value, variable_index);
} // end of BalasAdditiveBranchAndBoundBinaryRecursive function

/* This function writes to disk the final BILP array values */
void WriteFinalBILPValues(unsigned int number_of_variables, long long *optimal_objective_function_value, long long **optimal_variable_values)
{
	unsigned int i;

	FILE *outfile_final_BILP_optimal_objective_function_value = fopen("outputs/final_BILP_optimal_objective_function_value.txt", "w"); // write only
	fprintf(outfile_final_BILP_optimal_objective_function_value, "%lf\n", (double)optimal_objective_function_value[0] / optimal_objective_function_value[1]);
	fclose(outfile_final_BILP_optimal_objective_function_value);

	FILE *outfile_final_BILP_optimal_variable_values = fopen("outputs/final_BILP_optimal_variable_values.txt", "w"); // write only
	for (i = 0; i < number_of_variables; i++)
	{
		fprintf(outfile_final_BILP_optimal_variable_values, "%lf\n", (double)optimal_variable_values[0][i] / optimal_variable_values[1][i]);
	} // end of i loop
	fclose(outfile_final_BILP_optimal_variable_values);
} // end of WriteFinalBILPValues function

/*********************************************************************************/
/******************************* HELPER FUNCTIONS ********************************/
/*********************************************************************************/

/* This function prints the optimal objective function value and the variable values */
void PrintOptimalResults(unsigned int number_of_variables, double optimal_objective_function_value, long long **optimal_variable_values)
{
	unsigned int i;

	printf("PrintOptimalResults: optimal_objective_function_value = %.16f\n", optimal_objective_function_value);
	printf("PrintOptimalResults: optimal_variable_values:\n");
	printf("variable_index\tvariable_value\n");
	for (i = 0; i < number_of_variables; i++)
	{
		printf("%d\t%.16f\n", i, (double)optimal_variable_values[0][i] / optimal_variable_values[1][i]);
	} // end of i loop
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
} // end of Realloc1DUnsignedInt function

/* This function reallocates more memory to the passed 1d array */
void Realloc1DLongLong(long long ***array1d, unsigned int oldsize1d, unsigned int newsize1d)
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
			temp[k][i] = (*array1d)[k][i];
		} // end of i loop
	} // end of k loop

	/* Free 1d array */
	for (k = 0; k < 2; k++)
	{
		free((*array1d)[k]);
	} // end of k loop
	free(*array1d);

	(*array1d) = malloc(sizeof(long long*) * 2);
	for (k = 0; k < 2; k++)
	{
		(*array1d)[k] = malloc(sizeof(long long) * newsize1d);
	} // end of k loop

	for (k = 0; k < 2; k++)
	{
		for (i = 0; i < oldsize1d; i++)
		{
			(*array1d)[k][i] = temp[k][i];
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
		(*array1d)[0][i] = 0;
		(*array1d)[1][i] = 1;
	} // end of i loop
} // end of Realloc1DLongLong function

/* This function reallocates more memory to the passed 2d array */
void Realloc2DLongLong(long long ****array2d, unsigned int oldsize1d, unsigned int newsize1d, unsigned int oldsize2d, unsigned int newsize2d)
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
					temp[k][i][j] = (*array2d)[k][i][j];
				} // end of j loop
			} // end of i loop
		} // end of k loop

		/* Free 2d array */
		for (k = 0; k < 2; k++)
		{
			for (i = 0; i < oldsize1d; i++)
			{
				free((*array2d)[k][i]);
			} // end of i loop
			free((*array2d)[k]);
		} // end of k loop
		free(*array2d);

		(*array2d) = malloc(sizeof(long long*) * 2);
		for (k = 0; k < 2; k++)
		{
			(*array2d)[k] = malloc(sizeof(long long*) * newsize1d);
			for (i = 0; i < newsize1d; i++)
			{
				(*array2d)[k][i] = malloc(sizeof(double) * newsize2d);
			} // end of i loop
		} // end of k loop

		for (k = 0; k < 2; k++)
		{
			for (i = 0; i < oldsize1d; i++)
			{
				for (j = 0; j < oldsize2d; j++)
				{
					(*array2d)[k][i][j] = temp[k][i][j];
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
				(*array2d)[0][i][j] = 0;
				(*array2d)[1][i][j] = 1;
			} // end of j loop
		} // end of i loop
	} // end of if only 1st dim's size has changed or both

	if ((newsize1d == oldsize1d && newsize2d != oldsize2d) || (newsize1d != oldsize1d && newsize2d != oldsize2d)) // if only 2nd dim's size has changed or both
	{
		for (i = 0; i < newsize1d; i++) // for all rows
		{
			for (j = oldsize2d; j < newsize2d; j++) // fill in just new columns
			{
				(*array2d)[0][i][j] = 0;
				(*array2d)[1][i][j] = 1;
			} // end of j loop
		} // end of i loop
	} // end of if only 2nd dim's size has changed or both
} // end of Realloc2DLongLong function

/* This function adds long long rationals and keeps their numerators and denominators each as small as possible */
void LongLongRationalAddition(long long A_numerator, long long A_denominator, long long B_numerator, long long B_denominator, long long *C_numerator, long long *C_denominator)
{
	long long temp_gcd;

	(*C_numerator) = A_numerator * B_denominator + A_denominator * B_numerator;
	(*C_denominator) = A_denominator * B_denominator;
	temp_gcd = GreastestCommonDenominator((*C_numerator), (*C_denominator));
	(*C_numerator) /= temp_gcd;
	(*C_denominator) /= temp_gcd;
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
} // end of LongLongRationalDivision function

/* Find the greatest common denominator between a and b */
long long GreastestCommonDenominator(long long a, long long b)
{
	long long t;
	while (b != 0)
	{
		t = b;
		b = a % b;
		a = t;
	}
	return abs(a);
} // end of GreastestCommonDenominator function
