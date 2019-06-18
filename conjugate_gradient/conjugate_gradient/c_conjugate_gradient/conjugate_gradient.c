/* This function uses the iterative conjugate gradient method that solves A * x = b */
void ConjugateGradientMethod(unsigned int n, double** A, double** b, double** x)
{
	/* A = [n, n] */
	/* b = [n, 1] */
	/* x = [n, 1] */
	
	unsigned int i, j;
	
	double** r;
	r = malloc(sizeof(double*) * n);
	for (i = 0; i < n; i++)
	{
		r[i] = malloc(sizeof(double) * 1);
	} // end of i loop
	
	MatrixMatrixMultiplication(n, 1, n, A, x, 0, 0, r);
	
	for (i = 0; i < n; i++)
	{
		r[i][0] = b[i][0] - r[i][0];
	} // end of i loop
	
	double** p;
	p = malloc(sizeof(double*) * n);
	for (i = 0; i < n; i++)
	{
		p[i] = malloc(sizeof(double) * 1);
	} // end of i loop
	
	double rsold = 0.0, rsnew = 0;
	rsold = VectorDotProductRank2(n, r, r, 1, 1);
	
	double** Ap;
	Ap = malloc(sizeof(double*) * n);
	for (i = 0; i < n; i++)
	{
		Ap[i] = malloc(sizeof(double) * 1);
	} // end of i loop
	
	double alpha;
	for (i = 0; i < n; i++)
	{
		MatrixMatrixMultiplication(n, 1, n, A, p, 0, 0, Ap);
		
		alpha = rsold / VectorDotProductRank2(n, p, Ap, 1, 1);
		
		for (j = 0; j < n; j++)
		{
			x[j][0] += alpha * p[j][0];
			r[j][0] -= alpha * Ap[j][0];
		} // end of j loop
		
		rsnew = VectorDotProductRank2(n, r, r, 1, 1);
		if (sqrt(rsnew) < 1e-10)
		{
			break;
		}

		for (j = 0; j < n; j++)
		{
			p[j][0] = r[j][0] + (rsnew / rsold) * p[j][0];
		} // end of j loop
		
        rsold = rsnew;
	} // end of i loop
	
	/* Free dynamic memory */
	for (i = 0; i < n; i++)
	{
		free(Ap[i]);
		free(p[i]);
		free(r[i]);
	} // end of i loop
	free(Ap);
	free(p);
	free(r);

	return;
} // end of ConjugateGradientMethod function