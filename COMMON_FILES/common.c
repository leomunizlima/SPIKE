/*----------------------------------------------------------------------------
 * ALGEBRA LINEAR OPERATIONS
 *--------------------------------------------------------------------------*/
#include "protos.h"

/*----------------------------------------------------------------------------
 * Time measure
 *--------------------------------------------------------------------------*/
double get_TIME ()
{
	struct timeval tv; gettimeofday(&tv, 0);
	return (double)(tv.tv_sec * 100.0 + tv.tv_usec / 10000.0);
}

/*----------------------------------------------------------------------------
 * Get the extension of a file
 *--------------------------------------------------------------------------*/
char* get_EXTENSION (char *fspec)
{
	char *e = strrchr (fspec, '.');
	if (e == NULL)
		e = "";
	return e;
}

/*----------------------------------------------------------------------------
 * Read the input file (matrix in .mtx format)
 *--------------------------------------------------------------------------*/
void READ_INPUT (MAT* A, double** b, char* entry, int* msgs)
{

	/* Reading data from file: matrix A and RSH b */
	FILE *f = fopen(entry,"r");
	if (!f) { 
		fprintf (stderr,"Can't open file %s. Exiting.. [READ_INPUT]\n", entry);
		exit(0);
	}
	
	if (strcmp(get_EXTENSION(entry),".mtx") == 0) 
	{
		MATRIX_readCSR (A,f);
		
		int    i;
		int    n  = A->n;
		(*b)      = calloc(n,sizeof(double));

		char matrixname[300],filename[300];	
		FILE *fb;

		sscanf(entry,"%[^.].mtx",matrixname);
		sprintf(filename,"%s_b.mtx",matrixname); 

		fprintf(stderr,"%s\n",filename);

		fb = fopen(filename,"r");

		if (fb){ //if there is a file <matrix>_b.mtx with b vector data
			READ_RHS ((*b), n, fb);
			fclose(fb);			
		}			
		else{ // otherwise vector b is obtained in such way to be the answer vector x 
			double* x = calloc(n,sizeof(double));
			for (i = 0; i < n; ++i)
				x[i] = 1.0;			
			MATRIX_matvec (A,x,(*b));
			free(x);
		}
	}
	else
	{
		fprintf (stderr,"File needs to be .mtx format\n");
		exit(0);
	}
	
	fclose(f);
	return;
}

/*----------------------------------------------------------------------------
 * Read right-hand-side from file in .mtx format
 *--------------------------------------------------------------------------*/
void READ_RHS (double* b, int n, FILE* f)
{
	int M, N;
	int i;
	double VAL;
	char line[1025], header[1025],test[1025];
	char *head0, *head1, *rep, *field, *symm;
	
	fgets(header,200,f);
	head0 = strtok(header," \n");
	head1 = strtok(NULL," \n");
	rep   = strtok(NULL," \n");
	field = strtok(NULL," \n");
	symm  = strtok(NULL," \n");
	
	if ((strcmp(symm,"symmetric") != 0 && strcmp(symm,"general") != 0))
		fprintf (stderr,"\n Error. Matrix format not supported. Exiting.. [MATRIX_readCSR]\n\n");
	
	do 
	{
		if (fgets(line,1025,f) == NULL) 
			exit(0);
	}while (line[0] == '%');   
	sscanf(line,"%d %d", &N, &M);
	
	for(i=0; i < N; i++)
		fscanf(f,"%lf",&(b[i]));
}




/*----------------------------------------------------------------------------
 * Assembly and write the global vector solution x
 *--------------------------------------------------------------------------*/
void WRITE_SOLUTION (double* x, char* entry, char* output, int n, int* P, int* K, double* dw, int RANK, int* msgs)
{
	int i, j;


	if (RANK==0){	
		double *s = calloc(n,sizeof(double));
		/* Knapsack retrieve */
		if (msgs[4])
		{
			for (i = 0; i < n; ++i)
				s[K[i]] = x[i];
				
			memcpy(&x[0],&s[0],n*sizeof(double));
		}
	
		/* Reordering retrieve */
		if (msgs[3])
		{
			for (i = 0; i < n; ++i)
				s[P[i]] = x[i];
			
			memcpy(&x[0],&s[0],n*sizeof(double));
		}
				
		/* Scaling retrieve */
		if (msgs[1]){
			for (i = 0; i < n; ++i)
				x[i] = x[i] * exp(dw[i]);
		}
		char* r = malloc(100*sizeof(char));
		      r = strcat(output, ".txt");
		FILE *f = fopen(r, "w");
		
		for (i = 0; i < n; ++i)
			fprintf(f,"(%d) %.3lf\n",i,x[i]);
		fclose(f);
                
                double sum = 0;
 	       	double max = -1e100;
 	       	double min = 1e100;
                for (i = 0; i < n; ++i)
 		{
                        sum += x[i];
 			if (x[i]>max)
 				max = x[i];
 			if (x[i]<min)
 				min = x[i];
 		}
                
 		fprintf(stderr, "  - Average Solution: %lf\n",sum/n);
 		fprintf(stderr, "  - Max Solution: %lf\n",max);
 		fprintf(stderr, "  - Min Solution: %lf\n",min);
 		free(s);
 	
	}	
}

/*----------------------------------------------------------------------------
 * Constant times a vector plus a vector
 *--------------------------------------------------------------------------*/
inline int daxpy (int n, double a, double *x, double *y)
{
	long int i, m;
	register double sa;

	sa = a;
	m = n-3;
	for (i = 0; i < m; i += 4){
		y[i]   += sa * x[i];
		y[i+1] += sa * x[i+1];
		y[i+2] += sa * x[i+2];
		y[i+3] += sa * x[i+3];
	}
	for ( ; i < n; ++i) /* clean-up loop */
		y[i]   += sa * x[i];
	
	return 0;
}

/*----------------------------------------------------------------------------
 * Copies a vector x to a vector y
 *--------------------------------------------------------------------------*/
inline int dcopy (int n, double *x, double *y)
{
	long int i, m;

	m = n-3;
	for (i = 0; i < m; i += 4){
		y[i]   = x[i];
		y[i+1] = x[i+1];
		y[i+2] = x[i+2];
		y[i+3] = x[i+3];
	}
	for ( ; i < n; ++i) /* clean-up loop */
		y[i]   = x[i];
	
	return 0;
}

/*----------------------------------------------------------------------------
 * Forms the dot product of two vectors
 *--------------------------------------------------------------------------*/
inline double ddot (int n, double *x, double *y)
{
	long int i, m;
	double stemp;

	stemp = 0.0;
	m = n-4;

	for (i = 0; i < m; i += 5)
		stemp += x[i] * y[i] + x[i+1] * y[i+1] + x[i+2] * y[i+2] + x[i+3] * y[i+3] + x[i+4] * y[i+4];

	for ( ; i < n; i++)        /* clean-up loop */
		stemp += x[i] * y[i];

	return stemp;
} 

/*----------------------------------------------------------------------------
 * Scales a vector by a constant
 *--------------------------------------------------------------------------*/
inline int dscal (int n, double a, double *x, double* y)
{
	long int i, m;
	register double sa;

	sa = a;
	m = n-3;
	for (i = 0; i < m; i += 4){
		y[i]   = sa * x[i];
		y[i+1] = sa * x[i+1];
		y[i+2] = sa * x[i+2];
		y[i+3] = sa * x[i+3];
	}
	for ( ; i < n; ++i) /* clean-up loop */
		y[i]   = sa * x[i];

	return 0;
}

/*----------------------------------------------------------------------------
 * Erase a vector of doubles
 *--------------------------------------------------------------------------*/
inline int ddiff (int n, double *a, double *x, double *y)
{
	long int i, m;

	m = n-3;
	for (i = 0; i < m; i += 4){
		y[i]   = a[i]   - x[i];
		y[i+1] = a[i+1] - x[i+1];
		y[i+2] = a[i+2] - x[i+2];
		y[i+3] = a[i+3] - x[i+3];
	}
	for ( ; i < n; ++i) 
		y[i]   = a[i]   - x[i];;

	return 0;
}

/*----------------------------------------------------------------------------
 * Erase a vector of integers
 *--------------------------------------------------------------------------*/
inline int izero (int n, int *v)
{
	long int i, m;

	m = n-3;
	for (i = 0; i < m; i += 4){
		v[i]   = 0;
		v[i+1] = 0;
		v[i+2] = 0;
		v[i+3] = 0;
	}
	for ( ; i < n; ++i) 
		v[i]   = 0;

	return 0;
}

/*----------------------------------------------------------------------------
 * Test if an element is zero
 *--------------------------------------------------------------------------*/
int iszero (double x)
{
	return (x < 1e-15 && x > -(1e-15));
}

