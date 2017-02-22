/*----------------------------------------------------------------------------
 * MATRIX FUNCTIONS: IN CSR FORMAT
 *--------------------------------------------------------------------------*/
#include "protos.h"

int COMPARE_array (const void * a, const void * b)
{
	if (((ARRAY*)a)->arr3 <  ((ARRAY*)b)->arr3) return -1;
	if (((ARRAY*)a)->arr3 >  ((ARRAY*)b)->arr3) return  1;
	if (((ARRAY*)a)->arr3 == ((ARRAY*)b)->arr3)
	{
		if (((ARRAY*)a)->arr2 < ((ARRAY*)b)->arr2) return -1;
		if (((ARRAY*)a)->arr2 > ((ARRAY*)b)->arr2) return  1;
	}
	return 0;
}

/*----------------------------------------------------------------------------
 * Read matrix from file in MM format to a CSR structure
 *--------------------------------------------------------------------------*/
void MATRIX_readCSR (MAT* A, FILE* f)
{
	int M, N, nz, new_nz;
	int i, j, k, I, J, elem = 0;
	double VAL;
	char line[1025], header[1025],test[1025];
	char *head0, *head1, *rep, *field, *symm;
	ARRAY* a;
	
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
	sscanf(line,"%d %d %d", &N, &M, &nz);
	
	if (strcmp(symm,"symmetric") == 0)
	{
		a = malloc ((2*nz)*sizeof(ARRAY));
		for (i = 0, k = 0; i < nz; ++i)
		{
			fscanf(f,"%d%d%lf",&I,&J,&VAL);
			a[k].arr1 = VAL;
			a[k].arr2 = J - 1;
			a[k].arr3 = I - 1;
			k++;
			if (I != J)
			{
				a[k].arr1 = VAL;
				a[k].arr2 = I - 1;
				a[k].arr3 = J - 1;
				k++;
			}
		}
		nz = k;
		a  = realloc (a,nz*sizeof(ARRAY));
		qsort(a,nz,sizeof(ARRAY),COMPARE_array);
	}
	
	if (strcmp(symm,"general") == 0)
	{
		a = malloc (nz*sizeof(ARRAY));
		for (i = 0; i < nz; ++i)
		{
			fscanf(f,"%d%d%lf",&I,&J,&VAL);
			a[i].arr1 = VAL;
			a[i].arr2 = J - 1;
			a[i].arr3 = I - 1;
		}
		qsort(a,nz,sizeof(ARRAY),COMPARE_array);
	}

	/* reseve memory for matrices */
	A->m   = M;
	A->n   = N;
	A->nz  = nz;
	
	A->AA  = (double *) calloc(nz, sizeof (double));
	A->D   = (double *) calloc(N , sizeof (double));
	A->JA  = (int    *) calloc(nz, sizeof (int));
	A->IA  = (int    *) calloc(N+1,sizeof (int));
	
	for (i = 0; i < nz; ++i)
	{
		A->AA[i]   = a[i].arr1;
		A->JA[i]   = a[i].arr2;
		elem      += 1;
		A->IA[a[i].arr3+1] = elem;
	}

	free(a);
	
	/* Adjusting IA array */
	for (i = 1; i < N + 1; ++i)
		if (A->IA[i] == 0) 
			A->IA[i] = A->IA[i-1];

	/* Diagonal */
	if (M == N) /* square matrix */
	{
		for (i = 0; i < A->n; ++i)
		{
			int k1 = A->IA[i];
			int k2 = A->IA[i+1]-1;
			for (j = k1; j <= k2; ++j)
				if (i == A->JA[j]) 
					A->D[i] = A->AA[j];
		}
	}
}

/*----------------------------------------------------------------------------
 * Write matrix stored in CSR for parallel computations 
 *--------------------------------------------------------------------------*/
void MATRIX_writeCSR (MAT* A, double* f, int* s, int nP, int bandwidth, int *msgs)
{
	int  i,ii,jj,k,k1,k2,k3,k4,k5,nnz;
	char file[100];
	FILE *p;
	int nnz_A_bef, nnz_C, nnz_A_mid, nnz_A_aft, nnz_B;
	int sc, m, r, knp; 

	sc = msgs[1]; 
	m = msgs[2]; 
	r = msgs[3]; 
	knp = msgs[4]; 

	
	ARRAY *A_bef, *C, *A_mid, *A_aft, *B;
	nnz_C	  = 0;
	nnz_A_mid = 0;
	nnz_B     = 0;
	
	for (i = 0; i < nP; ++i)
	{	
		
		for (ii = s[i]; ii < s[i+1]; ++ii)
		{
			for (jj = A->IA[ii]; jj <= A->IA[ii+1] - 1; ++jj)
			{
				if (i != 0 && A->JA[jj] >= s[i-1] && A->JA[jj] < s[i])			// A_bef
				{
					if (ii < s[i] + bandwidth && A->JA[jj] >= s[i] - bandwidth)				// C
						nnz_C++;
				}
				else if (A->JA[jj] >= s[i] && A->JA[jj] < s[i+1]) 			// A
					nnz_A_mid++;
				else if (i != nP - 1 && A->JA[jj] >= s[i+1] && A->JA[jj] < s[i+2])	// A_aft
				{	
					if (ii >= s[i+1] - bandwidth && A->JA[jj] < s[i+1] + bandwidth)				// B
						nnz_B++;
				}
					
			}
			
		}

	}

	for (i = 0; i < nP; ++i)
	{
		nnz_A_bef = 0;
		nnz_C     = 0;
		nnz_A_mid = 0;
		nnz_A_aft = 0;
		nnz_B     = 0;
		
		for (ii = s[i]; ii < s[i+1]; ++ii)
		{
			for (jj = A->IA[ii]; jj <= A->IA[ii+1] - 1; ++jj)
			{
				if (i != 0 && A->JA[jj] >= s[i-1] && A->JA[jj] < s[i])			// A_bef
				{
					nnz_A_bef++;
					if (ii < s[i] + bandwidth && A->JA[jj] >= s[i] - bandwidth)				// C
						nnz_C++;
				}
				else if (A->JA[jj] >= s[i] && A->JA[jj] < s[i+1]) 			// A
					nnz_A_mid++;
				else if (i != nP - 1 && A->JA[jj] >= s[i+1] && A->JA[jj] < s[i+2])	// A_aft
				{	
					nnz_A_aft++;
					if (ii >= s[i+1] - bandwidth && A->JA[jj] < s[i+1] + bandwidth)				// B
						nnz_B++;
				}
					
			}
			
		}

		A_mid = calloc(nnz_A_mid,sizeof(ARRAY));

		if (i != 0)
		{
			A_bef = calloc(nnz_A_bef,sizeof(ARRAY));
			if (nnz_C==0){
				nnz_C = 1;
				C     = calloc(nnz_C    ,sizeof(ARRAY));
				C[0].arr2 = 0;
				C[0].arr3 = bandwidth-1;
			}
			else
				C    = calloc(nnz_C    ,sizeof(ARRAY));
		}
		if (i != nP-1) 
		{
			A_aft = calloc(nnz_A_aft,sizeof(ARRAY));
			if (nnz_B == 0){
				nnz_B = 1;
				B     = calloc(nnz_B    ,sizeof(ARRAY));
				B[0].arr2 = s[i+1]-s[i]-1;
				B[0].arr3 = 0;
			}
			else
				B     = calloc(nnz_B    ,sizeof(ARRAY));


		}

		k1 = k2 = k3 = k4 = k5 = 0;
		for (ii = s[i]; ii < s[i+1]; ++ii)
		{
			
			for (jj = A->IA[ii]; jj <= A->IA[ii+1] - 1; ++jj)
			{
				
				if (i != 0 && A->JA[jj] >= s[i-1] && A->JA[jj] < s[i])			// A_bef
				{ 	
					A_bef[k1].arr1 = A->AA[jj];
					A_bef[k1].arr2 = ii - s[i];
					A_bef[k1].arr3 = A->JA[jj]-s[i-1];
					k1++;
					if (ii < s[i] + bandwidth && A->JA[jj] >= s[i] - bandwidth)				// C
					{
						C[k2].arr1 = A->AA[jj];
						C[k2].arr2 = ii - s[i];
						C[k2].arr3 = A->JA[jj]-s[i]+bandwidth;
						k2++;
					}
				}
	
				else if (A->JA[jj] >= s[i] && A->JA[jj] < s[i+1]) 			// A
				{	
					A_mid[k3].arr1 = A->AA[jj];
					A_mid[k3].arr2 = ii - s[i];
					A_mid[k3].arr3 = A->JA[jj]-s[i];
					k3++;
				}
				
				else if (i != nP - 1 && A->JA[jj] >= s[i+1] && A->JA[jj] < s[i+2])	// A_aft
				{	
					A_aft[k4].arr1 = A->AA[jj];
					A_aft[k4].arr2 = ii - s[i];
					A_aft[k4].arr3 = A->JA[jj]-s[i+1];
					k4++;
					if (ii >= s[i+1] - bandwidth && A->JA[jj] < s[i+1] + bandwidth)				// B
					{
						B[k5].arr1 = A->AA[jj];
						B[k5].arr2 = ii - s[i];
						B[k5].arr3 = A->JA[jj]-s[i+1];
						k5++;
					}
				}
					
			}
			
		}

		int PRECISION = 14;
		sprintf(file,"data0%d_%d_%d_%d_%d.mtx",i, sc, m, r, knp);
		p = fopen(file,"w");
				
		if (i != 0)
		{
			fprintf (p, "%%%%MatrixMarket matrix coordinate real general\n");
			fprintf (p, "%d %d %d\n", s[i+1] - s[i], s[i] - s[i-1], nnz_A_bef);
			for (k = 0; k < nnz_A_bef; ++k)
				fprintf (p, "%d %d %.*e\n", A_bef[k].arr2+1, A_bef[k].arr3+1, PRECISION, A_bef[k].arr1); 
			fprintf (p,"\n");
			fprintf (p, "%%%%MatrixMarket matrix coordinate real general\n");
			fprintf (p, "%d %d %d\n", s[i+1] - s[i], bandwidth, nnz_C);
			for (k = 0; k < nnz_C; ++k){
				fprintf (p, "%d %d %.*e\n", C[k].arr2+1, C[k].arr3+1, PRECISION, C[k].arr1); 
			}
			fprintf (p,"\n");
		}
		
		fprintf (p, "%%%%MatrixMarket matrix coordinate real general\n");
		fprintf (p, "%d %d %d\n", s[i+1] - s[i], s[i+1] - s[i], nnz_A_mid);
		for (k = 0; k < nnz_A_mid; ++k){
			fprintf (p, "%d %d %.*e\n", A_mid[k].arr2+1, A_mid[k].arr3+1,PRECISION, A_mid[k].arr1); 
		}
		fprintf (p,"\n");
		
		if (i != nP - 1)
		{
			fprintf (p, "%%%%MatrixMarket matrix coordinate real general\n");
			fprintf (p, "%d %d %d\n", s[i+1] - s[i], s[i+2] - s[i+1], nnz_A_aft);
			for (k = 0; k < nnz_A_aft; ++k)
				fprintf (p, "%d %d %le\n", A_aft[k].arr2+1, A_aft[k].arr3+1, A_aft[k].arr1); 
			fprintf (p,"\n");
			fprintf (p, "%%%%MatrixMarket matrix coordinate real general\n");
			fprintf (p, "%d %d %d\n", s[i+1] - s[i], bandwidth, nnz_B);
			for (k = 0; k < nnz_B; ++k){
				fprintf (p, "%d %d %.*e\n", B[k].arr2+1, B[k].arr3+1, PRECISION, B[k].arr1);
			}	
			fprintf (p,"\n");
		}
		
		for (k = s[i]; k < s[i+1]; ++k)
			fprintf (p,"%.*e\n",PRECISION, f[k]);
		fprintf (p,"\n");
		
		fclose(p);
		
		free(A_mid);
		if (i != 0)
		{
			free(A_bef);
			free(C);
		}
		if (i != nP - 1)
		{
			free(A_aft);
			free(B);
		}
	}
}

long int MATRIX_bandwidth (MAT* A)
{
	int i;
	long int bandl, bandr;
	int n = A->n;
	unsigned long int bandwidth=0;
	
	bandl = 0;
	bandr = 0;
	for (i = 0; i < n ; i++)
	{
		if (i - A->JA[A->IA[i]] > bandl)
			bandl= i - A->JA[A->IA[i]];
		if (A->JA[A->IA[i+1]-1]-i > bandr);	
			bandr = A->JA[A->IA[i+1]-1]-i; 
	}
	if (bandl>bandr)
		bandwidth = bandl;
	else
		bandwidth = bandr;

	return bandwidth;
}



/*----------------------------------------------------------------------------
 * Compute the matrix-vector operation Ax = b
 *--------------------------------------------------------------------------*/
void MATRIX_matvec (MAT* A, double* x, double* b)
{
	int i, j, k1, k2;
	int n = A->n;
	double soma = 0;
	
	for (i = 0; i < n; ++i)
	{
		soma = 0;
		k1 = A->IA[i];
		k2 = A->IA[i+1]-1;
		for (j = k1; j <= k2; ++j)
			soma = soma + (A->AA[j]) * x[A->JA[j]];
		b[i] = soma;
	}  
}


/*----------------------------------------------------------------------------
 * Perform the colunm permutation
 *--------------------------------------------------------------------------*/
void MATRIX_COL_permutation (MAT* A, int* p)
{
	int i, j, k;
	int n   = A->n;
	int m   = A->m;
	int nz  = A->nz;  

	ARRAY* a = calloc (nz,sizeof(ARRAY));
	int*   q = calloc (n ,sizeof(int));
	
	for (i = 0; i < n; ++i) 
		q[p[i]] = i; 

	k = 0;
	for (i = 0; i < n; ++i)
	{
		for (j = A->IA[i]; j <= A->IA[i+1] - 1; ++j)
		{
			a[k].arr1 = A->AA[j];
			a[k].arr2 = q[A->JA[j]];
			a[k].arr3 = i;
				k = k + 1;
		}
		A->IA[i+1] = k;    
	}

	qsort(a,nz,sizeof(ARRAY),COMPARE_array);
	
	for (i = 0; i < nz; ++i)
	{
		A->AA[i] = a[i].arr1;
		A->JA[i] = a[i].arr2;
	}

	free(a);
	free(q);
}

/*----------------------------------------------------------------------------
 * Perform the colunm permutation
 *--------------------------------------------------------------------------*/
void MATRIX_ROW_permutation (MAT* A, int* p)
{
	int i, j, k;
	int n   = A->n;
	int m   = A->m;
	int nz  = A->nz;  

	double* auxAA = malloc( nz  *sizeof(double));
	int*    auxJA = malloc( nz  *sizeof(int));
	int*    auxIA = malloc((n+1)*sizeof(int));
  
	auxIA[0] = 0;
	k = 0;
	for (i = 0; i < n; ++i)
	{
		for (j = A->IA[p[i]]; j <= A->IA[p[i]+1] - 1; ++j)
		{
			auxAA[k] = A->AA[j];
			auxJA[k] = A->JA[j];
			      k  = k + 1;
		}
		auxIA[i+1] = k;    
	}

	memcpy(&A->AA[0],&auxAA[0],nz*sizeof(double));
	memcpy(&A->JA[0],&auxJA[0],nz*sizeof(int));
	memcpy(&A->IA[0],&auxIA[0],(n+1)*sizeof(int));

	free(auxAA);
	free(auxJA);
	free(auxIA);
}

/*----------------------------------------------------------------------------
 * Get the element ij stored as CSR
 *--------------------------------------------------------------------------*/
double MATRIX_aij (MAT* A, int i, int j)
{
	int k;
	int k1 = A->IA[i];
	int k2 = A->IA[i+1]-1;
	for (k = k1; k <= k2; ++k)
		if (A->JA[k] == j)
			return A->AA[k];
	return 0;
}


/*----------------------------------------------------------------------------
 *  Convert a matrix in CSR format to full format
 *--------------------------------------------------------------------------*/
void MATRIX_CSRto_FULL (MAT* A, double* a)
{
	int i, j, k1, k2;
	
	int n = A->n;
	int m = A->m;
	
	for (i = 0; i < n; ++i)
	{
		k1 = A->IA[i];
		k2 = A->IA[i+1]-1;
		for (j = k1; j <= k2; ++j)
			a[i + A->JA[j]*n] = A->AA[j];
	}
}

/*----------------------------------------------------------------------------
 * Free up memory allocated for MAT struct
 *--------------------------------------------------------------------------*/
void MATRIX_clean (MAT* A)
{
	free(A->AA);
	free(A->JA);
	free(A->IA);
	free(A->D);
	free(A);  
}


