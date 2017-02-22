/*----------------------------------------------------------------------------
 * MIN_MAX PARTITIONING
 *--------------------------------------------------------------------------*/
#include "../COMMON_FILES/protos.h"

# define KRED "\x1B[31m"
# define KGRN "\x1B[32m"
# define KBLU "\x1B[34m"
# define RESET "\033[0m"

/*----------------------------------------------------------------------------
 * MIN_MAX partitioning
 *--------------------------------------------------------------------------*/
void PARTITIONING_MIN_MAX (MAT* A, int nP, int** Fs, int** FL, int* msgs) 
{
	/* -------------------------------------------------------------------------------------
	 w: vector of n weights (nnz by line)
	 s: vector of k+1 separator (each partition i=0,k-1 has from s[i] to s[i+1]-1)
	 L: vector of k loads (total weight for each partition)
	--------------------------------------------------------------------------------------*/
	double time;
	/*---START TIME---------------> */ time =  get_TIME(); 
	
	int  i, j, B;
	int  n = A->n;	
	int* w = (int *) malloc(n*sizeof(int));
	int* s = (int *) malloc((nP+1)*sizeof(int));
	int* L = (int *) malloc(nP*sizeof(int));
		
	for (i = 0; i < n; ++i)
		w[i] = A->IA[i+1] - A->IA[i];
	
	for (i = 0; i < nP; ++i) 
		s[i] = i;
	
	s[nP] = n;
	for (i = 0; i < nP; ++i)
	{
		L[i] = 0;
		for (j = s[i]; j < s[i+1]; ++j)
			L[i] += w[j];
	}
	do {
		for (i = 0, j = 1; j < nP; ++j) 
			if (L[j] > L[i]) i = j;
		B = L[i];
		j = i;
		
		if (s[i+1] == s[i]+1) break;
		
		while (L[i] >= B && i > 0)
		{
			L[i]   -= w[s[i]];
			L[i-1] += w[s[i]];
			s[i]    = s[i]+1;
			if (L[i] < B) i = i-1;
		}
	} while (L[0] < B);
	while (i != j)
	{
		if (L[i+1] < B) i = i+1;
		s[i]    = s[i]-1;
		L[i]   += w[s[i]];
		L[i-1] -= w[s[i]];
	}
	
	(*Fs) = s;
	(*FL) = L;
	free(w); 
	
	/*---FINAL TIME---------------> */ time = (get_TIME() - time)/100.0;
	int outside = PARTITIONING_nnzout(A,s,nP);
		
	if (msgs[0])
	{
		fprintf(stderr,"\n  [ PARTITIONING into %d machines ]\n", nP);	
		for (i = 0; i < nP; ++i)
			fprintf(stderr,"  - Partition %03d: size = %d; rows = %04d:%04d; nnz = %d\n", \
			i+1, s[i+1]-s[i], s[i]+1, s[i+1], L[i]);
		fprintf(stderr,"  - Elements outside of blocks = %d\n",outside);
		fprintf(stderr,"  - Elapsed time: %.6f sec\n", time);
	}
	
	if (outside != 0)
	{
		fprintf(stderr, "\n Error. %d elements outside of blocks (must be zero). Exiting.. [PARTITIONING_MIN_MAX]\n\n",outside);
		exit(0);
	}
	
	
	return;
}

/*----------------------------------------------------------------------------
 * Verify if there are elements outside de diagonal or coupling blocks
 *--------------------------------------------------------------------------*/
int PARTITIONING_nnzout (MAT *A, int *s, int nP)
{
	int i, ii, j, c=0;
	for (i = s[0]; i < s[1]; ++i) 
	{ 
		for (j = A->IA[i]; j < A->IA[i+1]; ++j)
			if (A->JA[j]>=s[2]) c++; 
	}
	for (i = 1; i < nP-1; ++i)
	{
		for (ii = s[i]; ii < s[i+1]; ++ii)
		{
			for (j = A->IA[ii]; j < A->IA[ii+1]; ++j)
				if (A->JA[j] < s[i-1] || A->JA[j] >= s[i+2]) c++;
		}
	}
	for (i = s[nP-1]; i < s[nP]; ++i)
	{ 
		for (j = A->IA[i]; j < A->IA[i+1]; ++j) 
			if (A->JA[j]<s[nP-2]) c++;		
	}
	
	return c;
}

/*----------------------------------------------------------------------------
 * Print Matrix A in FULL (dense) format
 *--------------------------------------------------------------------------*/
void PARTITIONING_printFULL (MAT *A, int *s, int nP, int r)
{
	int i, j, ii, jj, c = 0;
	printf("\n");
	for (ii = s[0]; ii < s[1]; ++ii)
	{
		for (j = 0; j < nP; ++j) 
		{
			for (jj = s[j]; jj < s[j+1]; ++jj)
			{
				if (!j) 
					printf(KRED " %.0f" RESET, MATRIX_aij(A,ii,jj));
				else if (ii < s[1] && ii >= s[1]-r && jj >= s[1] && jj < s[1]+r) 
					printf(KBLU " %.0f" RESET, MATRIX_aij(A,ii,jj));
				else if (!iszero(MATRIX_aij(A,ii,jj))) 
				{ 
					printf(KGRN " %.0f" RESET, MATRIX_aij(A,ii,jj));
					c++;
				} else 
					printf(" %.0f", MATRIX_aij(A,ii,jj));
			}
		}
		printf("\n");
	}
	for (i = 1; i < nP-1; i++)
	{
		for (ii = s[i]; ii < s[i+1]; ++ii)
		{
			for (j = 0; j < nP ; ++j)
			{
				for (jj = s[j]; jj < s[j+1]; ++jj)
				{
					if (i == j) 
						printf(KRED " %.0f" RESET, MATRIX_aij(A,ii,jj));
					else if (ii >= s[i] && ii < s[i]+r && jj < s[i] && jj >= s[i]-r) 
						printf(KBLU " %.0f" RESET, MATRIX_aij(A,ii,jj));
					else if (ii < s[i+1] && ii >= s[i+1]-r && jj >= s[i+1] && jj<s[i+1]+r) 
						printf(KBLU " %.0f" RESET, MATRIX_aij(A,ii,jj));
					else if (!iszero(MATRIX_aij(A,ii,jj))) 
					{ 
						printf(KGRN " %.0f" RESET, MATRIX_aij(A,ii,jj));
						c++; 
					} else 
						printf(" %.0f", MATRIX_aij(A,ii,jj));
				}
			}
			printf("\n");
		}
	}
	for (ii = s[nP-1]; ii < s[nP]; ++ii)
	{
		for (j = 0; j < nP; ++j)
		{
			for (jj = s[j]; jj < s[j+1]; ++jj)
			{
				if (j == (nP-1)) 
					printf(KRED " %.0f" RESET, MATRIX_aij(A,ii,jj));
				else if (ii >= s[nP-1] && ii < s[nP-1]+r && jj < s[nP-1] && jj >= s[nP-1]-r) 
					printf(KBLU " %.0f" RESET, MATRIX_aij(A,ii,jj));
				else if (!iszero(MATRIX_aij(A,ii,jj))) 
				{ 
					printf(KGRN " %.0f" RESET, MATRIX_aij(A,ii,jj));
					c++; 
				} else 
					printf(" %.0f", MATRIX_aij(A,ii,jj));
			}
		}
		printf("\n");
	}
	printf("\n"); 
}
