/*----------------------------------------------------------------------------
 * ESPECTRAL REORDERING SOLVER
 *--------------------------------------------------------------------------*/
#include "../COMMON_FILES/protos.h"

/*----------------------------------------------------------------------------
 * SPECTRAL reordering
 *--------------------------------------------------------------------------*/
void mc73_fiedler(int* n, int* lirn, int* irn, int* ip, int* list, double* fvector, int* info, double* a);

void REORDERING_SPECTRAL (MAT* A, double* b, int** P, int* msgs)
{
	double time;
	long int bandwidth;
	
	bandwidth = MATRIX_bandwidth(A);
	if (msgs[0]) fprintf(stderr,"\n  [ REORDERING SPECTRAL]\n");
	if (msgs[0]) fprintf(stderr,"  - Bandwidth before = %ld\n", bandwidth);
	
	/*---START TIME---------------> */ time =  get_TIME();
	int i;
	int n = A->n;
	
	double *fvector = calloc (n ,sizeof(double));
	int    *list    = calloc (n ,sizeof(int));
	int    *info    = calloc (10,sizeof(int));
	int    *p       = calloc (n ,sizeof(int));
	int    *ia = A->IA;
	int    *ja = A->JA;
	int   lirn = A->nz;
	int    nnz = A->nz;
	
	/* -------------------------------------------------------------------- */    
	/* Convert matrix from 0-based C-notation to Fortran 1-based notation   */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < nnz; i++) {
		ja[i] += 1;
	}
	for (i = 0; i < n + 1; i++) {
		ia[i] += 1;
	}

	mc73_fiedler(&n,&lirn,ja,ia,list,fvector,info,NULL);
	
	/* -------------------------------------------------------------------- */    
	/* Convert matrix back to 0-based C-notation.                           */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < nnz; i++) {
		ja[i] -= 1;
	}
	for (i = 0; i < n + 1; i++) {
		ia[i] -= 1;
	}
	
	ARRAY* R = malloc(n*sizeof(ARRAY));
	for (i = 0; i < n; ++i)
	{
		R[i].arr1 = fvector[i];
		R[i].arr2 = i;		
	}	
	
	qsort (R,n,sizeof(ARRAY),COMPARE_eig); 
	
	for (i = 0; i < n; ++i) 
		p[i] = R[i].arr2; 

	free(R);
	free(info);
	free(list);
	free(fvector);
	
	MATRIX_ROW_permutation(A,p);
	MATRIX_COL_permutation(A,p);
	
	/*---PERMUTE COEF. VECTOR------ */
	double* bp = calloc(A->n,sizeof(double));		
	for (i = 0; i < A->n; ++i) 
		bp[i] = b[p[i]];
	
	memcpy(&b[0],&bp[0],A->n*sizeof(double));
	free(bp);
	
	/*---FINAL TIME---------------> */ time = (get_TIME() - time)/100.0;
	
	bandwidth = (int) MATRIX_bandwidth(A);
	if (msgs[0]) fprintf(stderr,"  - Bandwidth after  = %ld\n", bandwidth);
	if (msgs[0]) fprintf(stderr,"  - Elapsed time: %.6f sec\n", time);	

	(*P) = p;
	
	return;
}

void REORDERING_SPECTRAL_WGT (MAT* A, double* b, int** P, int* msgs)
{
	double time;
	long int bandwidth;
	
	bandwidth = MATRIX_bandwidth(A);
	if (msgs[0]) fprintf(stderr,"\n  [ REORDERING WEIGHTED SPECTRAL]\n");
	if (msgs[0]) fprintf(stderr,"  - Bandwidth before = %ld\n", bandwidth);
	
	/*---START TIME---------------> */ time =  get_TIME();	
	int i;
	int n = A->n;
	
	double *fvector = calloc (n ,sizeof(double));
	int    *list    = calloc (n ,sizeof(int));
	int    *info    = calloc (10,sizeof(int));
	int    *p       = calloc (n ,sizeof(int));
	int    *ia = A->IA;
	int    *ja = A->JA;
	double *a  = A->AA;
	int   lirn = A->nz;
	int    nnz = A->nz;
	
	/* -------------------------------------------------------------------- */    
	/* Convert matrix from 0-based C-notation to Fortran 1-based notation   */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < nnz; i++) {
		ja[i] += 1;
	}
	for (i = 0; i < n + 1; i++) {
		ia[i] += 1;
	}

	mc73_fiedler(&n,&lirn,ja,ia,list,fvector,info,a);

	/* -------------------------------------------------------------------- */    
	/* Convert matrix back to 0-based C-notation.                           */
	/* -------------------------------------------------------------------- */
	for (i = 0; i < nnz; i++) {
		ja[i] -= 1;
	}
	for (i = 0; i < n + 1; i++) {
		ia[i] -= 1;
	}
	
	ARRAY* R = malloc(n*sizeof(ARRAY));
	for (i = 0; i < n; ++i)
	{
		R[i].arr1 = fvector[i];
		R[i].arr2 = i;		
	}
	
	qsort (R,n,sizeof(ARRAY),COMPARE_eig); 
	
	for (i = 0; i < n; ++i) 
		p[i] = R[i].arr2; 

	free(R);
	free(info);
	free(list);
	free(fvector);
	
	MATRIX_ROW_permutation(A,p);
	MATRIX_COL_permutation(A,p);
	
	/*---PERMUTE COEF. VECTOR------ */
	double* bp = calloc(A->n,sizeof(double));		
	for (i = 0; i < A->n; ++i) 
		bp[i] = b[p[i]];
	
	memcpy(&b[0],&bp[0],A->n*sizeof(double));
	free(bp);
	
	/*---FINAL TIME---------------> */ time = (get_TIME() - time)/100.0;
	
	bandwidth = (int) MATRIX_bandwidth(A);
	if (msgs[0]) fprintf(stderr,"  - Bandwidth after  = %ld\n", bandwidth);
	if (msgs[0]) fprintf(stderr,"  - Elapsed time: %.6f sec\n", time);	

	(*P) = p;
	
	return;
}

int COMPARE_eig (const void * a, const void * b)
{
	if (((ARRAY*)a)->arr1 > ((ARRAY*)b)->arr1) return  1;
	if (((ARRAY*)a)->arr1 < ((ARRAY*)b)->arr1) return -1;
	return 0;
}


