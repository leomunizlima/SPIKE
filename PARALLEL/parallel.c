/*----------------------------------------------------------------------------
 * PARALLEL FUNCTIONS
 *--------------------------------------------------------------------------*/
#include "../COMMON_FILES/protos.h"

/*----------------------------------------------------------------------------
 * Forms the dot product of two vectors in parallel
 *--------------------------------------------------------------------------*/
inline double PARALLEL_ddot (int n, double *x, double *y)
{
	long int i, m;
	double stemp,sfinal;

	stemp = 0.0;
	m = n-4;

	for (i = 0; i < m; i += 5)
		stemp += x[i] * y[i] + x[i+1] * y[i+1] + x[i+2] * y[i+2] + x[i+3] * y[i+3] + x[i+4] * y[i+4];

	for ( ; i < n; i++)        /* clean-up loop */
		stemp += x[i] * y[i];

	MPI_Allreduce(&stemp, &sfinal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	return sfinal;
}

/*----------------------------------------------------------------------------
 * Parallel GMRES algorithm
 *--------------------------------------------------------------------------*/
double PARALLEL_GMRES_algorithm (MAT* A_bef, MAT* A, MAT* A_aft, double *x, double *b, int n, int k, int *iter, double eps, int RANK, int NPROC)
{ 
	int i, j, hz = k + 1;  

	double * ax = calloc(n,sizeof(double));
	double * r  = calloc(n,sizeof(double));
	double * p1 = calloc(n,sizeof(double));
	double * p2 = calloc(n,sizeof(double));
  
	double e0;
	double beta;
	double ret;

	PARALLEL_e_MATRIX_matvec (A_bef,A,A_aft,x,ax,RANK,NPROC);
	
	ddiff(n,b,ax,r);
  
	e0 = sqrt (PARALLEL_ddot(n,r,r));

	if (e0 <= eps)
	{
		free(ax);
		free(r);
		free(p1);
		free(p2);
		return e0;
	} 
 
	double** q = malloc(hz*sizeof(double*));
	for (i = 0; i < hz; ++i) 
		q[i] = malloc(n*sizeof(double));
  
	double** h = calloc(hz,sizeof(double*));
	for (i = 0; i < hz; ++i) 
		h[i] = calloc(hz,sizeof(double));
  
	double* e = malloc(hz*sizeof(double)); 
	double* s = malloc(hz*sizeof(double));
	double* c = malloc(hz*sizeof(double));
  
	dscal (n,1.0/e0,r,q[0]); 
  
	e[0]  = e0;
	*iter = 0;
  
	for (j = 0; e0 > eps && j < k; ++j)
	{
		*iter= *iter +1;

		PARALLEL_e_MATRIX_matvec (A_bef,A,A_aft,q[j],ax,RANK,NPROC);
		
		double nr1 = sqrt (PARALLEL_ddot(n,ax,ax));

		for (i = 0; i <= j; ++i)
		{
			h[i][j] = PARALLEL_ddot (n,q[i],ax);
			daxpy (n,-h[i][j],q[i],ax);
		}
    
		h[j+1][j] = sqrt(PARALLEL_ddot (n,ax,ax));
		double nr2 = 0.001 * h[j+1][j] + nr1;
    
		if (fabs(nr2  - nr1) < eps)
		{
			for (i = 0; i <= j; ++i)
			{
				double hr = PARALLEL_ddot (n,q[i],ax);
				h[i][j] += hr;
				daxpy (n,-hr,q[i],ax);
			}
			h[j+1][j] = sqrt(PARALLEL_ddot (n,ax,ax));
		}

		for (i = 0; i <= j - 1; ++i)
		{
			double x  = h[i][j];
			double y  = h[i+1][j];
			h[i][j]   = x * c[i+1] + y * s[i+1];
			h[i+1][j] = x * s[i+1] - y * c[i+1];
		}

		beta    = sqrt(h[j][j]*h[j][j] + h[j+1][j]*h[j+1][j]);
		s[j+1]  = h[j+1][j] / beta;
		c[j+1]  = h[j][j]   / beta;
		h[j][j] = beta;
		e[j+1]  = s[j+1] * e[j];
		e[j]    = c[j+1] * e[j];
		e0      = e[j+1];
   
		dscal(n,1.0 / h[j+1][j],ax,q[j+1]); 
	}
	--j;
  
	ret = e[j+1];
	double * y = malloc(hz * sizeof(double));
	for (i = j; i >= 0; --i) 
	{
		double sum = 0.0;
		for (k = i + 1; k <= j; ++k)
			sum += h[i][k] * y[k];
		y[i] = (e[i] - sum) / h[i][i];
	}
	for (i = 0; i <= j; ++i) 
		daxpy (n,y[i],q[i],x);
	free(y);
  
	/* Free up memory */
	for (i=0;i<hz;++i) 
		free(q[i]); free(q); 
	for (i=0;i<hz;++i) 
		free(h[i]); free(h);
	free(ax);
	free(r); 
	free(e);
	free(s); 
	free(c);
	free(p1);
	free(p2);
	return ret;
}

/*----------------------------------------------------------------------------
 * Parallel GMRES algorithm call
 *--------------------------------------------------------------------------*/
void PARALLEL_GMRES (MAT* A_bef, MAT* A, MAT* A_aft, double *x, double *b, double tol, int restart, int maxiter, int RANK, int NPROC, int* msgs)
{  
	double time;
	/*---START TIME---------------> */ time =  get_TIME(); 
	
	int     i, it;  
	int     n = A->n;

	double bn = sqrt(PARALLEL_ddot(n,b,b));
	double e;
	
	if (RANK == 0 && msgs[0]) 
	{
		fprintf(stderr,"\n  [ PARALLEL GMRES ]");
		fprintf(stderr, "\n  - ||b|| : %le\n", bn);
	}

	if (bn < tol) 
		return;

	for (i = 0; i < maxiter; i++)
	{
		e  = PARALLEL_GMRES_algorithm (A_bef, A, A_aft, x, b, n, restart, &it, tol * bn, RANK, NPROC);
		e /= bn;
		if (RANK == 0 && msgs[0]) 
			fprintf(stderr, "  - ITERS : %d - %d, eps = %le\n", i+1,i*restart + it, e);
		if (e < tol)
			break;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	/*---FINAL TIME---------------> */ time = (get_TIME() - time)/100.0;
	if (msgs[0]) fprintf (stderr,"    --RANK (%d): Elapsed Time = %lf sec\n",RANK,time);
	
}

/*----------------------------------------------------------------------------
 * Parallel Preconditioned GMRES algorithm
 *--------------------------------------------------------------------------*/
double PARALLEL_PGMRES_algorithm (MAT* A_bef, MAT* A, MAT* A_aft, MAT* C, MAT* B, MAT *S, void* A_pt[64], int A_iparm[64], double A_dparm[64], void* S_pt[64], int S_iparm[64], double S_dparm[64], double *x, double *b, int n, int k, int *iter, double eps, int RANK, int NPROC)
{ 
	int i, j, hz = k + 1;  

	double * ax = calloc(n,sizeof(double));
	double * r  = calloc(n,sizeof(double));
	double * p1 = calloc(n,sizeof(double));
	double * p2 = calloc(n,sizeof(double));
  
	double e0;
	double beta;
	double ret;

	PARALLEL_e_MATRIX_matvec (A_bef,A,A_aft,x,ax,RANK,NPROC);
	ddiff(n,b,ax,r);
	PARALLEL_SPIKE_algorithm (C,A,B,S,A_pt,A_iparm,A_dparm,S_pt,S_iparm,S_dparm,r,r,RANK,NPROC);
	
	e0 = sqrt (PARALLEL_ddot(n,r,r));

	if (e0 <= eps)
	{
		free(ax);
		free(r);
		free(p1);
		free(p2);
		return e0;
	} 
 
	double** q = malloc(hz*sizeof(double*));
	for (i = 0; i < hz; ++i) 
		q[i] = malloc(n*sizeof(double));
  
	double** h = calloc(hz,sizeof(double*));
	for (i = 0; i < hz; ++i) 
		h[i] = calloc(hz,sizeof(double));
  
	double* e = malloc(hz*sizeof(double)); 
	double* s = malloc(hz*sizeof(double));
	double* c = malloc(hz*sizeof(double));
  
	dscal (n,1.0/e0,r,q[0]); 
  
	e[0]  = e0;
	*iter = 0;
  
	for (j = 0; e0 > eps && j < k; ++j)
	{
		*iter= *iter +1;

		PARALLEL_e_MATRIX_matvec (A_bef,A,A_aft,q[j],ax,RANK,NPROC);
		PARALLEL_SPIKE_algorithm (C,A,B,S,A_pt,A_iparm,A_dparm,S_pt,S_iparm,S_dparm,ax,ax,RANK,NPROC);
		
		double nr1 = sqrt (PARALLEL_ddot(n,ax,ax));

		for (i = 0; i <= j; ++i)
		{
			h[i][j] = PARALLEL_ddot (n,q[i],ax);
			daxpy (n,-h[i][j],q[i],ax);
		}
    
		h[j+1][j] = sqrt(PARALLEL_ddot (n,ax,ax));
		double nr2 = 0.001 * h[j+1][j] + nr1;
    
		if (fabs(nr2  - nr1) < eps)
		{
			for (i = 0; i <= j; ++i)
			{
				double hr = PARALLEL_ddot (n,q[i],ax);
				h[i][j] += hr;
				daxpy (n,-hr,q[i],ax);
			}
			h[j+1][j] = sqrt(PARALLEL_ddot (n,ax,ax));
		}

		for (i = 0; i <= j - 1; ++i)
		{
			double x  = h[i][j];
			double y  = h[i+1][j];
			h[i][j]   = x * c[i+1] + y * s[i+1];
			h[i+1][j] = x * s[i+1] - y * c[i+1];
		}

		beta    = sqrt(h[j][j]*h[j][j] + h[j+1][j]*h[j+1][j]);
		s[j+1]  = h[j+1][j] / beta;
		c[j+1]  = h[j][j]   / beta;
		h[j][j] = beta;
		e[j+1]  = s[j+1] * e[j];
		e[j]    = c[j+1] * e[j];
		e0      = e[j+1];
   
		dscal(n,1.0 / h[j+1][j],ax,q[j+1]); 
	}
	--j;
  
	ret = e[j+1];
	double * y = malloc(hz * sizeof(double));
	for (i = j; i >= 0; --i) 
	{
		double sum = 0.0;
		for (k = i + 1; k <= j; ++k)
			sum += h[i][k] * y[k];
		y[i] = (e[i] - sum) / h[i][i];
	}
	for (i = 0; i <= j; ++i) 
		daxpy (n,y[i],q[i],x);
	free(y);
  
	/* Free up memory */
	for (i=0;i<hz;++i) 
		free(q[i]); free(q); 
	for (i=0;i<hz;++i) 
		free(h[i]); free(h);
	free(ax);
	free(r); 
	free(e);
	free(s); 
	free(c);
	free(p1);
	free(p2);
	return ret;
}

/*----------------------------------------------------------------------------
 * Parallel Preconditioned GMRES algorithm call
 *--------------------------------------------------------------------------*/
void PARALLEL_PGMRES (MAT* A_bef, MAT* A, MAT* A_aft, MAT* C, MAT* B, double *x, double *b, double tol, int restart, int maxiter, int RANK, int NPROC, int* msgs)
{  
	double time, time_max, time_min,time_tot;
	int j, COUNT=1;

	time_min=1e10;
	time_max = 0;
	time_tot = 0;

	for (j=0; j<COUNT;j++)
	{
	/*---START TIME---------------> */ time =  get_TIME(); 

		
	int     i, it;  
	int     n = A->n;
	
	memset(x,0,n*sizeof(double));
	double bn = sqrt(PARALLEL_ddot(n,b,b));
	double e;
	
	MPI_Request request;
	
	if (bn < tol) 
		return;

	MAT* S = malloc (sizeof(MAT));
	
	double *fB, *fC;
	double *V, *W, *WT;
	
	/* Pardiso parameters. */
	void*       A_pt[64];
	int         A_iparm[64];
	double      A_dparm[64];
	void*       S_pt[64];
	int         S_iparm[64];
	double      S_dparm[64];
		
	if (RANK == 0)
	{
		int size = B->m*B->n;
		int    k = B->m;
		fB = calloc (size,sizeof(double));
		V  = calloc (k*k ,sizeof(double));
		WT = calloc (k*k ,sizeof(double));

		if(!fB) printf("OPS FALHOU EM fB! \n");

		MATRIX_CSRto_FULL   (B,fB);
		PARDISO_bottom_tips (A,fB,V,k);

		MPI_Recv(WT,k*k,MPI_DOUBLE,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

		PARALLEL_REDUCED_system (V,WT,k,S,RANK);

		PARDISO_numerical_factorization (A,A_pt,A_iparm,A_dparm,1);
		PARDISO_numerical_factorization (S,S_pt,S_iparm,S_dparm,1);
	}
	else if (RANK == NPROC - 1)
	{
		int size = C->m*C->n;
		int    k = C->m;		
		fC = calloc (size,sizeof(double));
		W  = calloc (k*k ,sizeof(double));

		if(!fB) printf("OPS FALHOU EM fB! \n");

		MATRIX_CSRto_FULL (C,fC);
		PARDISO_top_tips  (A,fC,W,k);

		MPI_Isend(W,k*k,MPI_DOUBLE,RANK - 1,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);
		
		PARDISO_numerical_factorization (A,A_pt,A_iparm,A_dparm,1);
	}
	else 
	{	  
		int size = C->m*C->n;
		int    k = C->m;
		fC = calloc (size,sizeof(double));
		fB = calloc (size,sizeof(double));
		V  = calloc (k*k ,sizeof(double));
		W  = calloc (k*k ,sizeof(double));
		WT = calloc (k*k ,sizeof(double));

		MATRIX_CSRto_FULL   (C,fC);
		MATRIX_CSRto_FULL   (B,fB);
		PARDISO_bottom_tips (A,fB,V,k);
		PARDISO_top_tips    (A,fC,W,k);

		MPI_Recv (WT,k*k,MPI_DOUBLE,RANK + 1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Isend(W ,k*k,MPI_DOUBLE,RANK - 1,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);

		PARALLEL_REDUCED_system (V,WT,k,S,RANK);
		
		PARDISO_numerical_factorization (A,A_pt,A_iparm,A_dparm,1);
		PARDISO_numerical_factorization (S,S_pt,S_iparm,S_dparm,1);
	}
	
	if (RANK == 0 && msgs[0])
	{
		fprintf(stderr,"\n  [ PARALLEL GMRES WITH SPIKE PRECONDITIONER ]");
		fprintf(stderr, "\n  - ||b|| : %le\n", bn);
	}
	
// 	MPI_Barrier(MPI_COMM_WORLD);

	for (i = 0; i < maxiter; i++)
	{
		e  = PARALLEL_PGMRES_algorithm (A_bef,A,A_aft,C,B,S,A_pt,A_iparm,A_dparm,S_pt,S_iparm,S_dparm,x,b,n,restart,&it,tol * bn,RANK,NPROC);
		MPI_Barrier(MPI_COMM_WORLD);
		e /= bn;
		if (RANK == 0 && msgs[0]) 
			fprintf(stderr, "  - ITERS : %d - %d, eps = %le\n", i+1,i*restart + it, e);
		if (e < tol)
		{
			if (RANK == 0)
			{
				PARDISO_release_memory (A,A_pt,A_iparm,A_dparm,1);
				PARDISO_release_memory (S,S_pt,S_iparm,S_dparm,1);
				free(fB);
				free(V); 
				free(WT);
				MATRIX_clean(S);
				
			}
			else if (RANK == NPROC - 1)
			{
				PARDISO_release_memory (A,A_pt,A_iparm,A_dparm,1);
				free(fC);
				free(W);
				free(S);
			}
			else 
			{
				PARDISO_release_memory (A,A_pt,A_iparm,A_dparm,1);
				PARDISO_release_memory (S,S_pt,S_iparm,S_dparm,1);
				free(fB);
				free(fC);
				free(V);
				free(W);
				free(WT);
				MATRIX_clean(S);
				
			}
			
			break;
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	/*---FINAL TIME---------------> */ time = (get_TIME() - time)/100.0;
	if (time<time_min)
		time_min = time;
	if (time>time_max)
		time_max = time;
	time_tot += time;
	}

	time = time_tot/COUNT;

	if (msgs[0]) fprintf (stderr,"    --RANK (%d): Elapsed Time = %lf sec (Time min=%lf sec) (Time max=%lf sec)\n",RANK,time, time_min, time_max);
}

/*----------------------------------------------------------------------------
 * Storing the matrices in CSR format
 *--------------------------------------------------------------------------*/
void PARALLEL_MATRIX_readCSR (MAT* A_bef, MAT* A, MAT* A_aft, MAT* C, MAT* B, double** b, int bandwidth, int RANK, int NPROC, int *msgs)
{
	char file[50];
	int i;
	char header[1024];
	int sc, m, r, knp; 

	sc = msgs[1]; 
	m = msgs[2]; 
	r = msgs[3]; 
	knp = msgs[4]; 



	if (RANK == 0)
	{
		sprintf(file,"data0%d_%d_%d_%d_%d.mtx",RANK, sc, m, r, knp);

		FILE* f = fopen(file,"r");		
		if (!f){ fprintf (stderr,"Can't open %s\n. Exiting.. [PARALLEL_MATRIX_readCSR]",file); exit(0); }

		MATRIX_readCSR (A, f);     fgets(header,200,f); fgets(header,200,f);
		MATRIX_readCSR (A_aft, f); fgets(header,200,f); fgets(header,200,f);
		MATRIX_readCSR (B, f);     fgets(header,200,f); fgets(header,200,f);
		
		int n = A->n;
		   *b = (double*) malloc (n*sizeof(double));
		for (i = 0; i < n; ++i)
			fscanf(f,"%lf\n",&(*b)[i]);
		
		int max = 0;
		for (i = 0; i < A_aft->n; ++i)
                        if (A_aft->IA[i+1] != 0)
			if (max < A_aft->JA[A_aft->IA[i+1] - 1])
				max = A_aft->JA[A_aft->IA[i+1] - 1];
		
		A_aft->bandwidth[0] = max + 1;
		MPI_Send(&(A_aft->bandwidth[0]),1,MPI_INT,1,0,MPI_COMM_WORLD);			 // Sending A_aft_size to RANK 1
		MPI_Recv(&(A_aft->bandwidth[1]),1,MPI_INT,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE); // Receiving A_bef_size from RANK 1
		
		fclose(f);
	}
	else if (RANK == NPROC - 1)
	{
		sprintf(file,"data0%d_%d_%d_%d_%d.mtx",NPROC - 1, sc, m, r, knp);

		FILE* f = fopen(file,"r");		
		if (!f){ fprintf (stderr,"Can't open %s\n. Exiting.. [PARALLEL_MATRIX_readCSR]",file); exit(0); }

		MATRIX_readCSR (A_bef, f); fgets(header,200,f); fgets(header,200,f);
		MATRIX_readCSR (C, f);     fgets(header,200,f); fgets(header,200,f);
		MATRIX_readCSR (A, f);     fgets(header,200,f); fgets(header,200,f);
		
		int n = A->n;
		   *b = (double*) malloc (n*sizeof(double));
		for (i = 0; i < n; ++i)
			fscanf(f,"%lf\n",&(*b)[i]);
		
		int min = n;
		for (i = 0; i < A_bef->n; ++i)
                        if (A_bef->IA[i+1] != A_bef->nz)
			if (min > A_bef->JA[A_bef->IA[i]] && A_bef->IA[i] != A_bef->IA[i+1])
				min = A_bef->JA[A_bef->IA[i]];

		A_bef->bandwidth[0] = A_bef->m - min;
		MPI_Send(&(A_bef->bandwidth[0]),1,MPI_INT,RANK - 1,0,MPI_COMM_WORLD);			// Sending A_bef_size to RANK p - 1
		MPI_Recv(&(A_bef->bandwidth[1]),1,MPI_INT,RANK - 1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE); // Receiving A_aft_size from RANK p - 1
		
		fclose(f);
	}
	else
	{
		sprintf(file,"data0%d_%d_%d_%d_%d.mtx",RANK, sc, m, r, knp);


		FILE* f = fopen(file,"r");
		if (!f){ fprintf (stderr,"Can't open %s\n. Exiting.. [PARALLEL_MATRIX_readCSR]",file); exit(0); }

		MATRIX_readCSR (A_bef, f); fgets(header,200,f); fgets(header,200,f);
		MATRIX_readCSR (C, f);     fgets(header,200,f); fgets(header,200,f);
		MATRIX_readCSR (A, f);     fgets(header,200,f); fgets(header,200,f);
		MATRIX_readCSR (A_aft, f); fgets(header,200,f); fgets(header,200,f);
		MATRIX_readCSR (B, f);     fgets(header,200,f); fgets(header,200,f);

		int n = A->n;
		   *b = (double*) malloc (n*sizeof(double));
		for (i = 0; i < n; ++i)
			fscanf(f,"%lf\n",&(*b)[i]);
		
		int max = 0;
		for (i = 0; i < A_aft->n; ++i)
                        if (A_aft->IA[i+1] != 0)
			if (max < A_aft->JA[A_aft->IA[i+1] - 1])
				max = A_aft->JA[A_aft->IA[i+1] - 1];
		int min = n;
		for (i = 0; i < A_bef->n; ++i)
                        if (A_bef->IA[i+1] != A_bef->nz)
			if (min > A_bef->JA[A_bef->IA[i]] && A_bef->IA[i] != A_bef->IA[i+1])
				min = A_bef->JA[A_bef->IA[i]];
		
		A_aft->bandwidth[0] = max + 1;
		MPI_Send(&(A_aft->bandwidth[0]),1,MPI_INT,RANK + 1,0,MPI_COMM_WORLD);			 // Sending A_aft_size to RANK p + 1
		MPI_Recv(&(A_aft->bandwidth[1]),1,MPI_INT,RANK + 1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);  // Receiving A_bef_size from RANK p + 1
		
		A_bef->bandwidth[0] = A_bef->m - min;
		MPI_Send(&(A_bef->bandwidth[0]),1,MPI_INT,RANK - 1,0,MPI_COMM_WORLD);			// Sending A_bef_size to RANK p - 1
		MPI_Recv(&(A_bef->bandwidth[1]),1,MPI_INT,RANK - 1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE); // Receiving A_aft_size from RANK p - 1
		
		fclose(f);
	}
}

/*----------------------------------------------------------------------------
 * Parallel matrix vector product
 *--------------------------------------------------------------------------*/
void PARALLEL_MATRIX_matvec (MAT *A_bef, MAT* A, MAT* A_aft, double* x, double* r, int RANK, int NPROC)
{
	int i;
	MPI_Request request;
	
	if (RANK == 0)
	{	  
		int           n = A->n;
		int      mA_aft = A_aft->m;
		double* x_aft   = calloc (mA_aft,sizeof(double));
		double* b_temp1 = calloc (n ,sizeof(double));
		double* b_temp2 = calloc (n ,sizeof(double));

		MPI_Isend(x    ,n     ,MPI_DOUBLE,1,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);
		MPI_Recv (x_aft,mA_aft,MPI_DOUBLE,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		
		
		MATRIX_matvec (A    ,x    ,b_temp1);	/* A(i)*x(i)   */	
		MATRIX_matvec (A_aft,x_aft,b_temp2);	/* A_aft(i)*x(i+1) */
		
		for (i = 0; i < n; ++i)			/* assembling results*/
			r[i] = b_temp1[i] + b_temp2[i]; 
		
// 		for (i = 0; i < n; ++i)	
// 			fprintf (stderr,"A: (%d) %lf\n",i, r[i]);
		
		free(x_aft);
		free(b_temp1);
		free(b_temp2);
	}
	else if (RANK == NPROC - 1)
	{
		int           n = A->n;
		int      mA_bef = A_bef->m;
		double* x_bef   = calloc (mA_bef,sizeof(double));
		double* b_temp1 = calloc (n ,sizeof(double));
		double* b_temp2 = calloc (n ,sizeof(double));
		
		MPI_Isend(x    ,n     ,MPI_DOUBLE,RANK - 1,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);
		MPI_Recv (x_bef,mA_bef,MPI_DOUBLE,RANK - 1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

		MATRIX_matvec (A_bef,x_bef,b_temp1);	/* C(i)*x(i-1) */
		MATRIX_matvec (A    ,x    ,b_temp2);	/* A(i)*x(i)   */
		
		for (i = 0; i < n; ++i)			/* assembling results*/
			r[i] = b_temp1[i] + b_temp2[i];

/*		for (i = 0; i < n; ++i)	
			fprintf (stderr,"A: (%d) %lf\n",i, r[i]);*/	
		
		free(x_bef);
		free(b_temp1);
		free(b_temp2);	
	}
	else
	{
		int           n = A->n;
		int      mA_aft = A_aft->m;
		int      mA_bef = A_bef->m;
		double* x_aft   = calloc (mA_aft,sizeof(double));
		double* x_bef   = calloc (mA_bef,sizeof(double));
		double* b_temp1 = calloc (n ,sizeof(double));
		double* b_temp2 = calloc (n ,sizeof(double));
		double* b_temp3 = calloc (n ,sizeof(double));
		         
		MPI_Isend(x    ,n     ,MPI_DOUBLE,RANK - 1,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);
		MPI_Recv (x_bef,mA_bef,MPI_DOUBLE,RANK - 1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Isend(x    ,n     ,MPI_DOUBLE,RANK + 1,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);
		MPI_Recv (x_aft,mA_aft,MPI_DOUBLE,RANK + 1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
		
		MATRIX_matvec (A_bef,x_bef,b_temp1);	/* C(i)*x(i-1) */
		MATRIX_matvec (A    ,x    ,b_temp2);	/* A(i)*x(i)   */	
		MATRIX_matvec (A_aft,x_aft,b_temp3);	/* B(i)*x(i+1) */
		
		for (i = 0; i < n; ++i)			/* assembling results*/
			r[i] = b_temp1[i] + b_temp2[i] + b_temp3[i];
		
// 		for (i = 0; i < n; ++i)	
// 			00fprintf (stderr,"A: (%d) %lf\n",i, r[i]);
		
		free(x_aft);
		free(x_bef);
		free(b_temp1);
		free(b_temp2);
		free(b_temp3);
	}
}

/*----------------------------------------------------------------------------
 * Parallel matrix vector product
 *--------------------------------------------------------------------------*/
void PARALLEL_e_MATRIX_matvec (MAT *A_bef, MAT* A, MAT* A_aft, double* x, double* r, int RANK, int NPROC)
{
	int i;
	MPI_Request request;
	
	if (RANK == 0)
	{	  
		int           n = A->n;
		double* x_aft   = calloc (A_aft->m,sizeof(double));
		double* b_temp1 = calloc (n ,sizeof(double));
		double* b_temp2 = calloc (n ,sizeof(double));

		MPI_Isend(&x[n-A_aft->bandwidth[1]],A_aft->bandwidth[1],MPI_DOUBLE,1,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);
		MPI_Recv (&x_aft[0]                ,A_aft->bandwidth[0],MPI_DOUBLE,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
		MATRIX_matvec (A    ,x    ,b_temp1);	/* A(i)*x(i)   */	
		MATRIX_matvec (A_aft,x_aft,b_temp2);	/* A_aft(i)*x(i+1) */
		
		for (i = 0; i < n; ++i)			/* assembling results*/
			r[i] = b_temp1[i] + b_temp2[i]; 
		
// 		for (i = 0; i < n; ++i)	
// 			fprintf (stderr,"NOVO: (%d) %lf\n",RANK, r[i]);
		
		free(x_aft);
		free(b_temp1);
		free(b_temp2);
	}
	else if (RANK == NPROC - 1)
	{
		int           n = A->n;
		double* x_bef   = calloc (A_bef->m,sizeof(double));
		double* b_temp1 = calloc (n ,sizeof(double));
		double* b_temp2 = calloc (n ,sizeof(double));
		
		MPI_Isend(&x[0]                               ,A_bef->bandwidth[1],MPI_DOUBLE,RANK - 1,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);
		MPI_Recv (&x_bef[A_bef->m-A_bef->bandwidth[0]],A_bef->bandwidth[0],MPI_DOUBLE,RANK - 1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
		MATRIX_matvec (A_bef,x_bef,b_temp1);	/* C(i)*x(i-1) */
		MATRIX_matvec (A    ,x    ,b_temp2);	/* A(i)*x(i)   */
		
		for (i = 0; i < n; ++i)			/* assembling results*/
			r[i] = b_temp1[i] + b_temp2[i];
/*		
		for (i = 0; i < n; ++i)	
			fprintf (stderr,"NOVO: (%d) %lf\n",RANK, r[i]);*/

		free(x_bef);
		free(b_temp1);
		free(b_temp2);	
	}
	else
	{
		int           n = A->n;
		double* x_aft   = calloc (A_aft->m,sizeof(double));
		double* x_bef   = calloc (A_bef->m,sizeof(double));
		double* b_temp1 = calloc (n ,sizeof(double));
		double* b_temp2 = calloc (n ,sizeof(double));
		double* b_temp3 = calloc (n ,sizeof(double));
		         
		MPI_Isend(&x[n-A_aft->bandwidth[1]]           ,A_aft->bandwidth[1],MPI_DOUBLE,RANK + 1,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);
		MPI_Isend(&x[0]                               ,A_bef->bandwidth[1],MPI_DOUBLE,RANK - 1,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);
		MPI_Recv (&x_aft[0]                           ,A_aft->bandwidth[0],MPI_DOUBLE,RANK + 1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		MPI_Recv (&x_bef[A_bef->m-A_bef->bandwidth[0]],A_bef->bandwidth[0],MPI_DOUBLE,RANK - 1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
		MATRIX_matvec (A_bef,x_bef,b_temp1);	/* C(i)*x(i-1) */
		MATRIX_matvec (A    ,x    ,b_temp2);	/* A(i)*x(i)   */	
		MATRIX_matvec (A_aft,x_aft,b_temp3);	/* B(i)*x(i+1) */
		
		for (i = 0; i < n; ++i)			/* assembling results*/
			r[i] = b_temp1[i] + b_temp2[i] + b_temp3[i];
		
// 		for (i = 0; i < n; ++i)	
// 			fprintf (stderr,"NOVO: (%d) %lf\n",RANK, r[i]);
// 		
		free(x_aft);
		free(x_bef);
		free(b_temp1);
		free(b_temp2);
		free(b_temp3);
	}
}

/*----------------------------------------------------------------------------
 * Construct reduced system S for SPIKE preconditioner
 *--------------------------------------------------------------------------*/
void PARALLEL_REDUCED_system (double* V, double* W, int k, MAT* S, int RANK)
{
	int i,j,aa,ja;
	
	int nz = (2*k)+(2*k*k);
	int n  =  2*k;
	
	S->AA  = calloc(nz ,sizeof(double));
	S->D   = calloc(n  ,sizeof(double));
	S->JA  = calloc(nz ,sizeof(int   ));
	S->IA  = calloc(n+1,sizeof(int   ));
	S->n   = n;
	S->m   = n;
	S->nz  = nz;
	
	/* IA */
	for (i = 1; i < 2*k + 1; ++i)
		S->IA[i] = S->IA[i-1] + k + 1;

	/* AA, JA, g */
	aa = 0;
	ja = 0;
	for (i = 0; i < k; ++i)
	{
		S->AA[aa++] = 1.0;
		S->JA[ja++] = i;
		for (j = 0; j < k; ++j)
		{
			S->AA[aa++] = V[i*k + j];
			S->JA[ja++] = j + k;
		}
	}
	
	for (i = k; i < 2*k; ++i)
	{
		for (j = 0; j < k; ++j)
		{
			S->AA[aa++] = W[(i-k)*k + j];
			S->JA[ja++] = j;
		}
		S->AA[aa++] = 1.0;
		S->JA[ja++] = i;
	}
}

/*----------------------------------------------------------------------------
 * Parallel matrix vector product
 *--------------------------------------------------------------------------*/
void PARALLEL_SPIKE_algorithm (MAT *C, MAT* A, MAT* B, MAT *S, void* A_pt[64], int A_iparm[64], double A_dparm[64], void* S_pt[64], int S_iparm[64], double S_dparm[64], double* f, double* y, int RANK, int NPROC)
{
	int i,j;
	MPI_Request request;
	
	double *g, *gT, *gS;
	double *x, *xB, *xT, *xTp;
	double *r1,*r2,*p;
	
	if (RANK == 0)
	{
		int n = A->n;
		int k = B->m;
		g  = calloc (n  ,sizeof(double));
		gT = calloc (k  ,sizeof(double));
		gS = calloc (2*k,sizeof(double));
		x  = calloc (2*k,sizeof(double));
		xB = calloc (k  ,sizeof(double));
		xT = calloc (k  ,sizeof(double));
		r1 = calloc (n  ,sizeof(double));
		p  = calloc (n  ,sizeof(double));
		
		PARDISO_back_substitution (A,A_pt,A_iparm,A_dparm,f,g,1);
		
		MPI_Recv(gT,k,MPI_DOUBLE,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
		memcpy(&gS[0],&g[n-k],k*sizeof(double));
		memcpy(&gS[k],&gT[0] ,k*sizeof(double));
		
		PARDISO_back_substitution (S,S_pt,S_iparm,S_dparm,gS,x,1);
		
		memcpy(&xB[0],&x[0],k*sizeof(double));
		memcpy(&xT[0],&x[k],k*sizeof(double));
		
		MPI_Isend(xB,k,MPI_DOUBLE,1,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);
		
		MATRIX_matvec (B,xT,r1);
		for (i = 0; i < n; ++i)
			p[i] = f[i] - r1[i];
		
		PARDISO_back_substitution (A,A_pt,A_iparm,A_dparm,p,y,1);
	}
	else if (RANK == NPROC - 1)
	{
		int      n   = A->n;
		int      k   = C->m;
		g  = calloc (n,sizeof(double));
		xB = calloc (k,sizeof(double));
		r1 = calloc (n,sizeof(double));
		p  = calloc (n,sizeof(double));
		
		PARDISO_back_substitution (A,A_pt,A_iparm,A_dparm,f,g,1);
		
		MPI_Isend(g, k,MPI_DOUBLE,RANK - 1,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);
		MPI_Recv (xB,k,MPI_DOUBLE,RANK - 1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
		MATRIX_matvec (C,xB,r1);
		for (i = 0; i < n; ++i)
			p[i] = f[i] - r1[i];
		
		PARDISO_back_substitution (A,A_pt,A_iparm,A_dparm,p,y,1);
	}
	else
	{
		int     n = A->n;
		int     k = B->m;
		g  = calloc (n  ,sizeof(double));
		gT = calloc (k  ,sizeof(double));
		gS = calloc (2*k,sizeof(double));
		x  = calloc (2*k,sizeof(double));
		xB = calloc (k  ,sizeof(double));
		xT = calloc (k  ,sizeof(double));
		xTp = calloc(k  ,sizeof(double));
		r1 = calloc (n  ,sizeof(double));
		r2 = calloc (n  ,sizeof(double));
		p  = calloc (n  ,sizeof(double));

		PARDISO_back_substitution (A,A_pt,A_iparm,A_dparm,f,g,1);
		
		MPI_Isend(g, k,MPI_DOUBLE,RANK - 1,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);
		MPI_Recv (gT,k,MPI_DOUBLE,RANK + 1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
		memcpy(&gS[0],&g[n-k],k*sizeof(double));
		memcpy(&gS[k],&gT[0] ,k*sizeof(double));
		
		PARDISO_back_substitution (S,S_pt,S_iparm,S_dparm,gS,x,1);
		
		memcpy(&xB[0],&x[0],k*sizeof(double));
		memcpy(&xT[0],&x[k],k*sizeof(double));		
	
		MPI_Isend(xB ,k,MPI_DOUBLE,RANK + 1,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);
		MPI_Recv (xTp,k,MPI_DOUBLE,RANK - 1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

		MATRIX_matvec (C,xTp,r1);
		MATRIX_matvec (B,xT,r2);
		for (i = 0; i < n; ++i)
			p[i] = f[i] - r1[i] - r2[i];
		
		PARDISO_back_substitution (A,A_pt,A_iparm,A_dparm,p,y,1);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if (RANK == 0)
	{
		free(g); 
		free(gT);
		free(gS);
		free(x);
		free(xB);
		free(xT);
		free(r1);
		free(p);
	}
	else if (RANK == NPROC - 1)
	{
		free(g);
		free(xB);
		free(r1);
		free(p);
	}
	else
	{
		free(g);  
		free(gT); 
		free(gS); 
		free(x);  
		free(xB); 
		free(xT); 
		free(xTp);
		free(r1); 
		free(r2); 
		free(p); 
	}
}

double* PARALLEL_assemble_solution (double* x, int n, int N, int RANK, int NPROC)
{
	MPI_Request request;
	int     *array_n;
	double **array_x;
	double  *assembled_x; 
	int i,j;
	
	array_n = calloc(NPROC, sizeof(int));
	MPI_Gather(&n,1,MPI_INT,&array_n[RANK],1,MPI_INT,0,MPI_COMM_WORLD);
	
		
	if (RANK == 0)
	{
		assembled_x = calloc(N,sizeof(double));
		array_x     = calloc(NPROC,sizeof(double*)); 
		for (i = 0; i < NPROC; ++i)
			array_x[i] = calloc(array_n[i],sizeof(double));
		
		for (i = 0; i < array_n[0]; ++i)
			array_x[0][i] = x[i];
		
		for (i = 1; i < NPROC; ++i)
			MPI_Recv(array_x[i],array_n[i],MPI_DOUBLE,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		
		int k = 0;
		for (i = 0; i < NPROC; ++i)
			for (j = 0; j < array_n[i]; ++j)
				assembled_x[k++] = array_x[i][j];
			
		for (i = 0; i < NPROC; ++i)
			free(array_x[i]);
	
		free(array_x);

	}else{
		MPI_Isend(x,n,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&request); MPI_Request_free(&request);		
	}
	MPI_Barrier(MPI_COMM_WORLD);
	free(array_n);	
	if (RANK == 0) 
		return assembled_x;

	return 0;		
}

