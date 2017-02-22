#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "./COMMON_FILES/protos.h"

int main (int argc, char* argv[])
{
	int RANK,NPROC;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &RANK);
	MPI_Comm_size(MPI_COMM_WORLD, &NPROC);
	
	if (argc != 13 && RANK == 0)
	{
		fprintf(stderr,"\n Error!"); 
		fprintf(stderr,"\n Usage: ./program <matrix.mtx> <GMRES_tol> <GMRES_k> <GMRES_maxit> <scaling> <matching> <reordering> <knapsack> <k> <result.txt>. Exiting... [MAIN]\n\n");
		exit(0);
	}
	
	/*-------------------------------------------------*/
	/*---SETTING VARIABLES-----------------------------*/
	/*-------------------------------------------------*/
	int    i, release, bandwidth, n;

	double *dw;	// Scaling factors
	int    *Q; 	// Matching reordering 
	int    *P; 	// Reordering to minimize bandwidth
	int    *K; 	// Quadratic Knapsack reordering 
	
	int    nP      = NPROC;			// number os processors
	char*  entry   = argv[1];		// file entry name
	double tol     = atof(argv[2]);		// GMRES tolerance
	int    restart = atoi(argv[3]);		// GMRES restart
	int    maxit   = atoi(argv[4]);		// GMRES max iterations
	
	int*    msgs   = malloc(5*sizeof(int));
	        msgs[0] = 1;			// printf messages statistics
	        msgs[1] = atoi(argv[5]);	// 0 = no scaling
	        msgs[2] = atoi(argv[6]);	// 0 = no matching
	        msgs[3] = atoi(argv[7]);	// 0 = no reordering, 1 = spectral, 2 = weighted spectral
	        msgs[4] = atoi(argv[8]);	// 0 = no knapsack problem
	
	int    k       = atoi(argv[9]);	// size of coupling blocks
	char*  output  = argv[10];		// file output name
	
	/*-------------------------------------------------------------------*/
	/*---SERIAL PART (PREPROCESSING)-------------------------------------*/
	/*-------------------------------------------------------------------*/
	if (RANK == 0)
	{	
		MAT    *A = malloc(sizeof(MAT));
		double *b;
		 
		/*-------------------------------------------------*/
		/*---READING ENTRY FILES---------------------------*/
		/*-------------------------------------------------*/
		READ_INPUT (A,&b,entry,msgs);
		n = A->n; // n means the global size of the matrix
        	        
		fprintf(stderr,"\n  [ PARAMETERS ]\n");
		fprintf(stderr,"  - matrix: %s | n: %d | nnz: %d | mpi: %d | openmp: %d | <s,m,r,k> : <%d,%d,%d,%d> | k: %d |", \
			entry,A->n,A->nz,nP,atoi(getenv("OMP_NUM_THREADS")),msgs[1],msgs[2],msgs[3],msgs[4],k);

		/*-------------------------------------------------*/
		/*---MATCHING AND SCALING--------------------------*/
		/*-------------------------------------------------*/
		if (msgs[1] || msgs[2])
			REORDERING_MATCHING (A,b,&Q,&dw,msgs);
                
		/*-------------------------------------------------*/
		/*---REORDERING------------------------------------*/
		/*-------------------------------------------------*/
		if (msgs[3] == 1) REORDERING_SPECTRAL     (A,b,&P,msgs);
		if (msgs[3] == 2) REORDERING_SPECTRAL_WGT (A,b,&P,msgs);
		bandwidth = MATRIX_bandwidth(A);
	
		/*-------------------------------------------------*/
		/*---PARTITIONING----------------------------------*/
		/*-------------------------------------------------*/
		int *s,*L;
		PARTITIONING_MIN_MAX (A,nP,&s,&L,msgs);


		/*-------------------------------------------------*/
		/*---KNAPSACK PROBLEM------------------------------*/
		/*-------------------------------------------------*/
		if (msgs[4]) REORDERING_DE_MIN (A,b,&K,s,nP,k,msgs);

		/*-------------------------------------------------*/
		/*---WRITE AND FINALIZE----------------------------*/
		/*-------------------------------------------------*/
		MATRIX_writeCSR(A,b,s,nP,k,msgs);
		free(s);
		free(L);
				
		
		/*-------------------------------------------------*/
		/*---WAITING..---------- --------------------------*/
		/*-------------------------------------------------*/
		for (i = 1; i < NPROC; ++i)
			MPI_Send (&release,1,MPI_INT,i,0,MPI_COMM_WORLD);
	}
	else 
	{ 
		MPI_Recv (&release,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);		
	}
	
	/*-------------------------------------------------------------------*/
	/*---PARALLEL PART (PROCESSING)--------------------------------------*/
	/*-------------------------------------------------------------------*/
	/*-------------------------------------------------*/
	/*---ALLOCATING STRUCTURES-------------------------*/
	/*-------------------------------------------------*/
	MAT *A,*A_bef,*C,*A_aft,*B;
	double *b;
	A = (MAT*) malloc(sizeof(MAT));
	if (RANK != 0)
	{
		A_bef = (MAT*) malloc(sizeof(MAT));
		C     = (MAT*) malloc(sizeof(MAT));
	}
	if (RANK != NPROC - 1)
	{
		A_aft = (MAT*) malloc(sizeof(MAT));
		B     = (MAT*) malloc(sizeof(MAT));
	}
	
	MPI_Bcast               (&bandwidth,1,MPI_INT,0,MPI_COMM_WORLD);
	PARALLEL_MATRIX_readCSR (A_bef,A,A_aft,C,B,&b,bandwidth,RANK,NPROC,msgs);

	double *x1 = calloc(A->n,sizeof(double));
	double *x2 = calloc(A->n,sizeof(double));


	/*-------------------------------------------------*/
	/*---PARALLEL GMRES--------------------------------*/
	/*-------------------------------------------------*/
 	if (atoi(argv[11])){
		 PARALLEL_GMRES (A_bef,A,A_aft,x1,b,tol,restart,maxit,RANK,NPROC,msgs);
			
		 double *x = PARALLEL_assemble_solution (x1,A->n,n,RANK,NPROC);
		
 		 WRITE_SOLUTION (x,entry,output,n,P,K,dw,RANK,msgs);
		 free(x);	
	}

	/*-------------------------------------------------*/
	/*---PARALLEL PRECONDITIONED GMRES-----------------*/
	/*-------------------------------------------------*/
	PARALLEL_PGMRES (A_bef,A,A_aft,C,B,x2,b,tol,restart,maxit,RANK,NPROC,msgs);

	/*-------------------------------------------------*/
	/*---ASSEMBLE AND WRITE SOLUTION-------------------*/
	/*-------------------------------------------------*/
	MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);	
	double *x = PARALLEL_assemble_solution (x2,A->n,n,RANK,NPROC);
	
	WRITE_SOLUTION (x,entry,output,n,P,K,dw,RANK,msgs);
	
	/*-------------------------------------------------*/
	/*---DEALLOCATE STRUTURES AND FINALIZE-------------*/
	/*-------------------------------------------------*/
	if (RANK == 0){
		free(x);
		if (msgs[1]) free(dw);
		if (msgs[2]) free(Q); 
		if (msgs[3]) free(P);
	}

	MATRIX_clean(A);
	if (RANK != 0)
	{
		MATRIX_clean(A_bef);
		MATRIX_clean(C);
	}
	if (RANK != NPROC - 1)
	{
		MATRIX_clean(A_aft);
		MATRIX_clean(B);
	}
	free(b);
	free(x1);
	free(x2);
	free(msgs);
	
	MPI_Finalize();
	return 0;
}


