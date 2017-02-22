/*----------------------------------------------------------------------------
 * QUADRATIC KNAPSACK PROBLEM REORDERING
 *--------------------------------------------------------------------------*/
#include "../COMMON_FILES/protos.h"


void calc_N (MAT **NR, MAT **NC, int m, int n, MAT *A, int s_i, int s_ip1, int *p, int *q)
{
	int i, ii, j, k, nz;
	Node *h;
	list **row, **col;
	*NR = (MAT *) malloc(sizeof(MAT)); *NC = (MAT *) malloc(sizeof(MAT));
	(*NR)->m = m; (*NC)->m = m;
	(*NR)->n = n; (*NC)->n = n;
	row = (list **) malloc(m*sizeof(list *)); for (i=0; i<m; i++) row[i] = create_l();
	col = (list **) malloc(n*sizeof(list *)); for (i=0; i<n; i++) col[i] = create_l();
	for (i=0; i<m; i++) {
		ii = p[s_i+i];
		for (k=A->IA[ii]; k<A->IA[ii+1]; k++) {
			if (q[A->JA[k]]>=s_ip1 && q[A->JA[k]]<s_ip1+n) {
				j = q[A->JA[k]] - s_ip1;
				insert_l(row[i], create_N(j, fabs(A->AA[k])));
				insert_l(col[j], create_N(i, fabs(A->AA[k])));
			}
		}
	}
	for (i=0; i<n; i++) {
		ii = p[s_ip1+i];
		for (k=A->IA[ii]; k<A->IA[ii+1]; k++) {
			if (q[A->JA[k]]>=s_i && q[A->JA[k]]<s_i+m) {
				j = q[A->JA[k]] - s_i;
				for (h=row[j]->head; h; h=h->next) { if (h->i==i) break; }
				if (h) {
					h->d += fabs(A->AA[k]);
					for (h=col[i]->head; h; h=h->next) {
						if (h->i==j) {
							h->d += fabs(A->AA[k]);
							break;
						}
					}
				} else {
					insert_l(row[j], create_N(i, fabs(A->AA[k])));
					insert_l(col[i], create_N(j, fabs(A->AA[k])));
				}
			}
		}
	}
	for (nz=0, i=0; i<m; i++) nz += row[i]->n;
	(*NR)->nz = nz; (*NC)->nz = nz;
	(*NR)->AA = (double *) malloc(nz*sizeof(double)); (*NC)->AA = (double *) malloc(nz*sizeof(double));
	(*NR)->D = NULL; (*NC)->D = NULL;
	(*NR)->IA = (int *) malloc((m+1)*sizeof(int)); (*NC)->IA = (int *) malloc(nz*sizeof(int));
	(*NR)->JA = (int *) malloc(nz*sizeof(int)); (*NC)->JA = (int *) malloc((n+1)*sizeof(int));
	(*NR)->IA[0] = 0;
	for (k=0, i=0; i<m; i++) {
		for (h=row[i]->head; h; h=h->next, k++) {
			(*NR)->JA[k] = h->i;
			(*NR)->AA[k] = h->d;
		}
		(*NR)->IA[i+1] = k;
	}
	(*NC)->JA[0] = 0;
	for (k=0, j=0; j<n; j++) {
		for (h=col[j]->head; h; h=h->next, k++) {
			(*NC)->IA[k] = h->i;
			(*NC)->AA[k] = h->d;
		}
		(*NC)->JA[j+1] = k;
	}
	for (i=0; i<m; i++) destroy_l(row[i]);
	for (i=0; i<n; i++) destroy_l(col[i]);
	free(row); free(col);
};

void swap_rowcol (int a, int b, int *p, int *q) 
{
	int aux = p[a];
	
	p[a] = p[b]; 
	p[b] = aux;
	
	q[p[a]] = a; 
	q[p[b]] = b;
}

/*----------------------------------------------------------------------------
 * QUADRATIC KNAPSACK PROBLEM reordering
 *--------------------------------------------------------------------------*/
void REORDERING_DE_MIN (MAT *A, double* b, int **K, int *s, int nP, int k, int *msgs) 
{
	double time;
	/*---START TIME---------------> */ time =  get_TIME(); 
	
	if (msgs[0])
	{
		fprintf(stderr,"\n  [ QUADRATIC KNAPSACK PROBLEM ]\n");	
		fprintf(stderr,"  - nnz outside of d+c blocks = %d (before): ", REORDERING_DE_MIN_nnzout(A,s,nP,k));
		fprintf(stderr,"\n");
	}
	
	int i, j, ind, ind2, nrow, ncol, *q, *p;
	double key;
	MAT *NR, *NC;
	node_fh *h, **rowN, **colN; fib_heap *rowH, *colH;
	rowH = make_heap(); 
	colH = make_heap();
	q = (int *) malloc(A->m*sizeof(int));
	p = (int *) malloc(A->m*sizeof(int));
	
	for (i = 0; i < A->m; i++) 
	{ 
		p[i] = i; 
		q[i] = i; 
	}
	for (ind = 0; ind < nP-1; ind++) 
	{
		nrow = s[ind+1]-s[ind]; 
		ncol = s[ind+2]-s[ind+1];
		
		calc_N(&NR, &NC, nrow, ncol, A, s[ind], s[ind+1], p, q);
		
		rowN = (node_fh **) malloc(nrow*sizeof(node_fh *)); 
		colN = (node_fh **) malloc(ncol*sizeof(node_fh *));
		
		for (i = 0; i < nrow; i++) 
		{
			for (key = 0.0, j = NR->IA[i]; j < NR->IA[i+1]; j++) 
				key += NR->AA[j];
			rowN[i] = insert_heap(rowH, i, key);
		}
		for (j = 0; j < ncol; j++) 
		{
			for (key = 0.0, i = NC->JA[j]; i < NC->JA[j+1]; i++) 
				key += NC->AA[i];
			colN[j] = insert_heap(colH, j, key);
		}
		i = nrow; 
		j = ncol;
		while (i>k || j>k)
		{
			if ((i > k && j > k && minimum_heap(rowH)->key <= minimum_heap(colH)->key) || (i > k && j==k))
			{
				h = extract_min_heap(rowH);
				for (ind2=NR->IA[h->v]; ind2<NR->IA[h->v+1]; ind2++)
				{
					if (colN[NR->JA[ind2]])
						decrease_key_heap(colH, colN[NR->JA[ind2]], colN[NR->JA[ind2]]->key-NR->AA[ind2]);
				}
				rowN[h->v] = NULL; 
				free(h);
				i--;
			} else 
			{
				h = extract_min_heap(colH);
				for (ind2=NC->JA[h->v]; ind2<NC->JA[h->v+1]; ind2++)
				{
					if (rowN[NC->IA[ind2]])
						decrease_key_heap(rowH, rowN[NC->IA[ind2]], rowN[NC->IA[ind2]]->key-NC->AA[ind2]);
				}
				colN[h->v] = NULL; 
				free(h);
				j--;
			}
		}
		i = nrow - k;
		while (rowH->n)
		{
			h = extract_min_heap(rowH);
			if (s[ind]+h->v<s[ind+1]-k)
			{
				while (rowN[i]) i++;
				swap_rowcol(s[ind]+h->v, s[ind]+i, p, q);
				i++;
			}
			free(h);
		}
		j = 0;
		while (colH->n)
		{
			h = extract_min_heap(colH);
			if (s[ind+1]+h->v>=s[ind+1]+k)
			{
				while (colN[j]) j++;
				swap_rowcol(s[ind+1]+j, s[ind+1]+h->v, p, q);
				j++;
			}
			free(h);
		}
		free(rowN); 
		free(colN);
		MATRIX_clean(NR); 
		MATRIX_clean(NC);
	}
	free(rowH); 
	free(colH);
	free(q);
	
	MATRIX_ROW_permutation(A,p);
	MATRIX_COL_permutation(A,p);
	
	/*---PERMUTE COEF. VECTOR------ */
	double* bp = calloc(A->n,sizeof(double));		
	for (i = 0; i < A->n; ++i) 
		bp[i] = b[p[i]];
	
	memcpy(&b[0],&bp[0],A->n*sizeof(double));
	free(bp);
	
	/*---FINAL TIME---------------> */ time = (get_TIME() - time)/100.0;
	
	if (msgs[0])
	{
		fprintf(stderr,"  - nnz outside of d+c blocks = %d (after): ", REORDERING_DE_MIN_nnzout(A,s,nP,k));
		fprintf(stderr,"\n");
		fprintf(stderr,"  - Elapsed time: %.6f sec\n", time);
	}
	
	(*K) = p;
	
}

int REORDERING_DE_MIN_nnzout (MAT *A, int *s, int nP, int k)
{
	int i, ii, j, c=0;
	for (i = s[0]; i < s[1]; i++)
	{
		for (j = A->IA[i]; j < A->IA[i+1]; j++)
			if ((i < s[1]-k && A->JA[j] >= s[1]) || (i >= s[1]-k && A->JA[j] >= s[1]+k)) c++;
	}
	for (i=1; i<nP-1; i++)
	{
		for (ii=s[i]; ii<s[i+1]; ii++)
		{
			for (j = A->IA[ii]; j < A->IA[ii+1]; j++)
				if ((ii < s[i]+k && A->JA[j] < s[i]-k) || (ii >= s[i]+k && A->JA[j] < s[i]) ||
				    (ii < s[i+1]-k && A->JA[j] >= s[i+1]) || (ii >= s[i+1]-k && A->JA[j] >= s[i+1]+k)) c++;
		}
	}
	for (i=s[nP-1]; i<s[nP]; i++)
	{
		for (j = A->IA[i]; j < A->IA[i+1]; j++)
			if ((i < s[nP-1]+k && A->JA[j] < s[nP-1]-k) || (i >= s[nP-1]+k && A->JA[j] < s[nP-1])) c++;
	}
	return c;
}
