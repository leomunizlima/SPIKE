# include <stdlib.h>
# include "AL_graph.h"

AL_graph * create_AL (int n) {
	int i;
	AL_graph *G = (AL_graph *) malloc(sizeof(AL_graph));
	G->n = n;
	G->adj_lst = (adj_lst **) malloc(n*sizeof(adj_lst *));
	for (i=0; i<n; i++) G->adj_lst[i] = create_adj_lst();
	return G;
};

void destroy_AL (AL_graph *G) {
	int i;
	for (i=0; i<G->n; i++) destroy_adj_lst(G->adj_lst[i]);
	free(G->adj_lst);
	free(G);
};

void insert_arc (AL_graph *G, int tail, int head, double cost) {
	insert_adj_lst(G->adj_lst[tail-1], head, cost);
};

node_adj_lst * remove_arc (AL_graph *G, int tail, int head) {
	return remove_adj_lst(G->adj_lst[tail-1], head);
};

