# include "adj_lst.h"

typedef struct { // adjacency list graph storage
	int n;
	adj_lst **adj_lst;
} AL_graph;

AL_graph * create_AL (int n);
void destroy_AL (AL_graph *G);
void insert_arc (AL_graph *G, int tail, int head, double cost);
node_adj_lst * remove_arc (AL_graph *G, int tail, int head);

