typedef struct node_adj_lst {
	int head; double cost;
	struct node_adj_lst *next;
} node_adj_lst;

typedef struct {
	node_adj_lst *head;
} adj_lst;

adj_lst * create_adj_lst ();
void destroy_adj_lst (adj_lst *l);
void insert_adj_lst (adj_lst *l, int head, double cost);
node_adj_lst * remove_adj_lst (adj_lst *l, int head);

