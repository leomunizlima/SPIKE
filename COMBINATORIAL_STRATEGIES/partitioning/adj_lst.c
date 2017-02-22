# include <stdlib.h>
# include "adj_lst.h"

adj_lst * create_adj_lst () {
	adj_lst *l = (adj_lst *) malloc(sizeof(adj_lst));
	l->head = NULL;
	return l;
};

void destroy_adj_lst (adj_lst *l) {
	node_adj_lst *h;
	while (l->head) {
		h = l->head->next;
		free(l->head);
		l->head = h;
	}
	free(l);
};

void insert_adj_lst (adj_lst *l, int head, double cost) {
	node_adj_lst *h = (node_adj_lst *) malloc(sizeof(node_adj_lst));
	h->head = head;
	h->cost = cost;
	h->next = l->head;
	l->head = h;
};

node_adj_lst * remove_adj_lst (adj_lst *l, int head) {
	node_adj_lst *h = l->head, *aux;
	if (h->head==head) { l->head = h->next; return h; }
	for (aux=l->head; aux->next->head!=head; aux=aux->next) {}
	h = aux->next; aux->next = h->next;
	return h;
};

