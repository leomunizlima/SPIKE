# include <stdlib.h>
# include "list.h"

list * create_l () {
	list *l = (list *) malloc(sizeof(list));
	l->n = 0;
	l->head = NULL;
	return l;
};

Node * create_N (int i, double d) {
	Node *h = (Node *) malloc(sizeof(Node));
	h->i = i; h->d = d;
	return h;
};

void destroy_l (list *l) {
	Node *h;
	while (l->head) {
		h = l->head->next;
		free(l->head);
		l->head = h;
	}
	free(l);
};

void insert_l (list *l, Node *h) {
	h->next = l->head;
	l->head = h;
	l->n++;
};

