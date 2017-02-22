# ifndef list_h
# define list_h

typedef struct Node {
	int i; double d;
	struct Node *next;
} Node;

typedef struct {
	int n;
	Node *head;
} list;

list * create_l ();
Node * create_N (int i, double d);
void destroy_l (list *l);
void insert_l (list *l, Node *h);

# endif

