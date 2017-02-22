# include <float.h>

typedef struct node_fh {
	int v, degree, mark; double key;
	struct node_fh *p, *child, *left, *right;
} node_fh;

typedef struct {
	int n;
	node_fh *min;
} fib_heap;

fib_heap * make_heap ();
node_fh * insert_heap (fib_heap *H, int v, double k);
node_fh * minimum_heap (fib_heap *H);
node_fh * extract_min_heap (fib_heap *H);
void decrease_key_heap (fib_heap *H, node_fh *x, double k);
void delete_heap (fib_heap *H, node_fh *x);

