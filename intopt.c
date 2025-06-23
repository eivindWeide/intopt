#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct list_t list_t;

struct simplex_t {
  int m, n, *var;
  double **a, *b, *c, *x, y;
};

typedef struct {
  int m, n, k, h;
  double xh, ak, bk, *min, *max, **a, *b, *c, *x, z;
} node_t;

struct list_t {
        list_t* pred;
        list_t* succ;
        node_t*  data;
};

void pivot(struct simplex_t *s, int row, int col);
double xsimplex(int m, int n, double **a, double *b, double *c, double *x, double y, int *var, int h);

void* xmalloc(size_t size)
{
        void* ptr = malloc(size);

        if (ptr == NULL) {
                //fprintf(stderr, "out of memory\n");
                exit(1);
        }

        return ptr;
}

void free_node(node_t **node)
{
        int i;
        node_t *p;
        
        p = *node;

        if (p == NULL)
                return;

        free(p->b);
        free(p->c);
        free(p->x);
        free(p->min);
        free(p->max);
        for(i = 0; i < p->m+1; ++i) 
                free(p->a[i]);
        free(p->a);
        free(p);
        
        *node = NULL;
}

list_t *new_list(void *data)
{
        list_t *list;

        list = xmalloc(sizeof(list_t));
        list->succ = list;
        list->pred = list;
        list->data = data;
        return list;
}

void free_list(list_t **list)
{
        list_t *p;
        list_t *q;

        p = *list;

        if (p == NULL)
                return;

        p->pred->succ = NULL;
        do {
                q = p->succ;
                free_node(&(p->data));
                free(p->data);
                p->data = NULL;
                free(p);
                p = q;
        } while (p != NULL);

        *list = NULL;
}

void insert_last(list_t **list, node_t* data)
{
        
        //printf("INSERT LAST\n");
        list_t* p = *list;
        list_t* q;

        if (p == NULL) {
                *list = new_list(data);
                return;
        }

        q = p;
        do {
                if ((node_t*)q->data == data){
                        //printf("Already in list\n");
                        return;
                        }
                q = q->succ;
        } while (q != p);

        q = xmalloc(sizeof(list_t));

        q->data = data;
        q->pred = p->pred;
        q->succ = p;
        p->pred->succ = q;
        p->pred = q;
}

node_t *pop(list_t **list)
{
        list_t *current, *next;
        node_t *d;
        
        current = *list;
        next = current->succ;
        d = current->data;

        if (current == current->succ) {
                free(current);
                *list = NULL;
        } else {
                current->succ->pred = current->pred;
                current->pred->succ = current->succ;
                free(current);
                *list = next;
        }      
        
        return d;
}

void delete_z_less_than(list_t **list, double z)
{
        list_t* first = *list;
        list_t* current = *list;
        list_t* next;

        int i, skip;


        if (current == NULL)
                return;

        

        
        do {    
                if (current->succ == current && current->data->z < z) {
                        free_node(&(current->data));
                        free(current);
                        *list = NULL;
                        return;
                }
\
                skip = 0;
                next = current->succ;
                if (current->data->z < z) {
                        current->pred->succ = current->succ;
                        current->succ->pred = current->pred;
                        if (current == first) {
                                first = current->succ;
                                skip = 1;
                        }
                        //printf("%p\n", current->data);
                        free_node(&(current->data));
                        free(current);
                }
                current = next;
                
                
        } while (first != next || skip);
        *list = current;
        return;
}

int init(struct simplex_t *s, int m, int n, double **a, 
        double *b, double *c, double *x, double y, int *var)
{
        //printf("init!\n");
        int i, k;

        s->m = m;
        s->n = n;
        s->a = a;
        s->b = b;
        s->c = c;
        s->x = x;
        s->y = y;
        s->var = var;

        

        if (s->var == NULL) {
                s->var = xmalloc((m + n + 1) * sizeof(int));
                for (i = 0; i < m + n; ++i)
                        s->var[i] = i;
        }
        for (k = 0, i = 1; i < m; ++i) 
                if (b[i] < b[k])
                        k = i;
        
        return k;
}

int select_nonbasic(struct simplex_t *s)
{
        int i;
        for (i = 0; i < s->n; ++i)
                if (s->c[i] > 1e-6) 
                        return i;
                
        return -1;
}

void prepare(struct simplex_t *s, int k)
{
        //printf("prepare!\n");
        int m = s->m;
        int n = s->n;
        int i;
        double *t;

        for (i = m + n; i > n; --i)
                s->var[i] = s->var[i-1];
        s->var[n] = m + n;
        n = n + 1;
        for (i = 0; i < m; ++i)
                s->a[i][n-1] = -1; 

        s->x = xmalloc((m + n) * sizeof(double));
        t = xmalloc((n - 1) * sizeof(double));

        for (i = 0; i < n + m; ++i) 
                s->x[i] = 0;

        for (i = 0; i < n - 1; ++i)
                t[i] = s->c[i];

        
        s->c = xmalloc(n * sizeof(double));

        for (i = 0; i < n - 1; ++i)
                s->c[i] = 0;

        free(t);
        s->c[n - 1] = -1;
        s->n = n;
        pivot(s, k, n-1);

}

int initial(struct simplex_t *s, int m, int n, double **a,
        double *b, double *c, double *x, double y, int *var)
{
        //printf("initial!\n");
        int i, j, k;
        double w;
        k = init(s, m, n, a, b, c, x, y, var);
        if (b[k] >= 0) 
                return 1; //feasible

        prepare(s, k);
        n = s->n;
        s->y = xsimplex(m, n, s->a, s->b, s->c, s->x, 0, s->var, 1);
        for (i = 0; i < m + n; ++i) {
                if (s->var[i] == m + n - 1) {
                        
                        if (fabs(s->x[i]) > 1e-6) {
                                free(s->x);
                                free(s->c);
                                return 0;
                        } else {
                                break;
                        }
                }
        }
        if (i >= n) {
                // x_n+m is basic. find good nonbasic.
                for (j = k = 0; k < n; k = k + 1) 
                        if (fabs(s->a[i-n][k]) > fabs(s->a[i - n][j])) 
                                j = k;
                        
                pivot(s, i - n, j);
                i = j;
        }
        if (i < n - 1) {
                k = s->var[i]; 
                s->var[i] = s->var[n-1]; 
                s->var[n-1] = k;
                for (k = 0; k < m; ++k) {
                        w = s->a[k][n-1];
                        s->a[k][n-1] = s->a[k][i];
                        s->a[k][i] = w;
                }
        }
        
        free(s->c);
        
        s->c = c;
        s->y = y;

        for (k = n-1; k < n+m-1; ++k)
                s->var[k] = s->var[k+1];
        s->n = s->n - 1;
        n = s->n;
        double *t;
        t = xmalloc(n * sizeof(double));

        for (i = 0; i < n; ++i) 
                t[i] = 0;

        for (k = 0; k < n; ++k) {
                for (j = 0; j < n; ++j) {
                        if (k == s->var[j]) {
                                // x_k is nonbasic. add c_k
                                t[j] = t[j] + s->c[k];
                                goto next_k;
                        }
                }
                for (j = 0; j < m; ++j)
                        if (k == s->var[n+j])
                                break;
                s->y = s->y + s->c[k] * s->b[j];
                for (i = 0; i < n; ++i) 
                        t[i] = t[i] - s->c[k] * s->a[j][i];
                next_k:;
        }
        
        for (i = 0; i < n; ++i) {
                
                s->c[i] = t[i];
        }
        free(t);
        free(s->x);
        return 1;
}

void pivot(struct simplex_t *s, int row, int col)
{
        //printf("pivot! to rows: %d, cols: %d\n", row, col);

        double **a = s->a;
        double *b = s->b;
        double *c = s->c;
        int m = s->m;
        int n = s->n;
        int i, j, t;
        

        t = s->var[col];
        s->var[col] = s->var[n+row];
        s->var[n+row] = t;
        s->y = s->y + c[col] * b[row] / a[row][col];

        for (i = 0; i < n; ++i) 
                if (i != col) 
                        c[i] = s->c[i] - s->c[col] * a[row][i] / a[row][col];
                
        
        c[col] = - c[col] / a[row][col];
        for (i = 0; i < m; ++i) 
                if (i != row) 
                        b[i] = b[i] - a[i][col] * b[row] / a[row][col];
                
        for (i = 0; i < m; ++i) 
                if (i != row) 
                        for (j = 0; j < n; ++j) 
                                if (j != col) 
                                        a[i][j] = a[i][j] - a[i][col] * a[row][j] / a[row][col];
        
        for (i = 0; i < m; ++i) 
                if (i != row)
                        a[i][col] = -a[i][col] / a[row][col];
        
        for (i = 0; i < n; i = ++i) 
                if (i != col)
                        a[row][i] = a[row][i] / a[row][col];

        b[row] = b[row] / a[row][col];
        a[row][col] = 1 / a[row][col];
}


double xsimplex(int m, int n, double **a, double *b, double *c, 
                double *x, double y, int *var, int h) 
{
        //printf("xsimplex!, m: %d, n: %d\n", m, n);
        struct simplex_t s;
        

        int i, row, col;
        if (!initial(&s, m, n, a, b, c, x, y, var)) {
                free(s.var);
                return NAN; 
        }
        
        while ((col = select_nonbasic(&s)) >= 0) {
                row = -1;
                for (i = 0; i < m; ++i) {
                        if ((a[i][col] > 1e-6) && 
                           (row < 0 || b[i] / a[i][col] < 
                           b[row] / a[row][col])) {
                                row = i;
                           }
                }
                if (row < 0) {
                        free(s.var);
                        return INFINITY; // unbounded
                }
                pivot(&s, row, col);
        }
        if (h == 0) {
                for (i = 0; i < n; ++i)
                        if (s.var[i] < n) 
                                x[s.var[i]] = 0;
                        
                for (i = 0; i < m; ++i) 
                        if (s.var[n+i] < n) 
                                x[s.var[n+i]] = s.b[i];
                        
                free(s.var);
                
        } else {
                for (i = 0; i < n; ++i) 
                        x[i] = 0;
                for (i = n; i < n + m; ++i) 
                        x[i] = s.b[i-n];
        }
        return s.y;
}

double simplex (int m, int n, double **a, double *b, double *c, double *x, double y)
{   
        return xsimplex(m, n, a, b, c, x, y, NULL, 0);
}

node_t *initial_node(int m, int n, double **a, double *b, double *c)
{
        //printf("initial_node!\n");
        int i, j;
        node_t *p;
        
        p = xmalloc(sizeof(node_t));
        p->a = xmalloc((m+1) * sizeof(double));
        for (i = 0; i < m+1; ++i)
                p->a[i] = xmalloc((n+1) * sizeof(double));
        p->b = xmalloc((m+1) * sizeof(double));
        p->c = xmalloc((n+1) * sizeof(double));
        p->x = xmalloc((n+1) * sizeof(double));
        p->min = xmalloc(n * sizeof(double));
        p->max = xmalloc(n * sizeof(double));
        p->m = m;
        p->n = n;

        for (i = 0; i < m; ++i){
                for (j = 0; j < n; ++j){
                        p->a[i][j] = a[i][j];
                }
                
        }
        
        for (i = 0; i < m; ++i)
                p->b[i] = b[i];

        for (i = 0; i < n; ++i){
                p->c[i] = c[i];
                p->min[i] = -INFINITY;
                p->max[i] = INFINITY;
        }
        return p;
}

node_t *extend(node_t *p, int m, int n, double **a, double *b, double *c, int k, double ak, double bk)
{
        //printf("extend!\n");
        node_t *q;
        int i, j;

        q = xmalloc(sizeof(node_t));

        q->k = k;
        q->ak = ak;
        q->bk = bk;

        if (ak > 1e-6 && p->max[k] < INFINITY) {
                q->m = p->m;
        } else if (ak < -1e-6 && p->min[k] > 1e-6) {
                q->m = p->m;
        } else {
                q->m = p->m + 1;
        }
        q->n = p->n;
        q->h = -1;
        q->a = xmalloc((q->m+1) * sizeof(double));
        for (i = 0; i < q->m+1; ++i)
                q->a[i] = xmalloc((q->n+1) * sizeof(double));
        q->b = xmalloc((q->m+1) * sizeof(double));
        q->c = xmalloc((q->n+1) * sizeof(double));
        q->x = xmalloc((q->n+1) * sizeof(double));
        q->min = xmalloc(n * sizeof(double));
        q->max = xmalloc(n * sizeof(double));


        for (i = 0; i < n; ++i) {
                q->c[i] = c[i];
                q->min[i] = p->min[i];
                q->max[i] = p->max[i];
        }

        for (i = 0; i < m; ++i){
                for (j = 0; j < n; ++j)
                        q->a[i][j] = a[i][j];
                
                for (j = n; j < q->n+1; ++j){
                        q->a[i][j] = 0;
                        q->c[j] = 0;
                }

                q->b[i] = b[i];
        }
        for (i = m; i < q->m+1; ++i) {
                for (j = 0; j < q->n+1; ++j) 
                        q->a[i][j] = 0;
        
                q->b[i] = 0;
        }
        for (i = n; i < q->n+1; ++i){
                q->x[i] = 0;        
                q->c[i] = 0;
        }
                
        if (ak > 1e-6) {
                if (q->max[k] == INFINITY || bk < q->max[k]){
                        q->max[k] = bk;
                }
        } else if (q->min[k] == -INFINITY || -bk > q->min[k]) {
                q->min[k] = -bk;
        }

        for (i = m, j = 0; j < n; ++j) {
                if (q->min[j] > -INFINITY) {
                        q->a[i][j] = -1;
                        q->b[i] = -q->min[j];
                        i = i + 1;
                }
                if (q->max[j] < INFINITY) {
                        q->a[i][j] = 1;
                        q->b[i] = q->max[j];
                        i = i + 1;
                }
        }
        return q;
}

int is_integer(double *xp)
{
        double x = *xp;
        double r = round(x);
        if (fabs(r - x) < 1e-6) {
                *xp = r;
                return 1;
        } else {
                return 0;
        }
}

int integer(node_t *p)
{
        int i;

        for (i = 0; i < p->n; ++i)
                if (!is_integer(&p->x[i]))
                        return 0;
        return 1;
}

void bound(node_t *p, list_t **h, double *zp, double *x)
{
        int i;
        //printf("bound\n");
        //printf("best z: %lf, %lf\n", p->z, *zp);
        if (p->z > *zp) {
                *zp = p->z;

                for (i = 0; i < p->n; ++i)
                        x[i] = p->x[i];
                delete_z_less_than(h, p->z); 
        }
}

int branch(node_t *q, double z)
{
        //printf("branch\n");
        double min, max;
        int h;

        if (q->z < z)
                return 0;

        for (h = 0; h < q->n; ++h) {
                if (!is_integer(&q->x[h])) {
                        if (q->min[h] == -INFINITY) {
                                min = 0;
                        } else {
                                min = q->min[h];
                        }
                        max = q->max[h];
                        
                        if (floor(q->x[h]) < min || ceil(q->x[h]) > max)
                                continue;
                                

                        //printf("litentest\n");
                        q->h = h;
                        q->xh = q->x[h];
                        //free_node(&q);
                        return 1;
                }
                        
        }
        return 0;
}

void succ(node_t *p, list_t **h, int m, int n, double **a, double *b, 
        double *c, int k, double ak, double bk, double *zp, double *x) 
{
        //printf("succ with: %lf\n", bk);

        int i;
        node_t *q = extend(p, m, n, a, b, c, k, ak, bk);

        if (q == NULL){
                //printf("q is NULL");
                return;
        }
        
        //printf("%d\n", q->m);

        q->z = simplex(q->m, q->n, q->a, q->b, q->c, q->x, 0);

        //printf("%lf\n", q->z);

        if (isfinite(q->z)) {
                if (integer(q))
                        bound(q, h, zp, x);
                else if (branch(q, *zp)) {
                        insert_last(h, q);
                        return;
                }
        }

        free_node(&q);
        return;
}

double intopt(int m, int n, double **a, double *b, double *c, double *x)
{
        int i;
        node_t *p = initial_node(m, n, a, b, c);
        list_t *h = new_list(p);
        double z = -INFINITY;
        p->z = simplex(p->m, p->n, p->a, p->b, p->c, p->x, 0);

        if (integer(p) || !isfinite(p->z)) {
                z = p->z;
                if (integer(p)) {
                       for (i = 0; i < n; ++i) 
                                x[i] = p->x[i];
                }
                free_list(&h);
                return z;
        }
        branch(p, z);
        while (h != NULL) {
                p = pop(&h);
                succ(p, &h, m, n, a, b, c, p->h, 1, floor(p->xh), &z, x);
                succ(p, &h, m, n, a, b, c, p->h, -1, -ceil(p->xh), &z, x);
                free_node(&p);
        }

        free_list(&h);

        if (z == -INFINITY)
                return NAN;
        else
                return z;
}


int main(void)
{       
        size_t i, j;
        size_t rows, cols;
        double *xvector, *cvector, *bvector, **matrix, y;
        
        scanf("%zu %zu", &rows, &cols);
        //allocate to heap because we don't know size at compile time
        xvector = xmalloc((cols + 1) * sizeof(double));
        cvector = xmalloc(cols * sizeof(double));
        bvector = xmalloc(rows * sizeof(double));
        matrix = xmalloc(rows * sizeof(double));

        for (i = 0; i < rows; ++i)
                matrix[i] = xmalloc((cols + 1) * sizeof(double));

        for (i = 0; i < cols; ++i) 
                scanf("%lf", &cvector[i]);

        for (i = 0; i < rows; ++i)
                for (j = 0; j < cols; ++j)
                        scanf("%lf", &matrix[i][j]);

        for (i = 0; i < rows; ++i)
                scanf("%lf", &bvector[i]);
        
        y = intopt(rows, cols, matrix, bvector, cvector, xvector);

        printf("result is: %lf\n", y);

        int loop;

        printf("Coefficients are: ");
        for(loop = 0; loop < cols; loop++)
                printf("%f*x%d ", xvector[loop], loop);
        printf("\n");

        free(xvector);
        free(cvector);
        free(bvector);
        
        for(i = 0; i < rows; ++i) 
                free(matrix[i]);
        free(matrix);

        return 0;
}