/*
 *  Main header file for percolation code.
 */

/*
 *  System size L
 */

#define L 288

/*
 *  Although overall system is square, i.e. size L x L, we will define
 *  different variables for the first and second dimensions. This is
 *  because, in the parallel code, the local arrays will not be
 *  square. For example, using a simple 1D decomposition over P
 *  processes, then M = L/P and N = L
 */

#define M L
#define N L


#define ndims 1
#define TRUE  1
#define FALSE 0


/*
 *  Prototypes for supplied functions
 */

/*
 *  Visualisation
 */

void percwrite(char *percfile, int map[M][N], int ncluster);

/*
 *  Random numbers
 */

void rinit(int ijkl);
float uni(void);



int **alloc_2d_int(int rows,int cols);
void display_matrix(int** matrix,int rows,int cols);