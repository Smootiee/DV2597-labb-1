/***************************************************************************
 *
 * Sequential version of Gaussian elimination
 *
 ***************************************************************************/

#include <stdio.h>
#include <pthread.h>

#define MAX_SIZE 4096
#define NR_OF_CPUS 8
typedef double matrix[MAX_SIZE][MAX_SIZE];

int	N;		/* matrix size		*/
int	maxnum;		/* max number of element*/
char	*Init;		/* matrix init type	*/
int	PRINT;		/* print switch		*/
matrix	A;		/* matrix A		*/
double	b[MAX_SIZE];	/* vector b             */
double	y[MAX_SIZE];	/* vector y             */
int nr_of_threads;
int v;

pthread_barrier_t sync_barrier, main_thread_sync_barrier;
pthread_mutex_t lock;


struct thread_args {
    int k;
    int thread_id;
};


/* forward declarations */
void work_seq(void);
void Init_Matrix(void);
void Print_Matrix(void);
void Init_Default(void);
int Read_Options(int, char **);
void work_par();
void work_par2();

int
main(int argc, char **argv)
{
    int i, timestart, timeend, iter;

    Init_Default();		/* Init default values	*/
    Read_Options(argc,argv);	/* Read arguments	*/
    Init_Matrix();		/* Init the matrix	*/
    if (v == 1) {
        work_seq();
    }
    else if (v == 3) {
        work_par2();
    }
    if (PRINT == 1)
	   Print_Matrix();
}

void work_seq(){
    int i, j, k;
    for (k = 0; k < N; ++k){

        for (j = k + 1; j < N; ++j) {
            A[k][j] = A[k][j] / A[k][k]; /* Division step */
        }


        y[k] = b[k] / A[k][k];
        A[k][k] = 1.0;
        for (i = k+1; i < N; ++i) {
            for (j = k+1; j < N; ++j)
                A[i][j] = A[i][j] - A[i][k]*A[k][j]; /* Elimination step */
            b[i] = b[i] - A[i][k]*y[k];
            A[i][k] = 0.0;
        }
    }
}

void *division_step(struct thread_args *args){
    int j, k = args->k;

    for (j = k+1+args->thread_id; j < N; j+=nr_of_threads) {
        A[k][j] = A[k][j] / A[k][k]; /* Division step */
    }

    pthread_barrier_wait(&main_thread_sync_barrier);
}

void *elimination_step(struct thread_args *args){
    int i, j;
    int k = args->k;

    /*
        The entire elimination step in the parallelized function where the work load is divided in a round-robin principle.
        The thread_id is the shift for the start point and each iteration for i needs to jump `nr_of_threads` steps.
    */
    for (i = k+1+args->thread_id; i < N; i+=nr_of_threads) {
        for (j = k+1; j < N; ++j) {
            A[i][j] = A[i][j] - A[i][k]*A[k][j];  /* Elimination step */
        }
        b[i] = b[i] - A[i][k]*y[k];
        A[i][k] = 0.0;
        
    }
    pthread_barrier_wait(&main_thread_sync_barrier);

}

void work_par2(){
    int i, j, k;
    pthread_t *threads = (pthread_t *)malloc(sizeof(pthread_t)*nr_of_threads);
    struct thread_args *args = (struct thread_args*)malloc(sizeof(struct thread_args)*nr_of_threads);

    /* 
        Initializing barrier `main_thread_sync_barrier` such that it will need to be called by nr_of_threads + 1
        threads. This is in effect a way to stop the main thread from going further or exiting.
    */
    pthread_barrier_init(&main_thread_sync_barrier, NULL, nr_of_threads + 1);

    for (k = 0; k < N; ++k){
        /*
            This is a change which does not matter. Moved the calculation of y_k to the beginning of the for-loop
            so that the for-loop can do the division step on all elements in the row vector A_k.
        */
        y[k] = b[k] / A[k][k];
        for (j = N - 1; j >= k; --j) {
            A[k][j] = A[k][j] / A[k][k]; /* Division step */
        }


        /*
            Spawning nr_of_threads amount of threads and passes the k-value and the thread_id-value to
            give each thread work based on the round-robin principle.
        */
        for (i = 0; i < nr_of_threads; ++i){
            args[i].k = k;
            args[i].thread_id = i;
            pthread_create(&threads[i], NULL, (void *)elimination_step, &(args[i]));
        }

        // Main thread waiting for all the child threads to complete their work.
        // This is needed cause we have certain dependencies through the calculations.
        pthread_barrier_wait(&main_thread_sync_barrier);
    }

    // Free what was allocated, avoiding memory leaks.
    free(args);
    free(threads);
}

void
Init_Matrix()
{
    int i, j;
    printf("1");
    printf("\nsize      = %dx%d ", N, N);
    printf("\nmaxnum    = %d \n", maxnum);
    printf("Init	  = %s \n", Init);
    printf("Initializing matrix...");
    printf("2");
    if (strcmp(Init,"rand") == 0) {
        for (i = 0; i < N; i++){
            for (j = 0; j < N; j++) {
                if (i == j) /* diagonal dominance */
                    A[i][j] = (double)(rand() % maxnum) + 5.0;
                else
                    A[i][j] = (double)(rand() % maxnum) + 1.0;
            }
        }
    }
    if (strcmp(Init,"fast") == 0) {
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                if (i == j) /* diagonal dominance */
                    A[i][j] = 5.0;
                else
                    A[i][j] = 2.0;
            }
        }
    }
    /* Initialize vectors b and y */
    for (i = 0; i < N; i++) {
        b[i] = 2.0;
        y[i] = 1.0;
    }

    printf("done \n\n");
    if (PRINT == 1)
        Print_Matrix();
}

void
Print_Matrix()
{
    int i, j;

    printf("Matrix A:\n");
    for (i = 0; i < N; i++) {
        printf("[");
        for (j = 0; j < N; j++)
            printf(" %5.2f,", A[i][j]);
        printf("]\n");
    }
    printf("Vector b:\n[");
    for (j = 0; j < N; j++)
        printf(" %5.2f,", b[j]);
    printf("]\n");
    printf("Vector y:\n[");
    for (j = 0; j < N; j++)
        printf(" %5.2f,", y[j]);
    printf("]\n");
    printf("\n\n");
}

void
Init_Default()
{
    N = 2048;
    Init = "rand";
    maxnum = 15.0;
    PRINT = 0;
    nr_of_threads = 0;
    v = 1;
}

int
Read_Options(int argc, char **argv)
{
    char    *prog;

    prog = *argv;
    while (++argv, --argc > 0)
        if (**argv == '-')
            switch ( *++*argv ) {
                case 'n':
                    --argc;
                    N = atoi(*++argv);
                    break;
                case 'h':
                    printf("\nHELP: try sor -u \n\n");
                    exit(0);
                    break;
                case 'u':
                    printf("\nUsage: gaussian [-n problemsize]\n");
                    printf("           [-D] show default values \n");
                    printf("           [-h] help \n");
                    printf("           [-I init_type] fast/rand \n");
                    printf("           [-m maxnum] max random no \n");
                    printf("           [-P print_switch] 0/1 \n");
                    printf("           [-t nr_of_threads] number of threads that should be used.\n");
                    printf("           [-v alg_number] which algorithm should be used\n");
                    exit(0);
                    break;
                case 'D':
                    printf("\nDefault:  n         = %d ", N);
                    printf("\n          Init      = rand" );
                    printf("\n          maxnum    = 5 ");
                    printf("\n          P         = 0 ");
                    printf("\n          nr_of_threads = 0");
                    printf("\n          v         = 1\n\n");
                    exit(0);
                    break;
                case 'I':
                    --argc;
                    Init = *++argv;
                    break;
                case 'm':
                    --argc;
                    maxnum = atoi(*++argv);
                    break;
                case 'P':
                    --argc;
                    PRINT = atoi(*++argv);
                    break;
                case 't':
                    --argc;
                    nr_of_threads = atoi(*++argv);
                    break;
                case 'v':
                    --argc;
                    v = atoi(*++argv);
                    break;
                default:
                    printf("%s: ignored option: -%s\n", prog, *argv);
                    printf("HELP: try %s -u \n\n", prog);
                    break;
            }
}
