/***************************************************************************
 *
 * Sequential version of Quick sort
 *
 ***************************************************************************/

#define _XOPEN_SOURCE_600
#define _POSIX_C_SOURCE 200112L

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#define KILO (1024)
#define MEGA (1024*1024)
#define MAX_ITEMS (64*MEGA)
#define swap(v, a, b) {unsigned tmp; tmp=v[a]; v[a]=v[b]; v[b]=tmp;}

int *v;
unsigned int numWorkers = 0;

struct workerArgs {
    unsigned int worker_count;
    int *arr;
    unsigned int low;
    unsigned int high;
};

static void
print_array(void)
{
    int i;
    for (i = 0; i < MAX_ITEMS; i++)
        printf("%d ", v[i]);
    printf("\n");
}

static void
init_array(void)
{
    int i;
    v = (int *) malloc(MAX_ITEMS*sizeof(int));
    for (i = 0; i < MAX_ITEMS; i++)
        v[i] = rand();
}

static unsigned
partition(int *v, unsigned low, unsigned high, unsigned pivot_index)
{
    /* move pivot to the bottom of the vector */
    if (pivot_index != low){
        swap(v, low, pivot_index);
    }
    pivot_index = low;
    low++;

    /* invariant:
     * v[i] for i less than low are less than or equal to pivot
     * v[i] for i greater than high are greater than pivot
     */

    /* move elements into place */
    while (low <= high) {
        if (v[low] <= v[pivot_index])
            low++;
        else if (v[high] > v[pivot_index])
            high--;
        else
            swap(v, low, high);
    }

    /* put pivot back between two groups */
    if (high != pivot_index)
        swap(v, pivot_index, high);
    return high;
}

void * quick_sort(void * orgArgs)
{
    struct workerArgs* args = (struct workerArgs*) orgArgs;
    unsigned pivot_index;
    // Creating two threads that will get one part of the array, a low part and a high part
    // Split in the middle.
    pthread_t worker_one, worker_two;

    /* no need to sort a vector of zero or one element */
    if (args->low >= args->high)
        return NULL;

    /* select the pivot value */
    pivot_index = (args->low+args->high)/2;

    pivot_index = partition(args->arr, args->low, args->high, pivot_index);
    /* partition the vector */

    // Two structs created since the high and low respectivly will be different due to pivot_index -1 or +1.
    struct workerArgs args_one, args_two;
    args_one.arr = args->arr;
    args_one.low = args->low;
    args_one.high = pivot_index-1;
    args_one.worker_count = args->worker_count + 2;


    args_two.arr = args->arr;
    args_two.low = pivot_index+1;
    args_two.high = args->high;
    args_two.worker_count = args->worker_count + 2;

    // Thread created in regards of how many workers were stated in the inparameters
    if (args->worker_count < numWorkers){
        if (args->low < pivot_index){
            pthread_create(&(worker_one), NULL, quick_sort, &args_one);
        }
        if (pivot_index < args->high){
            pthread_create(&(worker_two), NULL, quick_sort, &args_two);
        }
        // Joins to wait for the threads to be done, so that the main thread wont move on.
        pthread_join(worker_one, NULL);
        pthread_join(worker_two, NULL);

    }
    else{


        /* sort the two sub arrays */
        // Sequential sort of the array.
        if (args->low < pivot_index)
            quick_sort(&args_one);
        if (pivot_index < args->high)
            quick_sort(&args_two);
    }

}

// Just a quick function to check if they array is sorted or not.
char* validate_sorted(){
    for(int i = 0; i < MAX_ITEMS - 1; i++){
        if(v[i] > v[i+1]){
            return "False";

        }
    }
    return "True";
}

int
main(int argc, char **argv)
{
    // Inparameter of how many workers should be used.
    if (argc > 1)
        numWorkers = atoi(argv[1]);

    struct workerArgs orgArg;

    init_array();

    // In order to keep track of the data in a easy way and sending it into the thread.
    // This struct is created with the needed data.
    orgArg.worker_count = 0;
    orgArg.low = 0;
    orgArg.high = MAX_ITEMS - 1;
    orgArg.arr = v;

    //print_array();
    quick_sort(&orgArg);

    printf("Is sorted: %s\n", validate_sorted());
    //print_array();
}
