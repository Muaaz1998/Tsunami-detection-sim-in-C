#include "header.h"

void* satellite_altimeter(void* arguments){ 

    //Need to dereference the structure and it's variables
    arg_struct *args = arguments;
    int row = args->row;
    int col = args->col;
    int i = -1;                          //We will use this as a reference to the index of the array which we will be populating
    int sentinel = args->sentinel;
    
    //Seed for rng
    srand(time(NULL));
    
    //Get the time of day when the sea water column height was measured
    time_t  current_time;
    time(&current_time);
    
    //Until it receives a terminating condition from the base station the satellite altimter will keep on populating the global array
    while(*(args -> term) != sentinel){
        //Lock the mutex so that when the base station is reading from the data array and the satellite is simultaneously writing to the data array
        // we don't trigger any deadlocks or perform a read operation on an uninitialized area of memory
        pthread_mutex_lock(args->mutex);
        //Increment the index
        // By performing a modulo operation over the index by the size of the array we can simulate a FIFO behabiour
        i = (i + 1) % ARRAY_SIZE;
        
        args->array[i].height = RandomFloat(args->threshold + 500, args->threshold);
        args->array[i].coord[0] =  random_int(row);
        args->array[i].coord[1] = random_int(col);
        strcpy(args->array[i].sa_time, ctime(&current_time));
        pthread_mutex_unlock(args->mutex);
    }
    
    return NULL;

}
