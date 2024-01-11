#include "header.h"


/* UNCOMMENT THE PRINT STATEMENTS FOR DEBUGGING PURPOSES */

int main (int argc, char** argv){

    int size, rank, color, key, provided;
	int nrows, ncols;
	MPI_Comm new_comm;
	double threashold;
    int iterations;
    int sentinel;
    
    //Initializing MPI with multithreaded support 
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if(provided < MPI_THREAD_MULTIPLE){
        printf("The threading support level is lesser than that demanded.\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
	if(rank == 0) {
	    
	    //Prompt the user for the values of rows, columns 
	    printf("Enter the grid dimensions (row col) or enter 0 for default: ");
	    fflush(stdout);
	    scanf("%d", &nrows);
	    
	    if(nrows == 0){
	        nrows=ncols=(int)sqrt(size - 1);
	    }else{
            scanf("%d", &ncols);
        }
        
        printf("Enter the number of iterations for the base station to run: ");
        fflush(stdout);
        scanf("%d", &iterations);
        
        printf("Enter a sentinel value to terminate the programme(must be int): ");
        fflush(stdout);
        scanf("%d", &sentinel);
        
        
        //Fails safe conditions to ensure that the programme runs smoothly
        if( (nrows*ncols) > size - 1) {
		    printf("ERROR: nrows*ncols)=%d * %d = %d > %d. Not enough processes to create nodes. Terminating programme.\n"
		    , nrows, ncols, nrows*ncols, size);
		    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		    return 0;
	    }else if((nrows*ncols) == size){
	        printf("ERROR: nrows*ncols)=%d * %d = %d == %d. Must have at least one more process than the value of row x column. Terminating programme.\n"
		    , nrows, ncols, nrows*ncols, size - 1);
		    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		    return 0;
	    }else{
	        //Prompt the user for the values of sea column height threashold
            printf("Enter a threashold for sea column height between 5500 - 7000: ");
            fflush(stdout);
            scanf("%lf", &threashold);
         }
	    
	    if(threashold < 5500.0 || threashold > 7000.0){
	        printf("Must enter a threashold value between 5500 and 7000!. Terminating programme.\n");
	        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		    return 0;
	    }
	    
	    if(iterations > 1000)
	    {
	        printf("As a fail safe condition the number of iterations must be < 1000. You entered : %d\n", iterations);
	        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
		    return 0;
	    }
		
	
	}
	
	MPI_Bcast(&nrows, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ncols, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&threashold, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&iterations, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&sentinel, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	// n rows x n cols for the sensor nodes. 1 node for the base station. If the user specifies rows x cols s.t rows x cols < total processors - 2
	//There will be unutilized processes
	
	if(rank < nrows * ncols ){              //Sensor node processes

	    color = 0;
	    key = 0;
	
	}else if(rank == size - 1){              //Base station process
	    color = 1;
	    key = size - 1;                                                         //We will always know that worldsize - 1 is the master process
	
	}else{                                  //If we have extra processes we don't use them
	    color = MPI_UNDEFINED;              //We could've terminated them by calling MPI_Finalize() on these extra processes but then 
	    key = rank;                         //during collective communication processes it will trigger a deadlock
	}
	
	//Split the communicator, seperate colours for the base station and the sensor nodes
	MPI_Comm_split(MPI_COMM_WORLD, color, key, &new_comm);
    
    //printf("Rows : %d  Col : %d  From rank: %d\n", nrows, ncols, rank);
    
    if (rank == size - 1)                                                       //Only last process will be the base station            
	    base_station( MPI_COMM_WORLD, new_comm, nrows * ncols, nrows, ncols, threashold, iterations, sentinel);
    
    if(rank < nrows * ncols)
	    sensor_nodes( MPI_COMM_WORLD, new_comm, nrows, ncols, threashold, sentinel);         //The rest of the processes will be the sensor nodes
    
    MPI_Finalize();
    return 0;
}
