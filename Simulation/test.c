#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <pthread.h>
#include <omp.h>
#include <errno.h>
#include <netdb.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>


//Global constants
#define SIZE 10
#define SHIFT_COL 0
#define SHIFT_ROW 1
#define DISPLACEMENT 1
#define NA 4
#define CORD 2
#define ARRAY_SIZE 10

// basic structure for the logged data

typedef struct{
    double height;
    int coord[2];
    char sa_time[32];
} ARR;

typedef struct{

    ARR* array;
    int threshold;
    int row;
    int col;
    int* term;
    pthread_mutex_t* mutex;
} arg_struct;

typedef struct{
    bool log;
    int matches;
    int my_rank;
    
    char sn_time[32];
    
    int my_cord[CORD];
    double height;
    double na_heights[NA];
    int na_coord[NA][CORD];
    int na_ranks[NA];
    char ip_hostname[32];
    char ip[32];
    int tolerance;
    char na_ip_hostname[4][32];
    char na_ip[4][32];
    
} Data;

typedef struct{
    char hostbuff[32];
    char ipbuff[32];
}IP_info;

typedef struct{
    int coord[2];
    int matches;
    int mismatches;
}Node_matches;

//Funtion prototypes for base station, sensor nodes and satellilte altimeter
int base_station(MPI_Comm world_comm, MPI_Comm comm, int sensor_nodes, int row, int col, int threashold, int iter);
int sensor_nodes(MPI_Comm world_comm, MPI_Comm comm, int rows, int cols, double threashold);               //Slaves will be our tsunameter sensor
float RandomFloat(float min, float max);
int random_int(int val);
void* satellite_altimeter(void* arguments);

//Function to get the ip address of the local machine
void checkIPbuffer(char *IPbuffer);
void checkHostName(int hostname);
void checkHostEntry(struct hostent * hostentry);
int get_ip_addrs(IP_info* ipbuff);

//Funtion to write to file
int writeToFile(Node_matches* node_matches, int sensor_nodes);
                     
//writeToFile(data, iteration, time_logged, satellite_height, reproting_cord, 1, TOL, comm_time, match);

/* UNCOMMENT THE PRINT STATEMENTS FOR DEBUGGING PURPOSES */

int main (int argc, char** argv){

    int size, rank, color, key, provided;
	int nrows, ncols;
	MPI_Comm new_comm;
	double threashold;
    int iterations;
    
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
	    base_station( MPI_COMM_WORLD, new_comm, nrows * ncols, nrows, ncols, threashold, iterations);
    
    if(rank < nrows * ncols)
	    sensor_nodes( MPI_COMM_WORLD, new_comm, nrows, ncols, threashold);         //The rest of the processes will be the sensor nodes
    
    MPI_Finalize();
    return 0;
}

/* This is the base station */
int base_station(MPI_Comm world_comm, MPI_Comm comm, int sensor_nodes, int row, int col, int threashold, int iter)
{   

    int size, worldSize, my_rank, iterations;	
    Data data;                                                  //Receive buffer for the logged data from the sensor nodes
    ARR* array = malloc(10 * sizeof(Data));                                          //Used to measure the communication time
    const int TOL = 100;
    int termination = 0;                                        // termination = 0 we continue, termination = 1 we terminate the sensor nodes
    struct timespec startComp, endComp;
    double comm_time;
    Node_matches* coords = malloc(sensor_nodes * sizeof(Node_matches));
    
    //Initializing MPI                                           
	MPI_Comm_size(world_comm, &worldSize);                      // size of the world communicator
  	MPI_Comm_size(comm, &size);                                 // size of the slave communicator
	MPI_Comm_rank(comm, &my_rank);                              // rank of the master communicator
	MPI_Status recv_status[sensor_nodes];
	
	//Creating the thread that'll simulate the satellite altimeter
	pthread_t thread;

	
	//Initialize the mutex variable
	pthread_mutex_t mutex;
	pthread_mutex_init(&mutex, NULL);
	
	//Creating the structure to be passed on as argument to our thread function
    arg_struct args; 
    args.array = array;
    args.threshold = threashold;
    args.row = row;
    args.col = col;
    args.term = &termination;
    args.mutex = &mutex;

	//Creating the data structure for the logged data
	MPI_Datatype MPI_log_type;
    int lengths[14] = {1, 1, 1, 32, 2, 1, 4, 8, 4, 32, 32, 1, 128, 128};
    
    MPI_Aint displacements[14];
    Data dummy_data;
    MPI_Aint base_address;
    MPI_Get_address(&dummy_data, &base_address);
    
     
    //Get the addresses of the structure object's variables in memory
    MPI_Get_address(&dummy_data.log, &displacements[0]);
    MPI_Get_address(&dummy_data.matches, &displacements[1]);
    MPI_Get_address(&dummy_data.my_rank, &displacements[2]);
    MPI_Get_address(&dummy_data.sn_time, &displacements[3]);
    MPI_Get_address(&dummy_data.my_cord, &displacements[4]);
    MPI_Get_address(&dummy_data.height, &displacements[5]);
    MPI_Get_address(&dummy_data.na_heights, &displacements[6]);
    MPI_Get_address(&dummy_data.na_coord, &displacements[7]);
    MPI_Get_address(&dummy_data.na_ranks, &displacements[8]);
    MPI_Get_address(&dummy_data.ip_hostname, &displacements[9]);
    MPI_Get_address(&dummy_data.ip, &displacements[10]);
    MPI_Get_address(&dummy_data.tolerance, &displacements[11]);
    MPI_Get_address(&dummy_data.na_ip_hostname, &displacements[12]);
    MPI_Get_address(&dummy_data.na_ip, &displacements[13]);

    for(int i = 0 ; i < 14; i++){
        displacements[i] = MPI_Aint_diff(displacements[i], base_address);
    }
    
    //Now we initialize the structure for the MPI
    MPI_Datatype types[14] = {MPI_C_BOOL , MPI_INT, MPI_INT, MPI_CHAR, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT ,MPI_CHAR, MPI_CHAR, MPI_INT, MPI_CHAR, MPI_CHAR};                              
    MPI_Type_create_struct(14, lengths, displacements, types, &MPI_log_type);    //Now we create our MPI structure
    MPI_Type_commit(&MPI_log_type);
    
    //Base station will receive the co-ordinates of the sensor nodes and store it
    MPI_Status recv[sensor_nodes];
    
    for(int i = 0; i < sensor_nodes; i++){
        MPI_Recv(coords[i].coord, 2, MPI_INT, i, 0, world_comm, &recv[i]);
        coords[i].matches = 0;
        coords[i].mismatches = 0;
    
    }
    

	//Start the satellite altimeter
	if(pthread_create(&thread, NULL, &satellite_altimeter, (void*)&args) != 0){
        perror("Failed to create thread!\n");
    }

     //Base station will run for a fixed number of iterations specified by the user
	// Until we reach the last iteration the base station will keep sending non-terminating messages to the sensor nodes
	// Once we hit the last iteration we will send a terminating message to each of the sensor nodes and satellite altimeter to terminate
	for(iterations = 0 ; iterations < iter; iterations++){
	    
	    if(iterations < iter - 1){
	        //Tell the satellite altimeter to start executing
	    
	        //At each iteration of the base station it'll tell the nodes to carry on executing 
	        for(int j = 0 ; j < sensor_nodes ; j++){
	            MPI_Send(&termination, 1, MPI_INT, j, 0, world_comm);
	        }
	        
	        
	        //The base station has to check for incoming log messages
	        for(int j = 0 ; j < sensor_nodes ; j++){
	            
	            //Measure the communication time between successive send/recv calls between the sensor nodes and base station
	            //Start time
	            clock_gettime(CLOCK_MONOTONIC, &startComp);
	            
	            //receive logged message from the sensor nodes
	            MPI_Recv(&data, 1, MPI_log_type, j, MPI_ANY_TAG, world_comm, &recv_status[j]); 
	        	
	        	//End time
	        	clock_gettime(CLOCK_MONOTONIC, &endComp);                           
	        	
	        	//Total communication time to the nearest nanosecond
	        	comm_time = (endComp.tv_sec - startComp.tv_sec) * 1e9;
                comm_time  = (comm_time + (endComp.tv_nsec - startComp.tv_nsec)) * 1e-9;                                     //This is our communication time
	        	
	        	
	        	//Now we need to match the sending co-ordinates of the data and the height(within a given tolerance) with the satellite's
	        	//data
	            if(data.log){	
	                //we initialize a mutex lock so that when reading from the data array(which is simultaneously being written to by the satellite)
	                //we don't trigger any race conditions or illegal reads to memory
	            	pthread_mutex_lock(&mutex);                              //lock mutex to prevent race conditions          
                    for(int i = 0 ; i < ARRAY_SIZE ; i++){ // h 6178 s 5563
                        if( (array[i].coord[0] == data.my_cord[0]) && (array[i].coord[1] == data.my_cord[1]) ){
                            if((abs(array[i].height - data.height) <= TOL)){
                                for(int k = 0; k < sensor_nodes; k++){
                                    if(data.my_cord[0] == coords[k].coord[0] && data.my_cord[1] == coords[k].coord[1]){
                                        coords[k].matches++;
                                        break;
                                    }
                                }
                            }else{
                                 for(int k = 0; k < sensor_nodes; k++){
                                    if(data.my_cord[0] == coords[k].coord[0] && data.my_cord[1] == coords[k].coord[1]){
                                        coords[k].mismatches++;
                                        break;
                                 }
                                    
                                }
                            }
                        }
                    }
                   pthread_mutex_unlock(&mutex);
	            }
	        }
	    }else{
	        //Terminate the base station and satellilte altimeter gracefully
	        
	        termination = 1;
	        
	        //Join the thread once the array has been initialized 
	        if((pthread_join(thread,NULL)) != 0){
                perror("Failed to join thread\n");
            } 
	        
	        //Terminate the sensor nodes gracefully
	        for(int j = 0 ; j < sensor_nodes ; j++){
	            MPI_Send(&termination, 1, MPI_INT, j, 0, world_comm);
	        }
	    }   
	}
    
    //Now we write to file
    writeToFile(coords, sensor_nodes);
    
    //Destroy the mutex and free memory allocated for dynamic array
    pthread_mutex_destroy(&mutex);
    free(array);
	MPI_Type_free(&MPI_log_type);
    MPI_Comm_free( &comm);
    return 0;
}


/* 
This is the sensor nodes 
They can also receive messages from the base station
*/
int sensor_nodes(MPI_Comm world_comm, MPI_Comm comm, int rows, int cols, double threashold)
{
	int ndims=2, size, my_rank, reorder, my_cart_rank, ierr, worldSize, my_world_rank, err;
	int dims[ndims],coord[ndims];
	int wrap_around[ndims];
    int termination;                                            //Sent by base station. Dictates whether or not we terminate the process
    int heights[SIZE];                                            //Where we will store our heights for the moving average
    int counter = 0;                                            //Keeps track of how many elements we have on our moving average
    double moving_average = 0.0;
    int neighbours[4];
    enum DIRECTIONS {DOWN, UP, LEFT, RIGHT};
    const int num_neigbours = 4;
    double na_moving_avg[NA];                               //Stores neighbours receinved moving average. Each index corresponds to the moving 
    int matches = 0;                                                             //Keep track of the number of matches we got
    int na_cord[NA][CORD];
    const int tolerance = 500;
    const double lower = threashold - tolerance;
    const double upper = threashold + tolerance;
    int log_tag = 0;                                     //Will use as a tag when sending log messages to the base station
    Data data;
    char na_ip_name[4][32];
    char na_ip[4][32];
    int i;
    IP_info ip;

    
    //Array to store the matches of each neighbour
    int na_matches[NA]; 
             
    //MPI_calls
    MPI_Comm comm2D;
    MPI_Status status;
    MPI_Comm_size(world_comm, &worldSize);                              // size of the world communicator
  	MPI_Comm_size(comm, &size);                                         // size of the slave communicator
	MPI_Comm_rank(comm, &my_rank);                                      // rank of the slave communicator
	MPI_Comm_rank(world_comm, &my_world_rank);
	
	//Seed for the random number generator, based on the rank of each individual sensor nodes
	srand(pow(my_world_rank, 7));
    
    //Get the ip address and hostname
    get_ip_addrs(&ip);
	
	//Initializing the dimensions
	dims[0] = rows;
	dims[1] = cols;
    
    //Creating our virtual topology
	err = MPI_Dims_create(size, ndims, dims);                               //If the dimensions are set to 0 MPI_Dims_create will select the best dimensions based
	if(err != MPI_SUCCESS){
	    printf("Couldn't create virtual topology");
	    return 2;
	}            
	                                                           
    /* create cartesian mapping */
	wrap_around[0] = 0;
	wrap_around[1] = 0; /* periodic shift is .false. */
	reorder = 0;
	ierr = 0;
	ierr = MPI_Cart_create(comm, ndims, dims, wrap_around, reorder, &comm2D);
	
	if(ierr != 0){
	    printf("ERROR[%d] creating CART\n",ierr);
	}
	
	/* find my coordinates in the cartesian communicator group */
	MPI_Cart_coords(comm2D, my_rank, ndims, coord);                                                         // coordinated is returned into the coord array
	
	/* use my cartesian coordinates to find my rank in cartesian group*/
	MPI_Cart_rank(comm2D, coord, &my_cart_rank);
	
	
	//Get the rank of each of the neighbours
	MPI_Cart_shift( comm2D, SHIFT_COL, DISPLACEMENT, &neighbours[LEFT], &neighbours[RIGHT]);
	MPI_Cart_shift( comm2D, SHIFT_ROW, DISPLACEMENT, &neighbours[DOWN], &neighbours[UP]);
	
	
	//Each MPI process will send it's coordinates to the base
	MPI_Send(coord, 2, MPI_INT, worldSize - 1, 0, world_comm);
	
	//Get the neighbours co-ordinates and store it
	
	#pragma omp parallel for default(shared) private(i) schedule(static, 1) 
	for(i = 0 ; i < num_neigbours ; i++){
	    if(neighbours[i] != MPI_PROC_NULL){
	        MPI_Cart_coords(comm2D, neighbours[i], ndims, na_cord[i]);                                      // 0 = down, 1 = up, 2 = left, 3 = right
	    }                                  
	    else{
	        na_cord[i][0] = -1;
	        na_cord[i][1] = -1;
	    }             
	}
    
    MPI_Status receive_status;
    
    //get the ipv4 address of each of the neighbours
    #pragma omp parallel for default(shared) private(i) schedule(static, 1)
    for(i = DOWN; i<=RIGHT; i++){
        if(neighbours[i] != MPI_PROC_NULL){
            MPI_Send(ip.hostbuff, 32, MPI_CHAR, neighbours[i], 0, comm2D);
            MPI_Recv(na_ip_name[i], 32, MPI_CHAR, neighbours[i], 0, comm2D, &receive_status);
        }else{
            strcpy(na_ip[i], "Null");
        }
    }	
    
    
    #pragma omp parallel for default(shared) private(i) schedule(static, 1)
    for(i = DOWN; i<=RIGHT; i++){
        if(neighbours[i] != MPI_PROC_NULL){
            MPI_Send(ip.ipbuff, 32, MPI_CHAR, neighbours[i], 0, comm2D);
            MPI_Recv(na_ip[i], 32, MPI_CHAR, neighbours[i], 0, comm2D, &receive_status);
        }else{
            strcpy(na_ip[i], "Null");
        }
    }	


	
	//Creating the data structure for the logged data
	MPI_Datatype MPI_log_type;
    int lengths[14] = {1, 1, 1, 32, 2, 1, 4, 8, 4, 32, 32, 1, 128, 128};
    
    MPI_Aint displacements[14];
    Data dummy_data;
    MPI_Aint base_address;
    MPI_Get_address(&dummy_data, &base_address);
    
     
    //Get the addresses of the structure object's variables in memory
    MPI_Get_address(&dummy_data.log, &displacements[0]);
    MPI_Get_address(&dummy_data.matches, &displacements[1]);
    MPI_Get_address(&dummy_data.my_rank, &displacements[2]);
    MPI_Get_address(&dummy_data.sn_time, &displacements[3]);
    MPI_Get_address(&dummy_data.my_cord, &displacements[4]);
    MPI_Get_address(&dummy_data.height, &displacements[5]);
    MPI_Get_address(&dummy_data.na_heights, &displacements[6]);
    MPI_Get_address(&dummy_data.na_coord, &displacements[7]);
    MPI_Get_address(&dummy_data.na_ranks, &displacements[8]);
    MPI_Get_address(&dummy_data.ip_hostname, &displacements[9]);
    MPI_Get_address(&dummy_data.ip, &displacements[10]);
    MPI_Get_address(&dummy_data.tolerance, &displacements[11]);
    MPI_Get_address(&dummy_data.na_ip_hostname, &displacements[12]);
    MPI_Get_address(&dummy_data.na_ip, &displacements[13]);

    for(int i = 0 ; i < 14; i++){
        displacements[i] = MPI_Aint_diff(displacements[i], base_address);
    }
    
    //Now we initialize the structure for the MPI
    MPI_Datatype types[14] = {MPI_C_BOOL , MPI_INT, MPI_INT, MPI_CHAR, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT ,MPI_CHAR, MPI_CHAR, MPI_INT, MPI_CHAR, MPI_CHAR};                              
    MPI_Type_create_struct(14, lengths, displacements, types, &MPI_log_type);    //Now we create our MPI structure
    MPI_Type_commit(&MPI_log_type);

	//The MPI process will terminate upon receiving a termination message from the base station / master node
	do{  
	   //printf("\nReceiving terminating message\n");
       MPI_Recv(&termination, 1, MPI_INT, worldSize - 1, MPI_ANY_TAG, world_comm, &status);
       //printf("Received terminating message\n");
       //What do we do when we don't terminate
       if(termination == 0)
       {    
            //printf("\nCurrent process: %d\n", my_rank);
            //Generate random float values after a specified time, to simulate the sensor reading
            
            
            float h = RandomFloat(upper, lower);
            heights[counter] = h;
            counter++;
            //printf("height : %.2f\n", h);
            
            //When counter == size we compute the moving average 
            if(counter == SIZE){
                moving_average = 0;
                for(int i = 0; i < SIZE ; i++){
                    moving_average += heights[i];
                }
                
                moving_average = moving_average / SIZE;
                counter = 0;
            }
            //printf("checking if Moving average > threashold\n Moving average of process %d is %.2f\n", my_cart_rank, moving_average);
            
            if(moving_average > threashold){
                //printf("Moving average > threashold\n");
                //Get readings from neighbour nodes
         
                //Asynchronous communication
                
               	MPI_Request send_request[4];
                MPI_Request receive_request[4];
                MPI_Status send_status[4];
                MPI_Status receive_status[4];
                
                #pragma omp parallel for default(shared) schedule(static, 1) private(i)
                for(i = DOWN ; i <= RIGHT ; i++){
                    MPI_Isend(&moving_average, 1, MPI_DOUBLE, neighbours[i], 0, comm2D, &send_request[i]);
                    MPI_Irecv(&na_moving_avg[i], 1, MPI_DOUBLE, neighbours[i], 0, comm2D, &receive_request[i]);
                }
                
                MPI_Waitall(4, send_request, send_status);
	            MPI_Waitall(4, receive_request, receive_status);
                //Checking which p2p communication triggers a deadlock
                
                #pragma omp parallel for default(shared) schedule(static, 1) private(i)
                for(i = DOWN; i <= RIGHT; i++){
                    if(neighbours[i] < 0){
                        na_moving_avg[i] = 0;
                    }
                }
                
                //printf("Not triggering deadlock at line 427");
      
                //Now we need to check neighbour's moving average for a potential matches
                matches = 0;
                
                #pragma omp parallel for default(shared) private(i) schedule(static, 2) reduction(+:matches) 
                for(i = DOWN ; i <= RIGHT ; i++){
                    if(abs(na_moving_avg[i] - moving_average) <= 500){
                        matches += 1;
                    }
                }

                /*
                1. Check if the number of matches > 1
                2. If the number of matches is > 1 each of the neighbouring processes will send it's number of matches along with it's rank
                to each of it's neighbours
                3. We then check the rank of the process which had the most matches(that'll be the process that consists of the most information)
                4. Only that process will send the log details to the base station by creating a structure consisting of all the details.
                */
                
                if(matches >= 2){
                    
                    MPI_Request send_request_2[4];
                    MPI_Request receive_request_2[4];
                    MPI_Status send_status_2[4];
                    MPI_Status receive_status_2[4];
                    
                    #pragma omp parallel for default(shared) schedule(static, 1) private(i)
                    for(i = DOWN ; i <= RIGHT ; i++){
                        MPI_Isend(&matches, 1, MPI_INT, neighbours[i], 1, comm2D, &send_request_2[i]);
                        MPI_Irecv(&na_matches[i], 1, MPI_INT, neighbours[i], 1, comm2D, &receive_request_2[i]);
                    }                
                    

                    
                    MPI_Waitall(4, send_request, send_status_2);
	                MPI_Waitall(4, receive_request, receive_status_2);
	                
	                #pragma omp parallel for default(shared) private(i) schedule(static, 1)
                    for(i = DOWN; i <= RIGHT; i++){
                        if(neighbours[i] < 0){
                            na_matches[i] = 0;
                        }
                    }	            
	               
                   
                    //Details to include : Time of alert, reporting node coordinates, reporting node height, ipv4 address, adjacent node co-ordinates,
                    // height and ipv4 address,communication time, sending node ranks and neighbours ranks   
                    //First we need to get the current time and then create the structure and send all relevant values
                    
                    //get the current time to notify the base station of the time of occurence of the alert
                    time_t  current_time;
                    time(&current_time);

                    //Properties of our structure
                    data.log = true;                                                    //The data sent by this node will need to be logged
                    data.matches = matches;                                             //matches 
                    data.my_rank = my_rank;
                    
                    //Time details
                    strcpy(data.sn_time, ctime(&current_time));
                    
                    //Also need to send the rank of the neighbouring processes with respect to the world communicator  
                    memcpy(data.my_cord, coord, sizeof(coord));                                 //Co-ordinates of reporting rank
                    data.height = moving_average;                                               //Moving average of sending node
                    memcpy(data.na_heights, na_moving_avg, sizeof(na_moving_avg));              //Moving average of neighbours
                    memcpy(data.na_coord, na_cord, sizeof(na_cord));                            //Neighbours co-ordinates 
                    strcpy(data.ip_hostname, ip.hostbuff);
                    strcpy(data.ip, ip.ipbuff);                                                      //ip address
                    memcpy(data.na_ranks, neighbours, sizeof(neighbours));
                    data.tolerance  = tolerance; 
                    memcpy(data.na_ip_hostname, na_ip_name, sizeof(na_ip_name));
                    memcpy(data.na_ip, na_ip, sizeof(na_ip));
                }
            }else{
                 
                 //printf("If moving average < threasold\n");
                 
                 //if moving average < threashold
                 //Then we essentially want to tell the base station that we have nothing to report
                 //Also need to simulate a sending and receiving of the moving average so that it doesn't trigger a deadlock
                 //Also need to receive the number of matches from our neighbours
                                                    
                //printf("Reached line 571\n");
                
             	MPI_Request send_request[4];
                MPI_Request receive_request[4];
                MPI_Status send_status[4];
                MPI_Status receive_status[4];
                 
                
                #pragma omp parallel for default(shared) private(i) schedule(static , 1)
                for(i = DOWN ; i <= RIGHT ; i++){
                    MPI_Isend(&moving_average, 1, MPI_DOUBLE, neighbours[i], 0, comm2D, &send_request[i]);
                    MPI_Irecv(&na_moving_avg[i], 1, MPI_DOUBLE, neighbours[i], 0, comm2D, &receive_request[i]);
                }

                MPI_Waitall(4, send_request, send_status);
	            MPI_Waitall(4, receive_request, receive_status);
                
                #pragma omp parallel for default(shared) private(i) schedule(static, 1)
                for(i = DOWN; i <= RIGHT; i++){
                    if(neighbours[i] < 0){
                        na_moving_avg[i] = 0;
                    }
                }
                      
                //printf("Reached line 592\n");

                MPI_Request send_request_2[4];
                MPI_Request receive_request_2[4];
                MPI_Status send_status_2[4];
                MPI_Status receive_status_2[4];
                
                #pragma omp parallel for default(shared) private(i) schedule(static, 1)
                for(i = DOWN ; i <= RIGHT ; i++){
                    MPI_Isend(&matches, 1, MPI_INT, neighbours[i], 1, comm2D, &send_request_2[i]);
                    MPI_Irecv(&na_matches[i], 1, MPI_INT, neighbours[i], 1, comm2D, &receive_request_2[i]);
                }
                
                MPI_Waitall(4, send_request, send_status_2);
	            MPI_Waitall(4, receive_request, receive_status_2);
                
                #pragma omp parallel for default(shared) private(i) schedule(static, 1)
                for(i = DOWN; i <= RIGHT; i++){
                    if(neighbours[i] < 0){
                        na_matches[i] = 0;
                    }
                }
                	     
               data.log = false;                                  //A way of telling the base station would be that we have had no matches 
            
            
            }
        
            //One MPI call to the base station reporting the data to be logged(if any mathces have been found)
            //printf("Sending to base station\n");
            MPI_Barrier(comm2D);
            MPI_Send(&data, 1, MPI_log_type, worldSize - 1, log_tag, world_comm);
            log_tag++;
            //printf("Done sending to base station data\n");
       }else{
            //We have received a terminating message from the base station
            printf("Received a terminating message from the base station. The node with rank %d will now terminate.\n", my_rank);
       }
	}while(termination == 0);
	

    //free(na_matches);
    MPI_Type_free(&MPI_log_type);
    MPI_Comm_free( &comm2D );
    return 0;
	
}

void* satellite_altimeter(void* arguments){ 

    //Need to dereference the structure and it's variables
    arg_struct *args = arguments;
    int row = args->row;
    int col = args->col;
    int i = -1;                          //We will use this as a reference to the index of the array which we will be populating
    
    //Seed for rng
    srand(time(NULL));
    
    //Get the time of day when the sea water column height was measured
    time_t  current_time;
    time(&current_time);
    
    //Until it receives a terminating condition from the base station the satellite altimter will keep on populating the global array
    while(*(args -> term) != 1){
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

//Funtion to generate random float values
float RandomFloat(float max, float min){
   return ((max - min) * ((float)rand() / RAND_MAX)) + min;
}

int random_int(int val)
{
    return rand() % val;
}



// Returns hostname for the local computer
void checkHostName(int hostname)
{
    if (hostname == -1)
    {
        perror("gethostname");
        exit(1);
    }
}
  
// Returns host information corresponding to host name
void checkHostEntry(struct hostent * hostentry)
{
    if (hostentry == NULL)
    {
        perror("gethostbyname");
        exit(1);
    }
}
  
// Converts space-delimited IPv4 addresses
// to dotted-decimal format
void checkIPbuffer(char *IPbuffer)
{
    if (NULL == IPbuffer)
    {
        perror("inet_ntoa");
        exit(1);
    }
}
  
// Driver code
int get_ip_addrs(IP_info* ip_buff)
{
    char hostbuffer[32];
    char *IPbuffer;
    struct hostent *host_entry;
    int hostname;
  
    // To retrieve hostname
    hostname = gethostname(hostbuffer, sizeof(hostbuffer));
    checkHostName(hostname);
  
    // To retrieve host information
    host_entry = gethostbyname(hostbuffer);
    checkHostEntry(host_entry);
  
    // To convert an Internet network
    // address into ASCII string
    IPbuffer = inet_ntoa(*((struct in_addr*)
                           host_entry->h_addr_list[0]));
  
    strcpy(ip_buff->hostbuff, hostbuffer);
    strcpy(ip_buff->ipbuff, IPbuffer);
  
    return 0;
}


int writeToFile(Node_matches *node_matches, int sensor_nodes)
{
    FILE* log = fopen("Test.csv", "w");
    fprintf(log,"Coordinates,Matches,Mismatches\n");
    for(int i = 0; i < sensor_nodes; i++){
        fprintf(log, "(%d %d),%d,%d\n", node_matches[i].coord[0], node_matches[i].coord[1], node_matches[i].matches, node_matches[i].mismatches);
    }
    fclose(log);
                    
    return 1;                       
}


