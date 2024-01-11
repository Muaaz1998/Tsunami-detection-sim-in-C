#include "header.h"

/* This is the base station */
int base_station(MPI_Comm world_comm, MPI_Comm comm, int sensor_nodes, int row, int col, int threashold, int iter, int sentinel)
{   

    int size, worldSize, my_rank, iterations;	
    Data data;                                                  //Receive buffer for the logged data from the sensor nodes
    ARR* array = malloc(10 * sizeof(Data));                                          //Used to measure the communication time
    bool match = false;
    bool mismatch = false;
    const int TOL = 100;
    double satellite_height;
    int reporting_cord[2];
    int total_matches = 0;
    int total_mismatches = 0;
    char sa_time[32];
    int termination = 0;                                        // termination = 0 we continue, termination = 1 we terminate the sensor nodes
    struct timespec startComp, endComp;
    double comm_time;
    
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
    args.sentinel = sentinel;

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
                            memcpy(reporting_cord, array[i].coord, sizeof(array[i].coord));
                            if((abs(array[i].height - data.height) <= TOL)){
                                match = true;
                                total_matches++;
                                strcpy(sa_time, array[i].sa_time);
                                //printf("Array height: %.2f\n", array[i].height);
                                //printf("Data height ; %.2f\n", data.height);
                                satellite_height = array[i].height;
                                break;
                            }
                            else{
                                mismatch = true;
                                total_mismatches++;
                                strcpy(sa_time, array[i].sa_time);
                                //printf("Array height: %.2f\n", array[i].height);
                                //printf("Data height ; %.2f\n", data.height);
                                satellite_height = array[i].height;
                                break;
                            }
                        }
                    }
	            	pthread_mutex_unlock(&mutex);
	                //After receiving we check the log status of the data
	               if(match || mismatch){
	                    // Writing log details to file
	                    time_t  current_time;
	                    time(&current_time);
	                    char log_time[32];
	                    strcpy(log_time , ctime(&current_time));
	                    int res;
	                    if(match){
	                        res = writeToFile(data, iterations, log_time, satellite_height, reporting_cord, comm_time, true, total_matches, total_mismatches, sa_time);
	                        printf("Match!\n");
	                        }
	                    if(mismatch){   
	                        res = writeToFile(data, iterations, log_time, satellite_height, reporting_cord, comm_time, false, total_matches, total_mismatches,sa_time);
	                        printf("Mismatch!\n");
	                        }
                        if(res == 1)
                            printf("Written to file successfully. Total number of log messages %d\n", (total_matches + total_mismatches));
                        
                        match = false;
                        mismatch = false;
	                
	                }
                }
	        }
	    }else{
	        //Terminate the base station and satellilte altimeter gracefully
	        
	        termination = sentinel;
	        
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
    
    //Destroy the mutex and free memory allocated for dynamic array
    pthread_mutex_destroy(&mutex);
    free(array);
	MPI_Type_free(&MPI_log_type);
    MPI_Comm_free( &comm);
    return 0;
}

