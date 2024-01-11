#include "header.h"

/*This is the sensor node*/

int sensor_nodes(MPI_Comm world_comm, MPI_Comm comm, int rows, int cols, double threashold, int sentinel)
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
	                
	                /*
	                printf("Process %d\n", my_cart_rank);
                    for(int i = 0 ; i < num_neigbours; i++){
                        if(neighbours[i] >= 0){
                            printf("Neighbour %d Moving average : %.2f and matches : %d\n", i, na_moving_avg[i], na_matches[i]);
                        }
                    
                    }
	                */
	                
               
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
                     /*
                    for(int i = DOWN; i <= RIGHT ; i++){
                        printf("Data coordinates : (%d, %d)\nHeight: %lf\n", data.na_coord[i][0], data.na_coord[i][1], data.na_heights[i]);
                    }
                    */

                }
              
            }else{
                 
                 
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
	}while(termination != sentinel);
	

    //free(na_matches);
    MPI_Type_free(&MPI_log_type);
    MPI_Comm_free( &comm2D );
    return 0;
	
}

//Funtion to generate random float values
float RandomFloat(float max, float min){
   return ((max - min) * ((float)rand() / RAND_MAX)) + min;
}

int random_int(int val)
{
    return rand() % val;
}


