#ifndef HEADER_H
#define HEADER_H

//Preprocessor directives
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
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


//Built in data structures
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
    int sentinel;
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



//Funtion prototypes for base station, sensor nodes and satellilte altimeter
int base_station(MPI_Comm world_comm, MPI_Comm comm, int sensor_nodes, int row, int col, int threashold, int iter, int sentinel);
int sensor_nodes(MPI_Comm world_comm, MPI_Comm comm, int rows, int cols, double threashold, int sentinel);               //Slaves will be our tsunameter sensor
float RandomFloat(float min, float max);
int random_int(int val);
void* satellite_altimeter(void* arguments);

//Function to get the ip address of the local machine
void checkIPbuffer(char *IPbuffer);
void checkHostName(int hostname);
void checkHostEntry(struct hostent * hostentry);
int get_ip_addrs(IP_info* ipbuff);

//Funtion to write to file
int writeToFile(Data data, int iterations, char* log_time, double satellite_height, int* reporting_cord, double comm_time, bool match_status, 
int total_matches, int total_mismatches, char* satellite_time);
                     


#endif
