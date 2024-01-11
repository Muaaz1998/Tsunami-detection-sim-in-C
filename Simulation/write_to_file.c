#include "header.h"

//Writes the base station's logged data into a file called "log.txt" containing the logged data

int writeToFile(Data data, int iterations, char* log_time, double satellite_height, int* reporting_cord,double comm_time, bool match_status, int total_matches, int total_mismatches, char* sa_time)
{
    FILE* log;
    //When writing to the file for the first time we want to just open the file and write
    if(total_matches + total_mismatches == 1){
        log = fopen("log.txt", "w");
    }else{
       //For the rest of the writing we want to append to the file
       log = fopen("log.txt", "ab");
    } 
    fprintf(log, "Iteration: %d\n\n", iterations);
    fprintf(log, "Logged time:\t\t\t %s", log_time);
    fprintf(log, "Alert reported time:\t %s", data.sn_time);
    
    if(match_status)                //match status == true implies a match, match status == false implies a mismatch
        fprintf(log, "Alert type: \t\t\t Match(or True)\n\n");
    else
        fprintf(log, "Alert type: \t\t\t Mismatch(or False)\n\n");
        
    fprintf(log, "Reporting Node\t\tCoord\t\t\t\tHeight\t\t\t\tip\t\t\t\tHostname\n");
    fprintf(log, "%d\t\t\t\t\t(%d, %d)\t\t\t\t%.2f\t\t\t%s\t\t\t%s\n\n", data.my_rank, data.my_cord[0], data.my_cord[1], data.height, data.ip, data.ip_hostname);
    fprintf(log, "Adjacent Nodes\t\tCoord\t\t\t\tHeight\t\t\t\tip\t\t\t\tHostname\n");
    for(int i = 0 ; i < 4; i++){
        if(data.na_ranks[i] >= 0)
            fprintf(log, "%d\t\t\t\t\t(%d, %d)\t\t\t\t%.2f\t\t\t%s\t\t\t%s\n", data.na_ranks[i], data.na_coord[i][0], data.na_coord[i][1], data.na_heights[i], data.na_ip[i], data.na_ip_hostname[i]);
         if(i == 3)
            fprintf(log,"\n");
    }    
    fprintf(log, "Satellite altimeter reporting time: %s", sa_time);
    fprintf(log, "Satellite altimeter reporting height : %.2f\n", satellite_height);
    fprintf(log, "Satellite altimeter reporting Cord : (%d, %d)\n", reporting_cord[0], reporting_cord[1]);
    fprintf(log, "Communication time: %.4fs\n", comm_time);
    fprintf(log, "Number of adjacent matches to reporting node: %d\n", data.matches);
    fprintf(log, "Max. tolerance range between nodes readings (m): 500\n");
    fprintf(log, "Max. tolerance range between satellite altimeter and reporting node readings (m): 100\n");
    fprintf(log, "Total matches so far: %d / %d\n", total_matches, (total_matches + total_mismatches));
    fprintf(log, "Total mismatches so far: %d / %d\n", total_mismatches, (total_matches + total_mismatches));
    fprintf(log, "------------------------------------------------------------------------------------------------------------\n");
    fprintf(log, "------------------------------------------------------------------------------------------------------------\n");
    
    fclose(log);
                    
    return 1;                       
}
