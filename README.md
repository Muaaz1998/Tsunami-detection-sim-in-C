# Tsunami-detection-sim-in-C
A simulation of Wireless Sensory Networks using MPI virtual topologies to simulate a tsunami detection system

In order to run the programme please open the terminal --> navigate to the directory where you put this folder and type in the commands:

`$make out`

followed by

`$make run`

**I have also included my tester file "test.c" which I used to generate the number of matches/mismatches of the nodes in the mpi virtual topology(i.e.
the sensor nodes the results of which are included in the report).

** I have not included any Makefile for the "test.c" file. You'll have to compile and run it manually. 

**The makefile has default MPI commands of -oversubscribe -np 10 (10 processes, oversubscribed due to the lack of hardware in my local machine).

**To run more processes please modify the Makefile.

**The log.txt included is for a cartesian grid of 3 x 3 , with a threashold value of 6000 and base station iteration of 30.

**The CSV file "Test.csv" contains the number of matches/mismatches for each sensor node
of MPI cartesian topology created when the WSN is simulated with a cartesian 
grid of 3 x 3 , with a threashold value of 6000 and base station iteration at 30.
