all: run
out:
	mpicc -fopenmp -Wall main.c base_station.c satellite_altimeter.c ip_functions.c sensor_nodes.c write_to_file.c -o out -lm
run:out
	mpirun  -quiet  -oversubscribe  -np 10  out

.PHONY: all run

