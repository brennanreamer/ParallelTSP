#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <limits.h>

#define MAX_NUM_CITIES 312

typedef struct { // Create a struct to hold the cost matrix
    int arr[MAX_NUM_CITIES][MAX_NUM_CITIES];
} cities;

double elapsed; // Array to hold the time for each process

int tsp(cities cost, int n, int start, int end) { // Compute the shortest path using the brute force method
    int min_dist = INT_MAX; // initialize the minimum distance to a large number

    double temp_start = MPI_Wtime(); //Start time for parallelized code

    // Compute the shortest distance in each chunk
    int i;
    for (i = start; i < end; i++) { // iterate through each city i within chunk of cities
        int dist = 0;
        int j = (i + 1) % n; // Initialize j to the next element in the array, allowing it to wrap around to the beginning of the array
        while (j != i) {  // search through all cities j that are not city i
            dist += cost.arr[i][j]; // Add the cost from city i to city j to the distance
            j = (j + 1) % n; // Increment j by one, allowing it to wrap around to the beginning of the array
        }

        if (dist < min_dist) { // If the distance is less than the minimum distance
            min_dist = dist; // Set the minimum distance to the distance
        }
    }

    double temp_end = MPI_Wtime(); //end time for parallelized code

    elapsed = temp_end - temp_start;

    return min_dist;
}

int main(int argc, char **argv) {
    // Initialize variables
    int rank, size;
    double start_time;
    cities cost;
    int n;

    //Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    if (rank == 0) { // Rank 0 reads the data from the file and broadcasts it to all processes
        // Open file
        FILE *fp;
        fp = fopen("distances_312.txt", "r");
        if (fp == NULL) {
            printf("Error: File not found.\n");
            return -1;
        }

        // Read number of cities
        fscanf(fp, "%d", &n);

        // Read cost matrix
        int i, j;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                fscanf(fp, "%d", &cost.arr[i][j]);
            }
        }
        fclose(fp);

        start_time = MPI_Wtime();

        
    }
    MPI_Bcast(&cost, sizeof(cities), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Divide the cities into equal-sized chunks for each process
    int start = rank * n / size;
    int end = (rank + 1) * n / size;

    int min_dist = tsp(cost, n, start, end); // Compute the shortest path between all cities
    int tot_min_dist;
    MPI_Reduce(&min_dist, &tot_min_dist, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD); // Combine the results from all processes

    //Find the sum of elapsed time for parallel code
    double sum;
    MPI_Reduce(&elapsed, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // Combine the results from all processes
    
    if (rank == 0) {
        double end_time = MPI_Wtime();
        printf("%f\n", sum/size); // Print the average of the time taken to compute the parallel section of code

        printf("Minimum distance: %d\n", min_dist); // Print the minimum distance
        printf("Time taken: %f seconds\n", end_time - start_time); // Print the total time taken
        
        //printf("Minimum distance\tElapsed time\n");
        //printf("%d\t%f\n", tot_min_dist, end_time - start_time); //Print in good format for copy pasting data
    }

    MPI_Finalize();
    return 0;
}