#ifndef _MAIN_C_
#define _MAIN_C_

#include "LBDefinitions.h"
#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "initLB.h"

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>

#define sign(x) ((x > 0) - (x < 0))

/* Calculate timespec difference */
struct timespec ts_diff(struct timespec a, struct timespec b) {
	a.tv_sec = a.tv_sec - b.tv_sec;
	a.tv_nsec = a.tv_nsec - b.tv_nsec;

	a.tv_sec = abs(a.tv_sec) - 1 * ((sign(a.tv_sec) * sign(a.tv_nsec)) < 0);
	a.tv_nsec = abs(1000000000 * ((sign(a.tv_sec) * sign(a.tv_nsec)) < 0) - abs(a.tv_nsec));

	return a;
}

/* Calculate double value (as decimal seconds) */
double ts_to_double(struct timespec time) {
	return time.tv_sec + time.tv_nsec / 10e9;
}

int main(int argc, char *argv[]) {
	double *collideField = NULL;
	double *streamField = NULL;
	int *flagField = NULL;
	int xlength;
	double tau;
	double velocityWall[D];
	int timesteps;
	int timestepsPerPlotting;
    struct timespec begin, end, mb, me;

    /* Start measuring total execution time */
    clock_gettime(CLOCK_REALTIME, &begin);

    printf("LBM simulation by Krivokapic, Mody, Malcher - CFD Lab SS2014\n");
    printf("============================================================\n");
    printf("Reading in parameters...\n");
	if (readParameters(&xlength, &tau, velocityWall, &timesteps,
			&timestepsPerPlotting, argc, argv)) {
		printf("Reading in parameters failed. Aborting program!\n");
		exit(-1);
	}
    printf("...done\n");
    
    printf("Starting LBM...\n");
    /* (xlength+2)^D elements must be stored for all lattices including boundaries */
    const int xl_to3 = (xlength + 2) * (xlength + 2) * (xlength + 2);
    collideField = malloc(sizeof *collideField * Q * xl_to3);
    streamField = malloc(sizeof *streamField * Q * xl_to3);
    flagField = malloc(sizeof *flagField * xl_to3);
    
    /* Initialize pointers */
	initialiseFields(collideField, streamField, flagField, xlength);

    double mlups = 0.0;
    double mlups_time = 0.0;
    
    for (int t = 0; t < timesteps; ++t) {
		double *swap = NULL;
        /* Start measuring MLUPS here */
        clock_gettime(CLOCK_REALTIME, &mb);
		doStreaming(collideField, streamField, flagField, xlength);
		swap = collideField;
		collideField = streamField;
		streamField = swap;

		doCollision(collideField, flagField, &tau, xlength);
		treatBoundary(collideField, flagField, velocityWall, xlength);
        
        /* Stop measuring MLUPS here (for current timestep) */
        clock_gettime(CLOCK_REALTIME, &me);
        mlups_time += ts_to_double(ts_diff(mb, me));
        mlups += xl_to3/1000000.0;
        /* Write output to vtk file for postprocessing */
		if (t % timestepsPerPlotting == 0) {
			writeVtkOutput(collideField, flagField, "lbm_out", t, xlength);
		}
	}
    clock_gettime(CLOCK_REALTIME, &end);
    printf("============================================================\n");
    printf("LBM simulation completed for %d cells and %d timesteps\n", xl_to3, timesteps);
    printf("Total Execution Time: %lf seconds\n", ts_to_double(ts_diff(begin, end)));
    printf("MLUPS: %.3lf\n", mlups/mlups_time);
    printf("Freeing allocated memory...\n");
    free(collideField);
    free(streamField);
    free(flagField);
	return 0;
}

#endif

