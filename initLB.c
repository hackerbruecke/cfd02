#include "initLB.h"
#include "LBDefinitions.h"

int readParameters(int *xlength, double *tau, double *velocityWall, int *timesteps, int *timestepsPerPlotting, int argc, char *argv[]){
	/* Checking if there are exactly 2 parameters in the input. Else give Error message. */
	if (argc != 2) {
		printf("You must provide a filename!\n");
		return -1;
	}
	const char *filename = argv[1];
	/* Reading all the input parameters */
	READ_INT(filename, *xlength);
	READ_DOUBLE(filename, *tau);
    read_double(filename, "vwallx", &velocityWall[0]);
    read_double(filename, "vwally", &velocityWall[1]);
    read_double(filename, "vwallz", &velocityWall[2]);
	READ_INT(filename, *timesteps);
	READ_INT(filename, *timestepsPerPlotting);

	return 0;
}


void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength){
    /* Set all values of flagField to 0, i.e. FLUID state. We will apply boundary conditions
        in the following */
    memset(flagField, FLUID, (xlength + 2)*(xlength + 2)*(xlength + 2) * sizeof(*flagField));
    /* The values for Boundary on Z = 0 set to No_Slip and Zmax plane set to Moving_Wall */ 

    for (int x = 0; x < xlength + 2; ++x) {
        for (int y = 0; y < xlength + 2; ++y) {
            flagField[fidx(xlength, x, y, 0)] = NO_SLIP;
            flagField[fidx(xlength, x, y, xlength+1)] = MOVING_WALL;
        }
    }

/* The values for Boundary on Y = 0 and Ymax plane set to No_Slip */
    for (int x = 0; x < xlength + 2; ++x) {
        for (int z = 0; z < xlength + 2; ++z) {
            flagField[fidx(xlength, x, 0, z)] = NO_SLIP;
            flagField[fidx(xlength, x, xlength+1, z)] = NO_SLIP;
        }
    }
/* The values for Boundary on Y = 0 and Ymax plane set to No_Slip */
    for (int y = 0; y < xlength + 2; ++y) {
        for (int z = 0; z < xlength + 2; ++z) {
            flagField[fidx(xlength, 0, y, z)] = NO_SLIP;
            flagField[fidx(xlength, xlength+1, y, z)] = NO_SLIP;
        }
    }
/* Stream and Collide Fields are initialized to the respective Latticeweights of the Cell */
    for (int x = 0; x < xlength + 2; ++x) {
        for (int y = 0; y < xlength + 2; ++y) {
            for (int z = 0; z < xlength + 2; ++z) {
                for (int i = 0; i < Q; ++i) {
                    collideField[idx(xlength, x, y, z, i)] = LATTICEWEIGHTS[i];
                    streamField[idx(xlength, x, y, z, i)] = LATTICEWEIGHTS[i];
                }
            }
        }
    }
}

