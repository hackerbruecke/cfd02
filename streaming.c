#include "streaming.h"
#include "LBDefinitions.h"

void doStreaming(double *collideField, double *streamField,int *flagField,int xlength){

	int dx, dy, dz;
	double fi;
	/*Setting distribution function for each moving direction/lattice velocity of every particle*/
	for (int z = 1; z < xlength + 1; ++z) {
		for (int y = 1; y < xlength + 1; ++y) {
			for (int x = 1; x < xlength + 1; ++x) {
				for (int i = 0; i < Q; ++i) {

					/*dx = c_i_x*dt, dt = 1*/
					dx = LATTICEVELOCITIES[i][0];
					dy = LATTICEVELOCITIES[i][1];
					dz = LATTICEVELOCITIES[i][2];

					/*New value for our distribution function (DF) of the index 'i'

					(We set it to DF(i) of the next particle, whose i-th lattice velocity
					points towards considered particle (x,y,z))

					Position of that next particle is given by (x-dx, y-dy, z-dz)*/

					fi = collideField[Q * ((z-dz)*(xlength+2)*(xlength+2) + (y-dy)*(xlength+2) + x-dx) + i];
					streamField[Q * (z*(xlength+2)*(xlength+2) + y*(xlength+2) + x) + i] = fi;
				}
			}
		}
	}

}

