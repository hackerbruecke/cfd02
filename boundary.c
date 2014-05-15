#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

void treatBoundary(double *collideField, int* flagField,
		const double * const wallVelocity, int xlength) {
	int i, x, y, z, dx, dy, dz;
	int N = xlength; /*substitution for the sake of simplicity*/
	double *currentCell;
	double finv, cu = 0;
	double density;

	/*Going through the hole boundary domain*/

	/*For each boundary cell/particle we update all the distribution functions,
	that point towards fluid, according to no-slip and moving-wall boundary conditions*/

	/*Two opposite/parallel boundary planes*/
	for (int k = 0; k < 2; ++k) {
		x = k * (N + 1);
		for (y = 0; y < N + 2; ++y) {
			for (z = 0; z < N + 1; ++z) {
				for (i = 0; i < Q; ++i) {
					/*dx = c_i_x*dt, dt = 1*/
					dx = LATTICEVELOCITIES[i][0];
					dy = LATTICEVELOCITIES[i][1];
					dz = LATTICEVELOCITIES[i][2];

					/*Checking if the i-th directions is facing the fluid particles*/
					/*y-z plane*/
					if (x+dx > 0 && x+dx < N + 1 &&
						y+dy > 0 && y+dy < N + 1 &&
						z+dz > 0 && z+dz < N + 1) {

						/*i-th distribution fun of boundary cell becomes distribution fun of
						the inverse lattice velocity of the pointed particle/lattice */

						/*Inverse lattice velocity to c_i is accessed
						at Q-i-1-th position of (x+dx, y+dy, z+dz) cell*/

						finv = collideField[Q * ((z+dz)*(N+2)*(N+2) + (y+dy)*(N+2) + x+dx) + Q-i-1];
						collideField[Q * (z*(N+2)*(N+2) + y*(N+2) + x) + i] = finv;
					}
					/*x and y swapped, so that another pair of boundary planes may be accessed within the same iteration*/
					/*x-z plane*/
					if (y+dx > 0 && y+dx < N + 1 &&
						x+dy > 0 && x+dy < N + 1 &&
						z+dz > 0 && z+dz < N + 1) {

						finv = collideField[Q * ((z+dz)*(N+2)*(N+2) + (x+dy)*(N+2) + y+dx) + Q-i-1];
						collideField[Q * (z*(N+2)*(N+2) + x*(N+2) + y) + i] = finv;
					}
					/*x and z swapped, out of same reason*/
					/*x-y plane*/
					if (z+dx > 0 && z+dx < N + 1 &&
						y+dy > 0 && y+dy < N + 1 &&
						x+dz > 0 && x+dz < N + 1) {

						finv = collideField[Q * ((x+dz)*(N+2)*(N+2) + (y+dy)*(N+2) + z+dx) + Q-i-1];

						/*Additional term when we are considering moving-wall (top plane - k=1)*/
						/*Velocity of the moving-wall is taken into account*/
						if(k == 1) {
							currentCell = collideField + Q * ((x+dz)*(N+2)*(N+2) + (y+dy)*(N+2) + z+dx);
							computeDensity(currentCell, &density);
							cu = 0; /*dot product of c_i and velocity of wall*/
							for (int d = 0; d < D; ++d) {
								cu += LATTICEVELOCITIES[i][d] * wallVelocity[d];
							}
							finv += 2*LATTICEWEIGHTS[i]*density*cu/(C_S*C_S);
						}
						collideField[Q * (x*(N+2)*(N+2) + y*(N+2) + z) + i] = finv;
					}
				}
			}
		}
	}

}













