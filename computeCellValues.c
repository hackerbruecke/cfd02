#include <stdio.h>
#include "LBDefinitions.h"
#include "computeCellValues.h"
#include <math.h>
#include <string.h>

void computeDensity(const double *const currentCell, double *density){
    /* Sum of all fi in the Cell to compute density */
    *density = 0;
    for (int i=0; i<Q; ++i) {
        *density += currentCell[i];
    }
}

void computeVelocity(const double * const currentCell, const double * const density, double *velocity) {
    /* Initializing velocity vector to 0 */
    velocity[0] = 0.0;
    velocity[1] = 0.0;
    velocity[2] = 0.0;

    /* Compute velocity by momentum equation */
    for (int j = 0; j < D; ++j) {
        for (int i = 0; i < Q; ++i) {
           velocity[j] += currentCell[i] * LATTICEVELOCITIES[i][j];
        }
        velocity[j] /= *density;
    }    
}

void computeFeq(const double * const density, const double * const velocity, double *feq){
	for (int i=0; i<Q; ++i) {
        double c_dot_u = 0.0;
        double u_dot_u = 0.0;
		for (int j = 0; j < D; ++j) {
            /* Compute dot products needed for the computation of feq */
            c_dot_u += LATTICEVELOCITIES[i][j] * velocity[j];
            u_dot_u += velocity[j] * velocity[j];
		}
		feq[i] = LATTICEWEIGHTS[i] * (*density) * (1 + (c_dot_u)/(C_S*C_S) +
			(c_dot_u)*(c_dot_u)/(2*C_S*C_S*C_S*C_S) - u_dot_u/(2*C_S*C_S));
	}
}

