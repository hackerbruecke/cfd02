#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_

#define Q 19
#define D 3

static const int LATTICEVELOCITIES[19][3] =
{
    {  0, -1, -1 }, { -1,  0, -1 }, {  0,  0, -1 }, {  1,  0, -1 }, {  0,  1, -1 },
    { -1, -1,  0 }, {  0, -1,  0 }, {  1, -1,  0 }, { -1,  0,  0 }, {  0,  0,  0 },
    {  1,  0,  0 }, { -1,  1,  0 }, {  0,  1,  0 }, {  1,  1,  0 }, {  0, -1,  1 },
    { -1,  0,  1 }, {  0,  0,  1 }, {  1,  0,  1 }, {  0,  1,  1 }
};
static const double LATTICEWEIGHTS[19] = 
{
    1.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36,
    1.0/36, 2.0/36, 1.0/36, 2.0/36, 12.0/36,
    2.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36,
    1.0/36, 2.0/36, 1.0/36, 1.0/36
};
static const double C_S = 0.57735026919l; 

static const int FLUID = 0;
static const int NO_SLIP = 1;
static const int MOVING_WALL = 2;

static inline int idx(int xlength, int x, int y, int z, int i) {
    return Q * (z * (xlength+2) * (xlength+2) + y * (xlength+2) + x) + i;
}

static inline int fidx(int xlength, int x, int y, int z) {
    return z * (xlength+2) * (xlength+2) + y * (xlength+2) + x;
}

static inline int inv(int idx) {
    return 18 - idx;
}
#endif

