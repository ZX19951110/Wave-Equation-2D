#include <time.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>

#define GRID_SZ 5000
#define ARR_SZ GRID_SZ * GRID_SZ
#define PEAK_SZ 31

double *process_withomp();
double *process_withomp1();
double *process_withoutomp();
void sequential_update_withomp(double *data, double *olddata, double *newdata, double, double, double );
void sequential_update_withomp1(double *data, double *olddata, double *newdata, double, double, double );
void sequential_update_withoutomp(double *data, double *olddata, double *newdata, double C, double K, double dt );
