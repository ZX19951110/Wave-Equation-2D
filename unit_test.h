#include <time.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>

#define GRID_SZ 1000
#define ARR_SZ GRID_SZ * GRID_SZ
#define PEAK_SZ 31
#define TIMES 30

double *process_withompusingtask();
double *process_withompusingunrollingfor();
double *process_withompusingfor();

double *process_withoutomp();

void sequential_update_withompusingtask(double *data, double *olddata, double *newdata, double, double, double );
void sequential_update_withompusingunrollingfor(double *data, double *olddata, double *newdata, double, double, double );
void sequential_update_withompusingfor(double *data, double *olddata, double *newdata, double, double, double );

void sequential_update_withoutomp(double *data, double *olddata, double *newdata, double C, double K, double dt );
