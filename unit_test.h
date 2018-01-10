#include <time.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>

#define GRID_SZ 3000
#define ARR_SZ GRID_SZ * GRID_SZ
#define PEAK_SZ 31
#define TIMES 30

double *process_withompusingtask();
double *process_withompusingfor();
double *process_withoutomp();

void sequential_update_withoutomp(double *data, double *olddata, double *newdata, double C, double K, double dt );
