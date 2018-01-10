#include <stdio.h>
#include "unit_test.h"

int main(int argc, char *argv[]) {
	
	double *rightdata= process_withoutomp();
	double *testdatausingtask= process_withompusingtask();
	double *testdatausingfor= process_withompusingfor();
	
	for(int i = 0; i < GRID_SZ; i++) {
		for(int j = 0; j < GRID_SZ; j++) {
			if (rightdata[i*GRID_SZ+j]!=testdatausingtask[i*GRID_SZ+j]) {
				printf("wrong task");
				return 0;
			}
		}
	}
	for(int i = 0; i < GRID_SZ; i++) {
		for(int j = 0; j < GRID_SZ; j++) {
			if (rightdata[i*GRID_SZ+j]!=testdatausingfor[i*GRID_SZ+j]) {
				printf("wrong for");
				return 0;
			}
		}
	}
	printf("no wrong\n");
}
double *process_withompusingtask() {
		double start = omp_get_wtime();
		
		int i, j;
		double dt = 0.04, C = 16, K = 0.1, h = 6;
		double *data, *olddata, *newdata, *tmp;
		double x[PEAK_SZ][PEAK_SZ], linspace[PEAK_SZ], delta = 2.0/(PEAK_SZ-1.0);
		data = (double*)malloc(sizeof(double)*ARR_SZ);
		olddata = (double*)malloc(sizeof(double)*ARR_SZ);
		newdata = (double*)malloc(sizeof(double)*ARR_SZ);
		
		#pragma omp parallel
		{
			#pragma omp single
			{
				#pragma omp task
				for(i = 0; i < ARR_SZ; i++) {
					data[i] = 1.0;
				}
				
				#pragma omp taskwait
				
				#pragma omp task
				for(i = 0; i < PEAK_SZ; i++){
					linspace[i] = -1.0 + delta * i;
					for(j = 0; j < PEAK_SZ; j++) {
						x[i][j] = linspace[i];
					}
				}
				
				#pragma omp taskwait
				
				#pragma omp task
				for(i = 0; i < PEAK_SZ; i++){
					for(j = 0; j < PEAK_SZ; j++){
						data[(i+20)*GRID_SZ+j+20] += h * exp( -5 * (pow(x[i][j], 2 ) + pow(x[j][i], 2 )));
					}
				}
				#pragma omp taskwait
				
				#pragma omp task
				for(i = 0; i < ARR_SZ; i++){
					olddata[i] = data[i];
				}
				
			}
		}
		for(i = 0; i < TIMES; i++) {
			int n, m;
			double pot;
			#pragma omp parallel for private(n,m,pot) schedule(auto)
			for( n = 0; n < GRID_SZ; n++) {
				for( m = 0; m < GRID_SZ; m++) {
					pot = data[(n+1 >= GRID_SZ ? n : n+1)*GRID_SZ+m]+
						data[(n-1 < 0 ? 0 : n-1)*GRID_SZ+m]+
						data[(m+1 >= GRID_SZ ? m : m+1)+n*GRID_SZ]+
						data[(m-1 < 0 ? 0 : m-1)+n*GRID_SZ]
						-4*data[n*GRID_SZ+m];
					newdata[n * GRID_SZ + m] = 
									(pow(C * dt, 2) * pot * 2 + 4 * data[n * GRID_SZ + m] - olddata[n * GRID_SZ + m] 
															* (2 - K * dt)) 
															/ (2 + K * dt);
										
				}
			}
			tmp = olddata;
			olddata = data;
			data = newdata;
			newdata = tmp;
		}
			
		
		double end = omp_get_wtime();
		printf("with omp spend: %f\n",end-start);
		return data;
}
double *process_withompusingfor() {
		double start = omp_get_wtime();
		
		int i, j;
		double dt = 0.04, C = 16, K = 0.1, h = 6;
		double *data, *olddata, *newdata, *tmp;
		double x[PEAK_SZ][PEAK_SZ], linspace[PEAK_SZ], delta = 2.0/(PEAK_SZ-1.0);
		data = (double*)malloc(sizeof(double)*ARR_SZ);
		olddata = (double*)malloc(sizeof(double)*ARR_SZ);
		newdata = (double*)malloc(sizeof(double)*ARR_SZ);
		
		#pragma omp parallel
		{
			#pragma omp parallel for private(i) schedule(auto)
			for(i = 0; i < ARR_SZ; i++) {
					data[i] = 1.0;
			}
			
			#pragma omp for private(i,j) schedule(auto)
			for(i = 0; i < PEAK_SZ; i++){
				linspace[i] = -1.0 + delta * i;
				for(j = 0; j < PEAK_SZ; j++) {
					x[i][j] = linspace[i];
				}
			}
			
			
			#pragma omp for private(i,j) schedule(auto)
			for(i = 0; i < PEAK_SZ; i++){
				for(j = 0; j < PEAK_SZ; j++){
					data[(i+20)*GRID_SZ+j+20] += h * exp( -5 * (pow(x[i][j], 2 ) + pow(x[j][i], 2 )));
				}
			}
			
			#pragma omp for private(i) schedule(auto)
			for(i = 0; i < ARR_SZ; i++){
				olddata[i] = data[i];
			}
		}
		for(i = 0; i < TIMES; i++) {
			int n, m;
			double pot;
			#pragma omp parallel for private(n,m,pot) schedule(auto)
			for( n = 0; n < GRID_SZ; n++) {
				for( m = 0; m < GRID_SZ; m++) {
					pot = data[(n+1 >= GRID_SZ ? n : n+1)*GRID_SZ+m]+
													data[(n-1 < 0 ? 0 : n-1)*GRID_SZ+m]+
													data[(m+1 >= GRID_SZ ? m : m+1)+n*GRID_SZ]+
													data[(m-1 < 0 ? 0 : m-1)+n*GRID_SZ]
													-4*data[n*GRID_SZ+m];
					newdata[n * GRID_SZ + m] = 
									(pow(C * dt, 2) * pot * 2 + 4 * data[n * GRID_SZ + m] - olddata[n * GRID_SZ + m] 
															* (2 - K * dt)) 
															/ (2 + K * dt);
										
				}
			}
				
			tmp = olddata;
			olddata = data;
			data = newdata;
			newdata = tmp;
		}
		
		double end = omp_get_wtime();
		printf("with omp spend: %f\n",end-start);
		return data;
}

// WITHOUTOMP PART
double *process_withoutomp() {
		double start = omp_get_wtime();
		int i, j;
		double dt = 0.04, C = 16, K = 0.1, h = 6;
		double *data, *olddata, *newdata, *tmp;
		double x[PEAK_SZ][PEAK_SZ], linspace[PEAK_SZ], delta = 2.0/(PEAK_SZ-1.0);
		data = (double*)malloc(sizeof(double)*ARR_SZ);
		olddata = (double*)malloc(sizeof(double)*ARR_SZ);
		newdata = (double*)malloc(sizeof(double)*ARR_SZ);
	
		for(i = 0; i < ARR_SZ; i++){
			data[i] = 1.0;
		}
	 
		for(i = 0; i < PEAK_SZ; i++){
			linspace[i] = -1.0 + delta * i;
		}
	 
		for(i = 0; i < PEAK_SZ; i++){
			for(j = 0; j < PEAK_SZ; j++){
				x[i][j] = linspace[i];
			}
		}
		
		for(i = 0; i < PEAK_SZ; i++){
			for(j = 0; j < PEAK_SZ; j++){
				data[(i+20)*GRID_SZ+j+20] += h * exp( -5 * (pow(x[i][j], 2 ) + pow(x[j][i], 2 )));
			}
		}
		
		for(i = 0; i < ARR_SZ; i++){
			olddata[i] = data[i];
		}

		for(i = 0; i < TIMES; i++){
				sequential_update_withoutomp( data, olddata, newdata, C, K, dt);
				tmp = olddata;
				olddata = data;
				data = newdata;
				newdata = tmp;
		}
		double end = omp_get_wtime();
		printf("without omp spend: %f\n",end-start);
		
		return data;
}
void sequential_update_withoutomp(double *data, double *olddata, double *newdata, double C, double K, double dt ){
		int i, j, add_i, sub_i, add_j, sub_j;
		double pot;
		for( i = 0; i < GRID_SZ; i++){
				for( j = 0; j < GRID_SZ; j++){
						add_i = i+1 >= GRID_SZ ? i : i+1;
						add_j = j+1 >= GRID_SZ ? j : j+1;
						sub_i = i-1 < 0 ? 0 : i-1;
						sub_j = j-1 < 0 ? 0 : j-1;
						pot = data[add_i*GRID_SZ+j]+
									data[sub_i*GRID_SZ+j]+
									data[add_j+i*GRID_SZ]+
									data[sub_j+i*GRID_SZ]-
									4*data[i*GRID_SZ+j];
						newdata[i * GRID_SZ + j] = 
								( pow(C * dt, 2) * pot * 2 + 4 * data[i * GRID_SZ + j] - olddata[i * GRID_SZ + j] *(2 - K * dt) ) / (2 + K * dt);
				}
		}
}
// WITHOUTOMP PART END
