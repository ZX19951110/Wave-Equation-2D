#include <stdio.h>
#include "unit_test.h"

int sz;
int arr_sz;
int main(int argc, char *argv[]) {
	
	sz = (argc==2?atoi(argv[1]):GRID_SZ);
	sz = sz>=100?sz:GRID_SZ;
	arr_sz = sz * sz;
	
	double *rightdata= process_withoutomp();
	double *testdatausingtask= process_withompusingtask();
	double *testdatausingfor= process_withompusingfor();
	
	for(int i = 0; i < sz; i++) {
		for(int j = 0; j < sz; j++) {
			if (rightdata[i*sz+j]!=testdatausingtask[i*sz+j]) {
				printf("wrong task");
				return 0;
			}
		}
	}
	for(int i = 0; i < sz; i++) {
		for(int j = 0; j < sz; j++) {
			if (rightdata[i*sz+j]!=testdatausingfor[i*sz+j]) {
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
		data = (double*)malloc(sizeof(double)*arr_sz);
		olddata = (double*)malloc(sizeof(double)*arr_sz);
		newdata = (double*)malloc(sizeof(double)*arr_sz);
		
		#pragma omp parallel
		{
			#pragma omp single
			{
				#pragma omp task
				for(i = 0; i < arr_sz; i++) {
					olddata[i] = data[i] = 1.0;
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
						olddata[(i+20)*sz+j+20] = data[(i+20)*sz+j+20] += h * exp( -5 * (pow(x[i][j], 2 ) + pow(x[j][i], 2 )));
					}
				}
				#pragma omp taskwait
				
			}
		}
		for(i = 0; i < TIMES; i++) {
			int n, m;
			double pot;
			#pragma omp parallel for private(n,m,pot) schedule(auto)
			for( n = 0; n < sz; n++) {
				for( m = 0; m < sz; m++) {
					pot = data[(n+1 >= sz ? n : n+1)*sz+m]+
						data[(n-1 < 0 ? 0 : n-1)*sz+m]+
						data[(m+1 >= sz ? m : m+1)+n*sz]+
						data[(m-1 < 0 ? 0 : m-1)+n*sz]
						-4*data[n*sz+m];
					newdata[n * sz + m] = 
									(pow(C * dt, 2) * pot * 2 + 4 * data[n * sz + m] - olddata[n * sz + m] 
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
		printf("with omp task spend: %f\n",end-start);
		return data;
}
double *process_withompusingfor() {
		double start = omp_get_wtime();
		
		int i, j;
		double dt = 0.04, C = 16, K = 0.1, h = 6;
		double *data, *olddata, *newdata, *tmp;
		double x[PEAK_SZ][PEAK_SZ], linspace[PEAK_SZ], delta = 2.0/(PEAK_SZ-1.0);
		data = (double*)malloc(sizeof(double)*arr_sz);
		olddata = (double*)malloc(sizeof(double)*arr_sz);
		newdata = (double*)malloc(sizeof(double)*arr_sz);
		
		#pragma omp parallel
		{
			#pragma omp for private(i) schedule(auto)
			for(i = 0; i < arr_sz; i++) {
				olddata[i] = data[i] = 1.0;
			}
			
			
			#pragma omp for private(i,j) schedule(auto)
			for(i = 0; i < PEAK_SZ; i++){
				linspace[i] = -1.0 + delta * i;
				for(j = 0; j < PEAK_SZ; j++) {
					x[i][j] = linspace[i];
				}
			}
			
			#pragma omp barrier
			
			#pragma omp for private(i,j) schedule(auto)
			for(i = 0; i < PEAK_SZ; i++){
				for(j = 0; j < PEAK_SZ; j++){
					olddata[(i+20)*sz+j+20] = data[(i+20)*sz+j+20] += h * exp( -5 * (pow(x[i][j], 2 ) + pow(x[j][i], 2 )));
				}
			}
		}
		for(i = 0; i < TIMES; i++) {
			int n, m;
			double pot;
			#pragma omp parallel for private(n,m,pot) schedule(auto)
			for( n = 0; n < sz; n++) {
				for( m = 0; m < sz; m++) {
					pot = data[(n+1 >= sz ? n : n+1)*sz+m]+
													data[(n-1 < 0 ? 0 : n-1)*sz+m]+
													data[(m+1 >= sz ? m : m+1)+n*sz]+
													data[(m-1 < 0 ? 0 : m-1)+n*sz]
													-4*data[n*sz+m];
					newdata[n * sz + m] = 
									(pow(C * dt, 2) * pot * 2 + 4 * data[n * sz + m] - olddata[n * sz + m] 
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
		printf("with omp for spend: %f\n",end-start);
		return data;
}

// WITHOUTOMP PART
double *process_withoutomp() {
		double start = omp_get_wtime();
		int i, j;
		double dt = 0.04, C = 16, K = 0.1, h = 6;
		double *data, *olddata, *newdata, *tmp;
		double x[PEAK_SZ][PEAK_SZ], linspace[PEAK_SZ], delta = 2.0/(PEAK_SZ-1.0);
		data = (double*)malloc(sizeof(double)*arr_sz);
		olddata = (double*)malloc(sizeof(double)*arr_sz);
		newdata = (double*)malloc(sizeof(double)*arr_sz);
	
		for(i = 0; i < arr_sz; i++){
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
				data[(i+20)*sz+j+20] += h * exp( -5 * (pow(x[i][j], 2 ) + pow(x[j][i], 2 )));
			}
		}
		
		for(i = 0; i < arr_sz; i++){
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
		for( i = 0; i < sz; i++){
				for( j = 0; j < sz; j++){
						add_i = i+1 >= sz ? i : i+1;
						add_j = j+1 >= sz ? j : j+1;
						sub_i = i-1 < 0 ? 0 : i-1;
						sub_j = j-1 < 0 ? 0 : j-1;
						pot = data[add_i*sz+j]+
									data[sub_i*sz+j]+
									data[add_j+i*sz]+
									data[sub_j+i*sz]-
									4*data[i*sz+j];
						newdata[i * sz + j] = 
								( pow(C * dt, 2) * pot * 2 + 4 * data[i * sz + j] - olddata[i * sz + j] *(2 - K * dt) ) / (2 + K * dt);
				}
		}
}
// WITHOUTOMP PART END
