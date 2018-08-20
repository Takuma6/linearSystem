#include <iostream>
#include <stdlib.h>

double *a;
int N, Nx, Ny, im;

/*
Assume 2D Matrix A
index
   j 0 1 2  
	|0 1 2| 0 i 
A = |3 4 5| 1
	|6 7 8| 2
in the style of 1D vector, 
A = |0 1 2 3 4 5 6 7 8|
  => A[i,j] = i*Ny+j 

when shift in x in a backward style,
	|2 0 1| 
A = |5 3 4| 
	|8 6 7| 
A = |2 0 1 5 3 4 8 6 7|

in y,
	|6 7 8| 
A = |0 1 2| 
	|3 4 5| 
A = |6 7 8 0 1 2 3 4 5|

*/

void roll_fw(double *array){
	for(int i=0; i<N-1; i++){std::swap(array[i], array[i+1]);}
}
void roll_bw(double *array){
	for(int i=N-1; i>0; i--){std::swap(array[i-1], array[i]);}
}
void roll_fw_2dx(double *array){
	for(int i=0; i<Nx; i++){
		for(int j=0; j<Ny-1; j++){
			im = i*Ny+j;
			std::swap(array[im], array[im+1]);
		}
	}
}
void roll_bw_2dx(double *array){
	for(int i=0; i<Nx; i++){
		for(int j=Ny-1; j>0; j--){
			im = i*Ny+j;
			std::swap(array[im-1], array[im]);
		}
	}
}
void roll_fw_2dy(double *array){
	for(int j=0; j<Ny; j++){
		for(int i=0; i<Nx-1; i++){
			im = i*Ny+j;
			std::swap(array[im], array[im+Ny]);
		}
	}
}
void roll_bw_2dy(double *array){
	for(int j=0; j<Ny; j++){
		for(int i=Nx-1; i>0; i--){
			im = i*Ny+j;
			std::swap(array[im-Ny], array[im]);
		}
	}
}

void initialize(double *array){
	for(int i=0; i<N; i++){array[i]=(double)i;}
}
void initialize_2d(double *array){
	for(int i=0; i<Nx; i++){
		for(int j=0; j<Ny; j++){
			im = i*Ny+j;
			array[im]=(double)im;
		}
	}
}

void print_values(double *array){
	for(int i=0; i<N; i++){printf("array[%d] = %f\n", i, array[i]);}	
}
void print_values_2d(double *array){
	for(int i=0; i<Nx; i++){printf("%f %f %f\n", array[i*Ny], array[i*Ny+1], array[i*Ny+2]);fflush(stdout);}	
}

int main(){
	Nx=3;
	Ny=3;
	static const int MEMORY_ALIGNMENT = 128;
	double* p = (double*)malloc(sizeof(double) * Nx * Ny);
	initialize_2d(p);
	print_values_2d(p);
	roll_fw_2dx(p);
	print_values_2d(p);
	roll_bw_2dx(p);
	print_values_2d(p);
	roll_fw_2dy(p);
	print_values_2d(p);
	roll_bw_2dy(p);
	print_values_2d(p);
	return 0;
}