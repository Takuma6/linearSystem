#include <iostream>
#include <stdlib.h>

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