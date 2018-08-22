CXX = g++
EIGEN_DIR = /usr/local/include/eigen3
EIGEN_LINKS = -I $(EIGEN_DIR) 
HDF_LINKS = -I/usr/local/include -L/usr/local/lib -lhdf5_hl -lhdf5

all:
	$(CXX) -std=c++11 poisson.cpp -o poisson.out $(EIGEN_LINKS) $(HDF_LINKS)

clean:
	rm -f *~ *.out