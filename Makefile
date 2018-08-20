CXX = g++
EIGEN_DIR = /usr/local/include/eigen3

all:
	$(CXX) -I $(EIGEN_DIR) sample_poisson.cpp -o poisson.out

clean:
	rm -f *~ *.out