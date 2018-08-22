/*
This is the sample code to solve Poisson eq. for the electric potential.
The way to set variables A, x, b and solve Ax=b with GMRES is already developed.
What I need is to bring the function setting phi here (and maybe particle class).
*/

#include <iostream>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
//#include "roll.h"
#include "SimIO.hpp"

typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
typedef Eigen::Triplet<double> T;

const int DIM = 2; 
int Ns[DIM];
int &Nx=Ns[0];
int &Ny=Ns[1];
int n, m;
double *eps, *phi, *phi_sum; 

class electric_coefficient
{
public:
    double eps_f, eps_t, eps_b;
    double sig_f, sig_t, sig_b;
    electric_coefficient();    //コンストラクタ
    void print_value();
};
electric_coefficient::electric_coefficient()    //コンストラクタの定義
{
  eps_f = 1 ;
  eps_t = 10;
  eps_b = 10;
  sig_f = 1 ;
  sig_t = 1 ;
  sig_b = 1 ;
}
void electric_coefficient::print_value(){std::cout << eps_f << std::endl;}

void input(int n1, int n2){
  Ns[0] = 1<<n1;
  Ns[1] = 1<<n2;
}

void Make_epsilon(electric_coefficient *ep, double *eps, double *phi_sum){
  double eps_p = ep->eps_t;
  #pragma omp parallel for
    for(int i = 0; i < Nx; i++){
      int im;
      for(int j = 0; j < Ny; j++){
        im = i * Ny + j;
        eps[im] = ep->eps_f + (eps_p - ep->eps_f)*phi_sum[im];
      }
    }
}

/*
Functions for Poisson Eq.
*/
inline void insertCoefficient(int id, int i, int j, double w, std::vector<T>& coeffs, Eigen::VectorXd& b){
  int id1 = i*n+j;
        if(i==-1 || i==Nx){
        	i   = (i+Nx)%Nx;
        	id1 =  i*Ny+j ;
        }
  else  if(j==-1 || j==Ny){
    			j   = (j+Ny)%j;
    			id1 =  i*Ny+j ;
    		}
  coeffs.push_back(T(id,id1,w));              // unknown coefficient
}

void buildProblem(std::vector<T>& coefficients, Eigen::VectorXd& b, int m){
  //b.setZero();
  b = Eigen::VectorXd::Random(m);
  for(int i=0; i<Nx; i++){
    for(int j=0; j<Ny; j++){
      int id = i*Ny+j;
      insertCoefficient(id, i-1,j,   -1, coefficients, b);
      insertCoefficient(id, i+1,j,   -1, coefficients, b);
      insertCoefficient(id, i  ,j-1, -1, coefficients, b);
      insertCoefficient(id, i  ,j+1, -1, coefficients, b);
      insertCoefficient(id, i  ,j  ,  4, coefficients, b);
    }
  }
}

/*
void saveAsBitmap(const Eigen::VectorXd& x, int n, const char* filename){
  Eigen::Array<unsigned char,Eigen::Dynamic,Eigen::Dynamic> bits = (x*255).cast<unsigned char>();
  QImage img(bits.data(), n,n,QImage::Format_Indexed8);
  img.setColorCount(256);
  for(int i=0;i<256;i++) img.setColor(i,qRgb(i,i,i));
  img.save(filename);
}
*/

int main(){
  input(4, 5);
  n = 1<<5   ;
  m = Nx*Ny  ;  // number of unknows (=number of pixels)
  //electric_coefficient *ep = new electric_coefficient;
  //ep->print_value();
  //(*ep).print_values();
  /*
  // set phi, epsilon
  Reset_phi(phi);
  Reset_phi(phi_sum);
  Make_phi_particle_sum(phi, phi_sum, p);
  Make_epsilon(ep, phi_sum);
  */

	// Assembly:
	std::vector<T> coefficients;            // list of non-zeros coefficients
	Eigen::VectorXd b(m);                   // the right hand side-vector resulting from the constraints
	buildProblem(coefficients, b, m);
	SpMat A(m,m);
	A.setFromTriplets(coefficients.begin(), coefficients.end());
	

	// solve
	Eigen::GMRES<SpMat> gmres(A);
	Eigen::VectorXd x;
  gmres.setMaxIterations(3000);
  gmres.setTolerance(1e-5);
  gmres.solveWithGuess(b, Eigen::VectorXd::Zero(m));
  x = gmres.solve(b);
	std::cout << "#iterations:     " << gmres.iterations() << std::endl;
	std::cout << "estimated error: " << gmres.error()      << std::endl;
	// update b, and solve again
	// x = solver.solve(b)

  // save to HDF5
  double *dmy = new double [m];
  for(int i=0; i<m; i++)dmy[i] = (double)x[i];

  SimIO::SimIO fp;
  fp.create_file("sample.h5");
  fp.create_group("params");
  fp.write_attr("dt", 1);
  SimIO::io_space vectorSpace = SimIO::Space::create(Nx, Ny);
  fp.write_data("velocity", vectorSpace, dmy);
  fp.close_file();
	return 0;
}