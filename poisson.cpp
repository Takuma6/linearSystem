/*
This is the sample code to solve Poisson eq. for the electric potential.
The way to set variables A, x, b and solve Ax=b with GMRES is already developed.
What I need is to bring the function setting phi here (and maybe particle class).
*/

#include <iostream>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include "roll.h"

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
private:
    double eps_f, eps_t, eps_b;
    double sig_f, sig_t, sig_b;
public:
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

/*
Functions for phi, epsilon
*/

void Make_phi_particle_sum(double *phi, double* phi_sum, Particle *p, const double radius){
  int *nlattice;
  nlattice = Ns;
  Make_phi_particle_sum_primitive(phi, phi_sum, p, DX, NP_domain, Sekibun_cell, nlattice, radius);
}

inline void Make_phi_particle_sum_primitive(double *phi,
                                            double *phi_sum,
                                            Particle *p,
                                            const double &dx,
                                            const int &np_domain,
                                            int **sekibun_cell,
                                            const int Nlattice[DIM],
                                            const double radius){
#pragma omp parallel for 
  for(int n = 0; n < Particle_Number; n++){
    double xp[DIM];
    for(int d = 0; d < DIM; d++) xp[d] = p[n].x[d];

    int x_int[DIM];
    double residue[DIM];
    int sw_in_cell = Particle_cell(xp, dx, x_int, residue);
    sw_in_cell = 1;
    
    int r_mesh[DIM];
    double dmy, dmy_phi;
    double r[DIM], x[DIM];
    for(int mesh = 0; mesh < np_domain; mesh++){
      Relative_coord(sekibun_cell[mesh], x_int, residue, sw_in_cell, Nlattice, dx, r_mesh, r);

      for(int d = 0; d < DIM; d++) x[d] = r_mesh[d]*DX;

      dmy     = Distance(x, xp);
      dmy_phi = Phi(dmy, radius);

#pragma omp atomic
      phi_sum[(r_mesh[0] * NY * NZ_) + (r_mesh[1] * NZ_) + r_mesh[2]] += dmy_phi;
    }
  }

  {
#pragma omp parallel for
    for(int i = 0; i < NX; i++){
      int im;
      for(int j = 0; j < NY; j++){
        for(int k = 0; k < NZ; k++){
          im = (i * NY * NZ_) + (j * NZ_) + k;
          phi[im] = MIN(phi_sum[im], 1.0);
        }
      }
    }
  }
}

void Make_epsilon(electric_coefficient *ep, double *eps, double *phi_sum){
  #pragma omp parallel for
    for(int i = 0; i < NX; i++){
      int im;
      for(int j = 0; j < NY; j++){
        im = i * NY + j;
        eps[im] = ep.eps_f + (ep.eps_p - ep.eps_f)*phi_sum[im]
      }
    }
}

/*
Functions for Poisson Eq.
*/
inline void insertCoefficient(int id, int i, int j, double w, std::vector<T>& coeffs, Eigen::VectorXd& b, const int n){
  int id1 = i*n+j;
        if(i==-1 || i==n){
        	i   = (i+n)%n;
        	id1 =  i*n+j ;
        }
  else  if(j==-1 || j==n){
    			j   = (j+n)%j;
    			id1 =  i*n+j ;
    		}
  coeffs.push_back(T(id,id1,w));              // unknown coefficient
}

void buildProblem(std::vector<T>& coefficients, Eigen::VectorXd& b, int n, int m){
  //b.setZero();
  b = Eigen::VectorXd::Random(m);
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      int id = i*n+j;
      insertCoefficient(id, i-1,j,   -1, coefficients, b, n);
      insertCoefficient(id, i+1,j,   -1, coefficients, b, n);
      insertCoefficient(id, i  ,j-1, -1, coefficients, b, n);
      insertCoefficient(id, i  ,j+1, -1, coefficients, b, n);
      insertCoefficient(id, i  ,j  ,  4, coefficients, b, n);
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
  input(5, 5);
	n = 1 << 5;  // size of the grid
  m = n*n   ;  // number of unknows (=number of pixels)
  electric_coefficient *ep = new electric_coefficient;
  ep->print_value();
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
	buildProblem(coefficients, b, n, m);
	SpMat A(m,m);
	A.setFromTriplets(coefficients.begin(), coefficients.end());
	

	// fill A and b
	Eigen::GMRES<SpMat> gmres(A);
	Eigen::VectorXd x;
  gmres.setMaxIterations(3000);
  gmres.setTolerance(1e-5);
  gmres.solveWithGuess(b, Eigen::VectorXd::Zero(m));
  x = gmres.solve(b);
	std::cout << "#iterations:     " << gmres.iterations() << std::endl;
	std::cout << "estimated error: " << gmres.error()      << std::endl;
	// update b, and solve again
	// x = solver.solve(b);

	return 0;
}