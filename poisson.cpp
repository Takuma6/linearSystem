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
double Ext[DIM] = {1,0}; 
const double dx=0.5;
double *phi, *phi_sum; 
double *eps, *bc_eps ;

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
  eps_f = 2 ;
  eps_t = 80;
  eps_b = 80;
  sig_f = 1 ;
  sig_t = 1 ;
  sig_b = 1 ;
}
void electric_coefficient::print_value(){std::cout << eps_f << std::endl;}

void Make_epsilon(electric_coefficient *ep){
  double eps_p = ep->eps_t;
  //#pragma omp parallel for
    for(int i = 0; i < Nx; i++){
      int im;
      for(int j = 0; j < Ny; j++){
        im = i * Ny + j;
        eps[im] = ep->eps_f + (eps_p - ep->eps_f)*phi[im];
      }
    }
}

void Make_bc_eps(){
  int im, im_bc;
  int i , j    ;
  for(int i_bc=0; i_bc<Nx+2; i_bc++){
    for(int j_bc=0; j_bc<Ny+2; j_bc++){
      i = (i_bc+Nx-1)%Nx;
      j = (j_bc+Ny-1)%Ny;
      im    = i   * Ny   +j   ;
      im_bc = i_bc*(Ny+2)+j_bc;
      bc_eps[im_bc] = eps[im] ;
    }
  }
}

void initialize(int n1, int n2, int &m, electric_coefficient *ep){
  Ns[0] = 1<<n1;
  Ns[1] = 1<<n2;
  m     = Nx*Ny;
  phi   = new double [m];
  eps   = new double [m];
  bc_eps= new double [(Nx+2)*(Ny+2)];
  SimIO::SimIO fp_dmy;
  fp_dmy.open_file("phi.h5", "r");
  fp_dmy.read_data("phi", phi);
  fp_dmy.close_file();
  Make_epsilon(ep);
  Make_bc_eps();
}

/*
Functions for Poisson Eq.
*/
inline void insertCoefficient(int id, int i, int j, double w, std::vector<T>& coeffs){
  int i_star = (i+Nx)%Nx        ;
  int id1    = (i_star*Ny+j+m)%m;
  coeffs.push_back(T(id,id1,w/(2*dx*dx)));              // unknown coefficient
}

void buildProblem(std::vector<T>& coefficients, Eigen::VectorXd& b, int m){
  b.setZero();
  //b = Eigen::VectorXd::Random(m);
  for(int i=0; i<Nx; i++){
    for(int j=0; j<Ny; j++){
      int id    = i*Ny+j;
      int id_bc = (i+1)*(Ny+2)+(j+1);
      b[id] = (bc_eps[id_bc+(Ny+2)] - bc_eps[id_bc-(Ny+2)])*Ext[0]/(2*dx)
             +(bc_eps[id_bc+1     ] - bc_eps[id_bc-1     ])*Ext[1]/(2*dx);
      insertCoefficient(id, i-1,j,     bc_eps[id_bc       ] + bc_eps[id_bc-(Ny+2)], coefficients);
      insertCoefficient(id, i+1,j,     bc_eps[id_bc       ] + bc_eps[id_bc+(Ny+2)], coefficients);
      insertCoefficient(id, i  ,j-1,   bc_eps[id_bc       ] + bc_eps[id_bc-1     ], coefficients);
      insertCoefficient(id, i  ,j+1,   bc_eps[id_bc       ] + bc_eps[id_bc+1     ], coefficients);
      insertCoefficient(id, i  ,j  , -(bc_eps[id_bc-(Ny+2)] + bc_eps[id_bc+(Ny+2)]
                                     + bc_eps[id_bc-1     ] + bc_eps[id_bc+1     ]
                                     + bc_eps[id_bc       ] * 4)                  , coefficients);
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
  electric_coefficient *ep = new electric_coefficient;
  initialize(6, 6, m, ep);
  //ep->print_value();
  (*ep).print_value();
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
  //gmres.setTolerance(1e-10);
  //gmres.solveWithGuess(b, Eigen::VectorXd::Zero(m));
  x = gmres.solve(b);
	std::cout << "#iterations:     " << gmres.iterations() << std::endl;
	std::cout << "estimated error: " << gmres.error()      << std::endl;
	// update b, and solve again
	// x = solver.solve(b)
  double *dmy = new double [m];
  for(int i=0; i<m; i++)dmy[i]=x[i];

  // save to HDF5
  SimIO::SimIO fp;
  fp.create_file("sample.h5");
  fp.create_group("field");
  fp.write_attr("dt", 1);
  SimIO::io_space matrixSpace = SimIO::Space::create(Nx, Ny);
  fp.write_data("potential", matrixSpace, dmy);
  fp.close_file();
	return 0;
}