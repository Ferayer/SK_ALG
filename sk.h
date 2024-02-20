#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <bitset>
#include <Eigen/Dense>
#include <fstream>

using namespace std;
using namespace Eigen;


typedef Matrix< std::complex<double> , 2 , 2 > Matrix2cd;

struct u_mat{
    public:
        double a;
        double b;
        double c;
        double d;
        Matrix U;
        array<std::complex<double>, 4> axis;

        vector<int> arrpox_list; //store resulting construction
        
        u_mat(double a, double b, double c, double d, vector<int> _arrpox_list)
		: arrpox_list(_arrpox_list) 
    {
      U(0,0) = complex<double>(a, d);
      U(0,1) = complex<double>(c, b);
      U(1,0) = complex<double>(-c, a);
      U(1,1) = complex<double>(a, -d);
    }

    Matrix2cd random()
    {
        //generate random unitary gate
    }

    bool su2()
    {
        //cheeck whether the unitary is SU(2)
    }

    Matrix2cd to_su2() {
    std::complex<double> det = U.determinant();
    std::complex<double> one(1, 0);
    return (sqrt(one / det) * U);
    }

    double trace_dist(Matrix2cd v) {
    Matrix2cd u_v = U - v;
    Matrix2cd u_v_dagger = u_v.conjugate().transpose();
    Matrix2cd sqrt_uv = (u_v * u_v_dagger).array().sqrt().matrix();

    double p, q, r2;
	  p = norm(sqrt_uv(0,0)) + norm(sqrt_uv(1,0));
	  q = norm(sqrt_uv(0,1)) + norm(sqrt_uv(1,1));
	  r2 = norm(sqrt_uv(0,0) * conj(sqrt_uv(0,1)) + conj(sqrt_uv(1,0))*sqrt_uv(1,1));
	  return (p + q + sqrt(pow(p - q, 2) + 4 * r2))/2;
    }

    void u_to_bloch() {

    array<std::complex<double>, 4> _axis;
    double angle = (acos((U(0, 0) + U(1, 1)) / 2.0)).real();
    double Sin = sin(angle);
    if (Sin < pow(10, -10)) {
        std::complex<double> nx(0, 0);
        std::complex<double> ny(0, 0);
        std::complex<double> nz(1, 0);
        std::complex<double> angle_2(2 * angle, 0);

        _axis[0] = nx;
        _axis[1] = ny;
        _axis[2] = nz;
        _axis[3] = angle_2;
    } else {
        std::complex<double> j_2(0, 2);
        std::complex<double> nx = (u(0, 1) + u(1, 0)) / (Sin * j_2);
        std::complex<double> ny = (u(0, 1) - u(1, 0)) / (Sin * 2);
        std::complex<double> nz = (u(0, 0) - u(1, 1)) / (Sin * j_2);
        std::complex<double> angle_2(2 * angle, 0);

        _axis[0] = nx;
        _axis[1] = ny;
        _axis[2] = nz;
        _axis[3] = angle_2;
    }

    axis = _axis;}
    

    Matrix2cd diagonalize() {
    Eigen::ComplexEigenSolver<Matrix2cd> solver(U);
    Matrix2cd eigenvectors = solver.eigenvectors();
    
    return eigenvectors;
    }




};

class sk{
    
  public:
  sk(int _depth, int _recursion)
  : depth(_depth), : recursion(_recursion)
  {}

  Matrix2cd find_closest_u(vector<Matrix2cd> gate_list, Matrix2cd u); 

  vector<Matrix2cd> gc_decomp(Matrix2cd u);
  Matrix2cd sk_algo(vector<Matrix2cd> gate_list, Matrix2cd u, int n) ;

  vector<Matrix2cd> read_unitaries(const string &filename);
  //generate_gate_list
  //write gate list

  private:
    int depth;
    int recursion;

};