#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include "qpp.h"

using namespace std;
using namespace qpp;

cmat to_su2(cmat u)
{

  cplx det = u.determinant();
  cplx one(1,0);
  return (sqrt(one / det) * u);
}

double trace_dist(cmat u, cmat v)
{
  cmat u_v = u - v;
  cmat u_v_dagger = (u_v.conjugate()).transpose();
  cmat sqrt_uv = (u_v * u_v_dagger).array().sqrt().matrix();
  
  // return (0.5*trace(sqrt_uv)).real();
  double p, q, r2;
	p = norm(sqrt_uv(0,0)) + norm(sqrt_uv(1,0));
	q = norm(sqrt_uv(0,1)) + norm(sqrt_uv(1,1));
	r2 = norm(sqrt_uv(0,0) * conj(sqrt_uv(0,1)) + conj(sqrt_uv(1,0))*sqrt_uv(1,1));
	return (p + q + sqrt(pow(p - q, 2) + 4 * r2))/2;

}

vector<string> binary_prod(int n)
{
  
  vector<string> bin_list;
  for (int i = 1; i <= n; i++)
  {
    
    for (int j = 0;j < pow(2,i); j++)
    {
      boost::dynamic_bitset<> B(i, j);
      string B_string;
      to_string(B,B_string);
      bin_list.push_back(B_string);
    }
    
  }
  
  return bin_list;

}

vector<cmat> create_unitaries(vector<cmat> base, int limit)
{

  vector<cmat> gate_list;
  vector<string> bin_list = binary_prod(limit);
  for (int i = 0; i < bin_list.size(); i++)
  {
    cmat u = gt.Id2;
    string bits = bin_list[i];
    for (int j = 0; j < bits.length(); j++)
    {
        char bit = bits[j];
        int index = int(bit) - 48;
        u = u * base[index];
    }
    gate_list.push_back(u);
  }
  return gate_list;
}

cmat find_closest_u(vector<cmat> gate_list, cmat u)
{
  double min_dist = 10;
  cmat min_u = gt.Id2;
  int min_gate = 0;
  for (int i = 0; i<gate_list.size(); i++)
  {
    double tr_dist = trace_dist(gate_list[i],u);
    if (tr_dist < min_dist)
    {
      min_dist = tr_dist;
      min_gate = i;
      min_u = gate_list[i];
    }
    
  }
  return min_u;
}

array<cplx,4> u_to_bloch(cmat u)
{
  array<cplx,4> axis;
  double angle = (acos((u(0,0) + u(1,1))/2.0)).real();
  double Sin = sin(angle);
  if (Sin < pow(10,-10))
  {
    cplx nx(0,0);
    cplx ny(0,0);
    cplx nz(1,0);
    cplx angle_2(2*angle,0);
    
    axis[0] = nx;
    axis[1] = ny;
    axis[2] = nz;
    axis[3] = angle_2;
  }
  else
  {
    cplx j_2(0,2); 
    cplx nx = (u(0,1) + u(1,0)) / (Sin * j_2);
    cplx ny = (u(0,1) - u(1,0)) / (Sin * 2);
    cplx nz = (u(0,0) - u(1,1)) / (Sin * j_2);
    cplx angle_2(2*angle,0);
    
    axis[0] = nx;
    axis[1] = ny;
    axis[2] = nz;
    axis[3] = angle_2;
    
    
  }
  
  return axis;
}


cmat diagonalize(cmat u) {
    assert(u.rows() == 2 && u.cols() == 2);

    Eigen::ComplexEigenSolver<cmat> solver(u);

    cmat eigenvectors = solver.eigenvectors();
    
    // cmat res = eigenvectors.inverse()*u*eigenvectors;
    // cout << "1";
    //cout << disp(eigenvectors)<< endl ;
    // cout << "2";
    // cout << disp(res) << endl;

    //return eigenvectors;
    std::complex<double> det = u(0,0) * u(1,1) - u(1,0) * u(1,1) ;
    std::complex<double> trace = u(0,0) + u(1,1) ; 
    std::complex<double> lambda1 = (trace + sqrt(trace * trace - 4.0 * det)) / 2.0;
    std::complex<double> lambda2 = (trace - sqrt(trace * trace - 4.0 * det)) / 2.0;

    std::complex<double> v1 = lambda1 - u(1,1);
    std::complex<double> v2 = lambda2 - u(1,1);

    cmat diagonal(2,2);
    diagonal << lambda1, 0. ,0. , lambda2;
    //cout << disp(diagonal)<< endl ;

    return eigenvectors;
}

vector<cmat> gc_decomp(cmat u) 
{
    array<cplx,4> u_bloch = u_to_bloch(u);
    std::complex<double> result(0, 1);
    std::complex<double> neg(-1, 0);

    // The angle phi calculation
    std::complex<double> phi = 2.0 * asin(sqrt(sqrt(0.5 - 0.5 * cos(u_bloch[3] / 2.0))));
    cmat v(2,2) ;
    cmat w(2,2) ;
    v << cos(phi / 2.0),neg*result*sin(phi / 2.0) , neg*result*sin(phi / 2.0), cos(phi / 2.0);
    //if(M_PI < phi)
    if(u_bloch[2].real() > 0)
    { 
        w << cos((2 * M_PI - phi) / 2.0), neg*sin( (2 * M_PI - phi) / 2.0) , sin( (2 * M_PI - phi) / 2.0), cos( (2 * M_PI - phi) / 2.0);
    }
    else{
        w << cos(phi / 2.0), neg*sin(phi / 2.0) , sin(phi / 2.0), cos(phi / 2.0);
    }
    cmat ud = (diagonalize(u));
    cmat vwvdwd = (diagonalize(v*w*v.adjoint()*w.adjoint()));
    cmat s = to_su2(ud*vwvdwd.adjoint());
    cmat v_hat = s*v*s.adjoint();
    cmat w_hat = s*w*s.adjoint();
    vector<cmat> vw = {v_hat,w_hat};

    return vw;
}

cmat sk_algo(vector<cmat> gate_list, cmat u, int n) 
{
    if(n == 0){
      cmat re = find_closest_u(gate_list, u);
      //cout << "distance of recurrsion" << n << " : " << trace_dist(re, u)  << '\n' ;      
      return re;
    }
    else{
      cmat u_p = sk_algo(gate_list, u , n-1);
      //cmat uup = normalize(u*(u_p.adjoint()));
      //cout << "uup" << disp(u*(u_p.adjoint())) << endl;
      vector<cmat> vw = gc_decomp(u*(u_p.adjoint()));
      cmat v = vw[0];
      cmat w = vw[1];
      // cout << "v" << disp(v) << endl ;
      // cout << "w" << disp(w) << endl ;
      cmat v_p = sk_algo(gate_list ,v , n-1);
      cmat w_p = sk_algo(gate_list ,w , n-1);
      cmat re = v_p*w_p*v_p.adjoint()*w_p.adjoint()*u_p;
      //cout << "distance of recurrsion" << n << " : " << trace_dist(re, u)  << '\n' ;      

      return re;}
}

vector<cmat> read_unitaries(const string& filename)
{
    vector<cmat> unitaries;
    ifstream file(filename, ios::binary);

    if (!file.is_open())
    {
        cerr << "Error opening file: " << filename << endl;
        return unitaries;
    }

    int num_matrices;
    file.read(reinterpret_cast<char*>(&num_matrices), sizeof(int));

    while (num_matrices-- > 0)
    {
        cmat u(2, 2); // Assuming 2x2 matrices, adjust accordingly
        file.read(reinterpret_cast<char*>(u.data()), sizeof(complex<double>) * u.size());
        unitaries.push_back(u);
    }

    file.close();
    return unitaries;
}


int main() {

    //cout << "Please cin 4 complex number\n";
    //double a,b;

    /*for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2;j++)
        {
          cin >> a >> b;
          cplx z(a,b);
          u(i,j) = z;
        
        }   
        
    }
    //Need to check whether the input is unitary
    */
    int lim;
    cout << "depth of gate list: " << '\n' ;
    cin >> lim;
    int limit;
    cout << "Number of recursion: " << '\n' ;
    cin >> limit;
    vector<cmat> gate_list = read_unitaries("gate_list_" + std::to_string(lim) + ".dat");
    cmat u(2,2);
    cmat v(2,2);
    u = randU();
    u = to_su2(u); 
    cout << " u is  " << '\n';
    cout << disp(u) << '\n';
    cout << "run sk algorithm ....." << '\n';
    cmat appro_x = sk_algo(gate_list, u ,limit);
    cout << "result : " << '\n';
    cout << disp(appro_x) << '\n' ;
    cout << "distance of matrix: " << trace_dist(appro_x, u)  << '\n' ;
    // cout << i << "th sk run : "<< '\n';
    // cout <<  u << '\n';
    // array<cplx,4> axis = u_to_bloch(u);
    // cout << "Bloch sphere representation of u is ";
    // for (int i = 0; i < 4; i++)
    // {
    //   cout << axis[i] << ' ';
    // }

}