#include <sk.h>

  Matrix2cd sk::find_closest_u(vector<Matrix2cd> gate_list, Matrix2cd u) {
    double min_dist = 10;
    Matrix2cd min_u = Matrix2cd::Identity();
    int min_gate = 0;
    for (int i = 0; i < gate_list.size(); i++) {
        double tr_dist = trace_dist(gate_list[i], u);
        if (tr_dist < min_dist) {
            min_dist = tr_dist;
            min_gate = i;
            min_u = gate_list[i];
        }
    }
    return min_u;
  }

  vector<Matrix2cd> sk::gc_decomp(Matrix2cd u) {
    array<std::complex<double>, 4> u_bloch = u_to_bloch(u);
    std::complex<double> result(0, 1);
    std::complex<double> neg(-1, 0);

    std::complex<double> phi = 2.0 * asin(sqrt(sqrt(0.5 - 0.5 * cos(u_bloch[3] / 2.0))));
    Matrix2cd v, w;
    v << cos(phi / 2.0), neg * result * sin(phi / 2.0), neg * result * sin(phi / 2.0), cos(phi / 2.0);
    if (u_bloch[2].real() > 0) {
        w << cos((2 * M_PI - phi) / 2.0), neg * sin((2 * M_PI - phi) / 2.0), sin((2 * M_PI - phi) / 2.0),
            cos((2 * M_PI - phi) / 2.0);
    } else {
        w << cos(phi / 2.0), neg * sin(phi / 2.0), sin(phi / 2.0), cos(phi / 2.0);
    }
    Matrix2cd ud = diagonalize(u);
    Matrix2cd vwvdwd = diagonalize(v * w * v.adjoint() * w.adjoint());
    Matrix2cd s = to_su2(ud * vwvdwd.adjoint());
    Matrix2cd v_hat = s * v * s.adjoint();
    Matrix2cd w_hat = s * w * s.adjoint();
    vector<Matrix2cd> vw = {v_hat, w_hat};


    if(trace_dist(v_hat*w_hat*v_hat.adjoint()*w_hat.adjoint(),u) > 1e-10){
      cout << "dist error" << endl;
      cout << trace_dist(v_hat*w_hat*v_hat.adjoint()*w_hat.adjoint(),u) << endl ;
      cout << "vwvvwwvv" << endl << v_hat*w_hat*v_hat.adjoint()*w_hat.adjoint() << endl;
      cout << "u" << endl << u << endl;
      cout << "s" << endl << s << endl;
      cout << "=========================================" << endl;
    }
    // else {
    //   cout << u_bloch[0]  << "  " << u_bloch[1]  << "  " << u_bloch[2]  << "  " << u_bloch[3]  << "  " ;
    //   cout << "phi  " << phi << endl;
    //   cout << trace_dist(v_hat*w_hat*v_hat.adjoint()*w_hat.adjoint(),u) << endl ;
    //   cout << "ud" << endl << ud << endl;
    //   cout << "vwvdwd" << endl << vwvdwd << endl;
    //   cout << "s" << endl << s << endl;
    //   cout << "v_h" << endl << v_hat << endl;
    //   cout << "w_h" << endl << w_hat << endl;
    //   cout << "=========================================" << endl;
    // }
    return vw;
    }
    Matrix2cd sk::sk_algo(vector<Matrix2cd> gate_list, Matrix2cd u, int n) {
    if (n == 0) {
        Matrix2cd re = find_closest_u(gate_list, u);
        return re;
    } 
    else {
        Matrix2cd u_p = sk_algo(gate_list, u, n - 1);
        vector<Matrix2cd> vw = gc_decomp(u * (u_p.adjoint()));
        Matrix2cd v = vw[0];
        Matrix2cd w = vw[1];
        Matrix2cd v_p = sk_algo(gate_list, v, n - 1);
        Matrix2cd w_p = sk_algo(gate_list, w, n - 1);
        Matrix2cd re = v_p * w_p * v_p.adjoint() * w_p.adjoint() * u_p;
        return re;
    }
  }

  vector<Matrix2cd> sk::read_unitaries(const string &filename) {
    vector<Matrix2cd> unitaries;
    ifstream file(filename, ios::binary);

    if (!file.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return unitaries;
    }

    int num_matrices;
    file.read(reinterpret_cast<char *>(&num_matrices), sizeof(int));

    while (num_matrices-- > 0) {
        Matrix2cd u(2, 2); // Assuming 2x2 matrices, adjust accordingly
        file.read(reinterpret_cast<char *>(u.data()), sizeof(complex<double>) * u.size());
        unitaries.push_back(u);
    }

    file.close();
    return unitaries;
  }

