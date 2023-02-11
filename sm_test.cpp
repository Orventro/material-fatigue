#include <iostream>
#include <cmath>
#include <armadillo>

using namespace arma;

enum flatten_order{
    COLUMN, ROW
};

vec flatten(mat &m, flatten_order order){
    vec v(m.n_elem);
    if(order == ROW) {
        for (int i = 0; i < m.n_rows; i++)
            for (int j = 0; j < m.n_cols; j++)
                v(i*m.n_cols + j) = m(i,j);
    } else {
        for (int j = 0; j < m.n_cols; j++)
            for (int i = 0; i < m.n_rows; i++)
                v(j*m.n_rows + i) = m(i,j);
    }
    return v;
}

void next(mat &u, int n, mat &r, mat &th, double h, double tau, vec &phi, vec &psi) {
    mat b = reshape(u, 1, n*n);
    std::cout << b << std::endl;
    flatten(u, ROW);
}

int main() {
    int L = 5;
    sp_mat m(L, L);
    m(1, 2) = 5;
    vec a;
    mat u(L, L);
    next(u, L, u, u, 0, 0, a, a);
}
