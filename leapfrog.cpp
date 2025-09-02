#include <iostream>
#include <fstream>
#include <cmath>
int main() {
    int nx, nt;
    double dx, dt, m2, lambda;
    std::ifstream fin("input.txt");
    fin >> nx >> nt >> dx >> dt >> m2 >> lambda;
    fin.close();
    int *sup = new int[nx];
    int *sdn = new int[nx];
    double *phi = new double[nx];
    double *pi = new double[nx];
    double *psi = new double[nx];
    double *chi = new double[nx];
    dx = 2*M_PI*ceil((nx*dx)/(2*M_PI))/nx;
    for (int x = 0; x < nx; ++x) {
        if (x == 0) {
            sup[x] = x+1;
            sdn[x] = nx-1;
        }
        else if (x == nx-1) {
            sup[x] = 0;
            sdn[x] = x-1;
        }
        else {
            sup[x] = x+1;
            sdn[x] = x-1;
        }
        phi[x] = cos(x*dx);
        pi[x] = sin(x*dx);
    }
    std::ofstream fout;
    for (int t = 0; t < nt; ++t) {
        fout.open("output_"+std::to_string(t)+".txt");
        for (int x = 0; x < nx; ++x) {
            fout << x*dx << '\t' << phi[x] << '\n';
            psi[x] = (phi[sup[x]]-2*phi[x]+phi[sdn[x]])/(dx*dx)-m2*phi[x]-lambda*phi[x]*phi[x]*phi[x];
        }
        for (int x = 0; x < nx; ++x) {
            phi[x] = phi[x] + dt * pi[x] + dt*dt/2 * psi[x];
        }
        for (int x = 0; x < nx; ++x) {
            chi[x] = (phi[sup[x]]-2*phi[x]+phi[sdn[x]])/(dx*dx)-m2*phi[x]-lambda*phi[x]*phi[x]*phi[x];
        }
        for (int x = 0; x < nx; ++x) {
            pi[x] = pi[x] + dt/2 * (psi[x]+chi[x]);
        }
        fout.close();
    }
    delete[] psi;
    delete[] chi;
    delete[] pi;
    delete[] phi;
    delete[] sup;
    delete[] sdn;
}
