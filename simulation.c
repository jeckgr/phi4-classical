#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
int main() {
    int nx, nt;
    double dx, dt, mu, lambda;
    FILE *f = fopen("input.txt", "r");
    fscanf(f,"%d\n%d\n%lf\n%lf\n%lf\n%lf",&nx,&nt,&dx,&dt,&mu,&lambda);
    fclose(f);
    int *sup = (int*)malloc(nx*sizeof(int));
    int *sdn = (int*)malloc(nx*sizeof(int));
    double *phi = (double*)malloc(nx*sizeof(double));
    double *pi = (double*)malloc(nx*sizeof(double));
    double *psi = (double*)malloc(nx*sizeof(double));
    double *chi = (double*)malloc(nx*sizeof(double));
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
    for (int t = 0; t < nt; ++t) {
        char stra[20] = "output_";
        char strb[20];
        char strc[20] = ".txt";
        sprintf(strb,"%d",t);
        strcat(stra,strb);
        strcat(stra,strc);
        f = fopen(stra, "w");
        for (int x = 0; x < nx; ++x) {
            fprintf(f,"%f\t%f\n",x*dx,phi[x]);
            psi[x] = (phi[sup[x]]-2*phi[x]+phi[sdn[x]])/(dx*dx)-mu*phi[x]-lambda*phi[x]*phi[x]*phi[x];
        }
        for (int x = 0; x < nx; ++x) {
            phi[x] = phi[x] + dt * pi[x] + dt*dt/2 * psi[x];
        }
        for (int x = 0; x < nx; ++x) {
            chi[x] = (phi[sup[x]]-2*phi[x]+phi[sdn[x]])/(dx*dx)-mu*phi[x]-lambda*phi[x]*phi[x]*phi[x];
        }
        for (int x = 0; x < nx; ++x) {
            pi[x] = pi[x] + dt/2 * (psi[x]+chi[x]);
        }
        fclose(f);
    }
    free(psi);
    free(chi);
    free(pi);
    free(phi);
    free(sup);
    free(sdn);
}
