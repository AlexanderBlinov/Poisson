#include <stdlib.h>
#include <algorithm>
#include <math.h>

using namespace std;
// Neodnorodnost'
double F(double x, double y, double z) {
    return 3 * exp(x + y + z);
}

// Granichnoe uslovie pri x=0
double A0(double y, double z) {
    return exp(y + z);
}

// Granichnoe uslovie pri x=X
double A1(double y, double z, double X) {
    return exp(X + y + z);
}

// Granichnoe uslovie pri y=0
double B0(double x, double z) {
    return exp(x + z);
}

// Granichnoe uslovie pri y=Y
double B1(double x, double z, double Y) {
    return exp(x + Y + z);
}

// Granichnoe uslovie pri z=0
double C0(double x, double y) {
    return exp(x + y);
}

// Granichnoe uslovie pri z=Z
double C1(double x, double y, double Z) {
    return exp(x + y + Z);
}

int main() {
    const auto start = clock();

    double X = 1, Y = 1, Z = 1;
    double h1 = 0.01, h2 = 0.01, h3 = 0.01;

    int Nx = (int)(X / h1) - 1, Ny = (int)(Y / h2) - 1, Nz = (int)(Z / h3) - 1;

    int rit = 300;
    double w = 1.7;

    double *U = (double *)calloc((size_t) Nx * Ny * Nz, sizeof(double));

    for (int it = 0; it < rit; ++it) {
        for (int i = 0; i < Nx; ++i) {
            for (int j = 0; j < Ny; ++j) {
                for (int k = 0; k < Nz; ++k) {
                    double uim = i == 0 ? A0((j + 1) * h2, (k + 1) * h3) : U[((i - 1) * Ny + j) * Nz + k];
                    double uip = i == Nx - 1 ? A1((j + 1) * h2, (k + 1) * h3, X) : U[((i + 1) * Ny + j) * Nz + k];
                    double ujm = j == 0 ? B0((i + 1) * h1, (k + 1) * h3) : U[(i * Ny + j - 1) * Nz + k];
                    double ujp = j == Ny - 1 ? B1((i + 1) * h1, (k + 1) * h3, Y) : U[(i * Ny + j + 1) * Nz + k];
                    double ukm = k == 0 ? C0((i + 1) * h1, (j + 1) * h2) : U[(i * Ny + j) * Nz + k - 1];
                    double ukp = k == Nz - 1 ? C1((i + 1) * h1, (j + 1) * h2, Z) : U[(i * Ny + j) * Nz + k + 1];
                    double u = U[(i * Ny + j) * Nz + k];

                    U[(i * Ny + j) * Nz + k] = w * ((uip + uim) / (h1 * h1)
                                                    + (ujp + ujm) / (h2 * h2)
                                                    + (ukp + ukm) / (h3 * h3) - F((i + 1) * h1,
                                                                                  (j + 1) * h2,
                                                                                  (k + 1) * h3)) /
                                                   (2 / (h1 * h1) + 2 / (h2 * h2) + 2 / (h3 * h3)) + (1 - w) * u;
                }
            }
        }
    }

    const auto end = clock();
    printf("Time: %f\n", (double)(end - start) / CLOCKS_PER_SEC);

    FILE *f = fopen("output.txt", "w");
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            for (int k = 0; k < Nz; ++k) {
                fprintf(f, "%f ", U[(i * Ny + j) * Nz + k]);
            }
            fprintf(f, "\n");
        }
        fprintf(f, "-----------------------------------------------\n");
    }

    return 0;
}