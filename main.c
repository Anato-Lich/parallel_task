#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>  // Включаем OpenMP

#define IDX(i, j, k, Ny, Nz) ((i) * (Ny) * (Nz) + (j) * (Nz) + (k))

// Функция аналитического решения для варианта 2
double u_analytical(double x, double y, double z, double t, double Lx, double Ly, double Lz) {
    double at = M_PI * sqrt(1.0 / (Lx * Lx) + 1.0 / (Ly * Ly) + 4.0 / (Lz * Lz));
    return sin(M_PI * x / Lx) * sin(M_PI * y / Ly) * sin(2 * M_PI * z / Lz) * cos(at * t + 2 * M_PI);
}

// Инициализация начального состояния
void initialize(double *u, int Nx, int Ny, int Nz, double hx, double hy, double hz, double Lx, double Ly, double Lz) {
    #pragma omp parallel for collapse(3)  // Параллелизация трех вложенных циклов
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                double x = i * hx;
                double y = j * hy;
                double z = k * hz;
                u[IDX(i, j, k, Ny, Nz)] = u_analytical(x, y, z, 0, Lx, Ly, Lz); // начальные условия
            }
        }
    }
}

// Выделение памяти для 3D массива
double *allocate_3d_array(int Nx, int Ny, int Nz) {
    return (double *)malloc(Nx * Ny * Nz * sizeof(double));
}

// Освобождение памяти
void free_3d_array(double *array) {
    free(array);
}

// Основная функция
int main(int argc, char *argv[]) {
    if (argc < 5) {
        fprintf(stderr, "Usage: %s Lx Ly Lz N\n", argv[0]);
        return EXIT_FAILURE;
    }

    double Lx = atof(argv[1]);
    double Ly = atof(argv[2]);
    double Lz = atof(argv[3]);
    double N = atoi(argv[4]);
    int Nx = N, Ny = N, Nz = N; // размеры сетки
    double T = 1.0; // общее время моделирования
    double hx = Lx / (Nx - 1);
    double hy = Ly / (Ny - 1);
    double hz = Lz / (Nz - 1);
    double tau = 0.01; // временной шаг
    int time_steps = (int)(T / tau);

    // Выделение памяти для массивов
    double *u_curr = allocate_3d_array(Nx, Ny, Nz);
    double *u_prev = allocate_3d_array(Nx, Ny, Nz);
    double *u_next = allocate_3d_array(Nx, Ny, Nz);

    // Инициализация начальных условий
    initialize(u_curr, Nx, Ny, Nz, hx, hy, hz, Lx, Ly, Lz);
    initialize(u_prev, Nx, Ny, Nz, hx, hy, hz, Lx, Ly, Lz);

    // Основной цикл по времени
    for (int n = 1; n <= time_steps; n++) {
        #pragma omp parallel for collapse(3) // Параллелизация трех вложенных циклов
        for (int i = 1; i < Nx - 1; i++) {
            for (int j = 1; j < Ny - 1; j++) {
                for (int k = 1; k < Nz - 1; k++) {
                    // Разностная схема для внутренней части блока
                    u_next[IDX(i, j, k, Ny, Nz)] = 2 * u_curr[IDX(i, j, k, Ny, Nz)] - u_prev[IDX(i, j, k, Ny, Nz)] +
                                      tau * tau * (
                                        (u_curr[IDX(i + 1, j, k, Ny, Nz)] - 2 * u_curr[IDX(i, j, k, Ny, Nz)] + u_curr[IDX(i - 1, j, k, Ny, Nz)]) / (hx * hx) +
                                        (u_curr[IDX(i, j + 1, k, Ny, Nz)] - 2 * u_curr[IDX(i, j, k, Ny, Nz)] + u_curr[IDX(i, j - 1, k, Ny, Nz)]) / (hy * hy) +
                                        (u_curr[IDX(i, j, k + 1, Ny, Nz)] - 2 * u_curr[IDX(i, j, k, Ny, Nz)] + u_curr[IDX(i, j, k - 1, Ny, Nz)]) / (hz * hz)
                                      );
                }
            }
        }

        // Периодическое граничное условие для направления z
        for (int i = 1; i < Nx - 1; i++) {
            for (int j = 1; j < Ny - 1; j++) {
                // k = 0 и k = Nz-1 (периодическая граница)
                u_next[IDX(i, j, 0, Ny, Nz)] = u_next[IDX(i, j, Nz - 2, Ny, Nz)];
                u_next[IDX(i, j, Nz - 1, Ny, Nz)] = u_next[IDX(i, j, 1, Ny, Nz)];
            }
        }

        // Обновление массивов
        double *temp = u_prev;
        u_prev = u_curr;
        u_curr = u_next;
        u_next = temp;
    }

    // Освобождение памяти
    free_3d_array(u_curr);
    free_3d_array(u_prev);
    free_3d_array(u_next);

    return EXIT_SUCCESS;
}
