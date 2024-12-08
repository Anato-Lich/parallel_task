#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <omp.h> 

#ifndef M_PI
#define M_PI 3.1415926
#endif

int get_index(int i, int j, int k, int Ny, int Nz) 
{
    return i * Ny * Nz + j * Nz + k;
}

int get_global_index(int i, int j, int k, int Ny, int Nz, int x_offset, int y_offset, int z_offset) 
{
    int global_i = i + x_offset;
    int global_j = j + y_offset;
    int global_k = k + z_offset;

    return global_i * Ny * Nz + global_j * Nz + global_k;
}

double analytical_func(double x, double y, double z, double t, double Lx, double Ly, double Lz) 
{
    double at = M_PI * sqrt(1.0 / (Lx * Lx) + 1.0 / (Ly * Ly) + 4.0 / (Lz * Lz));
    return sin(M_PI * x / Lx) * sin(M_PI * y / Ly) * sin(2 * M_PI * z / Lz) * cos(at * t + 2 * M_PI);
}

void initialize(double* u, int Nx, int Ny, int Nz, double hx, double hy, double hz, 
                double Lx, double Ly, double Lz, double t, int x_offset, int y_offset, int z_offset) 
{   
    #pragma omp parallel for collapse(3)  // Параллелизация трех вложенных циклов
    for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int k = 1; k < Nz - 1; k++) {
                double x = (i - 1 + x_offset) * hx;
                double y = (j - 1 + y_offset) * hy;
                double z = (k - 1 + z_offset) * hz;
                u[get_index(i, j, k, Ny, Nz)] = analytical_func(x, y, z, t, Lx, Ly, Lz);
            }
        }
    }
}

void apply_boundary_conditions(double* u, int Nx, int Ny, int Nz, int* neighbors) 
{
    if (neighbors[0] == MPI_PROC_NULL)
    {
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                u[get_index(1, j, k, Ny, Nz)] = 0.0;
            }
        }
    }

    if (neighbors[1] == MPI_PROC_NULL)
    {
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < Ny; j++) {
            for (int k = 0; k < Nz; k++) {
                u[get_index(Nx - 2, j, k, Ny, Nz)] = 0.0;
            }
        }
    }

    if (neighbors[2] == MPI_PROC_NULL)
    {
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < Nx; i++) {
            for (int k = 0; k < Nz; k++) {
                u[get_index(i, 1, k, Ny, Nz)] = 0.0;
            }
        }
    }
    
    if (neighbors[3] == MPI_PROC_NULL)
    {
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < Nx; i++) {
            for (int k = 0; k < Nz; k++) {
                u[get_index(i, Ny-2, k, Ny, Nz)] = 0.0;
            }
        }
    }
}


void exchange_boundaries_3d(double *u, int local_Nx, int local_Ny, int local_Nz, MPI_Comm comm, int* neighbors) 
{
    MPI_Request requests[12];
    MPI_Status statuses[12];
    int req_count = 0;

    // Временные буферы для обмена
    double *send_buf_XL = (double*)malloc(local_Ny * local_Nz * sizeof(double)); // Отправка левому соседу по X
    double *recv_buf_XL = (double*)malloc(local_Ny * local_Nz * sizeof(double)); // Получение от левого соседа по X
    double *send_buf_XR = (double*)malloc(local_Ny * local_Nz * sizeof(double)); // Отправка правому соседу по X
    double *recv_buf_XR = (double*)malloc(local_Ny * local_Nz * sizeof(double)); // Получение от правого соседа по X

    double *send_buf_YL = (double*)malloc(local_Nx * local_Nz * sizeof(double)); // Отправка левому соседу по Y
    double *recv_buf_YL = (double*)malloc(local_Nx * local_Nz * sizeof(double)); // Получение от левого соседа по Y
    double *send_buf_YR = (double*)malloc(local_Nx * local_Nz * sizeof(double)); // Отправка правому соседу по Y
    double *recv_buf_YR = (double*)malloc(local_Nx * local_Nz * sizeof(double)); // Получение от правого соседа по Y

    double *send_buf_ZL = (double*)malloc(local_Nx * local_Ny * sizeof(double)); // Отправка левому соседу по Z
    double *recv_buf_ZL = (double*)malloc(local_Nx * local_Ny * sizeof(double)); // Получение от левого соседа по Z
    double *send_buf_ZR = (double*)malloc(local_Nx * local_Ny * sizeof(double)); // Отправка правому соседу по Z
    double *recv_buf_ZR = (double*)malloc(local_Nx * local_Ny * sizeof(double)); // Получение от правого соседа по Z

    // Заполнение буферов для обмена
    if (neighbors[0] != MPI_PROC_NULL) { // Левый сосед по X
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < local_Ny; j++) {
            for (int k = 0; k < local_Nz; k++) {
                send_buf_XL[j * local_Nz + k] = u[get_index(1, j, k, local_Ny, local_Nz)];
            }
        }
        MPI_Irecv(recv_buf_XL, local_Ny * local_Nz, MPI_DOUBLE, neighbors[0], 0, comm, &requests[req_count++]);
        MPI_Isend(send_buf_XL, local_Ny * local_Nz, MPI_DOUBLE, neighbors[0], 1, comm, &requests[req_count++]);
    }

    if (neighbors[1] != MPI_PROC_NULL) { // Правый сосед по X
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < local_Ny; j++) {
            for (int k = 0; k < local_Nz; k++) {
                send_buf_XR[j * local_Nz + k] = u[get_index(local_Nx - 2, j, k, local_Ny, local_Nz)];
            }
        }
        MPI_Irecv(recv_buf_XR, local_Ny * local_Nz, MPI_DOUBLE, neighbors[1], 1, comm, &requests[req_count++]);
        MPI_Isend(send_buf_XR, local_Ny * local_Nz, MPI_DOUBLE, neighbors[1], 0, comm, &requests[req_count++]);
    }

    if (neighbors[2] != MPI_PROC_NULL) { // Левый сосед по Y
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < local_Nx; i++) {
            for (int k = 0; k < local_Nz; k++) {
                send_buf_YL[i * local_Nz + k] = u[get_index(i, 1, k, local_Ny, local_Nz)];
            }
        }
        MPI_Irecv(recv_buf_YL, local_Nx * local_Nz, MPI_DOUBLE, neighbors[2], 2, comm, &requests[req_count++]);
        MPI_Isend(send_buf_YL, local_Nx * local_Nz, MPI_DOUBLE, neighbors[2], 3, comm, &requests[req_count++]);
    }

    if (neighbors[3] != MPI_PROC_NULL) { // Правый сосед по Y
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < local_Nx; i++) {
            for (int k = 0; k < local_Nz; k++) {
                send_buf_YR[i * local_Nz + k] = u[get_index(i, local_Ny - 2, k, local_Ny, local_Nz)];
            }
        }
        MPI_Irecv(recv_buf_YR, local_Nx * local_Nz, MPI_DOUBLE, neighbors[3], 3, comm, &requests[req_count++]);
        MPI_Isend(send_buf_YR, local_Nx * local_Nz, MPI_DOUBLE, neighbors[3], 2, comm, &requests[req_count++]);
    }

    if (neighbors[4] != MPI_PROC_NULL) { // Левый сосед по Z
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < local_Nx; i++) {
            for (int j = 0; j < local_Ny; j++) {
                send_buf_ZL[i * local_Ny + j] = u[get_index(i, j, 1, local_Ny, local_Nz)];
            }
        }
        MPI_Irecv(recv_buf_ZL, local_Nx * local_Ny, MPI_DOUBLE, neighbors[4], 4, comm, &requests[req_count++]);
        MPI_Isend(send_buf_ZL, local_Nx * local_Ny, MPI_DOUBLE, neighbors[4], 5, comm, &requests[req_count++]);
    }

    if (neighbors[5] != MPI_PROC_NULL) { // Правый сосед по Z
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < local_Nx; i++) {
            for (int j = 0; j < local_Ny; j++) {
                send_buf_ZR[i * local_Ny + j] = u[get_index(i, j, local_Nz - 2, local_Ny, local_Nz)];
            }
        }
        MPI_Irecv(recv_buf_ZR, local_Nx * local_Ny, MPI_DOUBLE, neighbors[5], 5, comm, &requests[req_count++]);
        MPI_Isend(send_buf_ZR, local_Nx * local_Ny, MPI_DOUBLE, neighbors[5], 4, comm, &requests[req_count++]);
    }

    // Ожидание завершения всех коммуникаций
    MPI_Waitall(req_count, requests, statuses);
    
    // Запись данных из recv_buf обратно в u
    if (neighbors[0] != MPI_PROC_NULL) { // Левый сосед по X
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < local_Ny; ++j) {
            for (int k = 0; k < local_Nz; ++k) {
                u[get_index(0, j, k, local_Ny, local_Nz)] = recv_buf_XL[j * local_Nz + k];
            }
        }
    }
    if (neighbors[1] != MPI_PROC_NULL) { // Правый сосед по X
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < local_Ny; ++j) {
            for (int k = 0; k < local_Nz; ++k) {
                u[get_index(local_Nx - 1, j, k, local_Ny, local_Nz)] = recv_buf_XR[j * local_Nz + k];
            }
        }
    }

    if (neighbors[2] != MPI_PROC_NULL) { // Левый сосед по Y
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < local_Nx; ++i) {
            for (int k = 0; k < local_Nz; ++k) {
                u[get_index(i, 0, k, local_Ny, local_Nz)] = recv_buf_YL[i * local_Nz + k];
            }
        }
    }
    if (neighbors[3] != MPI_PROC_NULL) { // Правый сосед по Y
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < local_Nx; ++i) {
            for (int k = 0; k < local_Nz; ++k) {
                u[get_index(i, local_Ny - 1, k, local_Ny, local_Nz)] = recv_buf_YR[i * local_Nz + k];
            }
        }
    }

    if (neighbors[4] != MPI_PROC_NULL) { // Левый сосед по Z
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < local_Nx; ++i) {
            for (int j = 0; j < local_Ny; ++j) {
                u[get_index(i, j, 0, local_Ny, local_Nz)] = recv_buf_ZL[i * local_Ny + j];
            }
        }
    } else { // Периодическая граница для Z (при отсутствии левого соседа)
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < local_Nx; ++i) {
            for (int j = 0; j < local_Ny; ++j) {
                u[get_index(i, j, 0, local_Ny, local_Nz)] = u[get_index(i, j, local_Nz - 2, local_Ny, local_Nz)];
            }
        }
    }

    if (neighbors[5] != MPI_PROC_NULL) { // Правый сосед по Z
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < local_Nx; ++i) {
            for (int j = 0; j < local_Ny; ++j) {
                u[get_index(i, j, local_Nz - 1, local_Ny, local_Nz)] = recv_buf_ZR[i * local_Ny + j];
            }
        }
    } else { // Периодическая граница для Z (при отсутствии правого соседа)
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < local_Nx; ++i) {
            for (int j = 0; j < local_Ny; ++j) {
                u[get_index(i, j, local_Nz - 1, local_Ny, local_Nz)] = u[get_index(i, j, 1, local_Ny, local_Nz)];
            }
        }
    }


    // Освобождение памяти временных буферов
    free(recv_buf_XL);
    free(recv_buf_XR);
    free(recv_buf_YL);
    free(recv_buf_YR);
    free(recv_buf_ZL);
    free(recv_buf_ZR);
}



 
void compute_function(int Nx, int Ny, int Nz, double* u_next, double* u_curr, double* u_prev,
                      double tau, double hx, double hy, double hz, MPI_Comm cart_comm, int* neigbors) {
    exchange_boundaries_3d(u_curr, Nx, Ny, Nz, cart_comm, neigbors);
    
    #pragma omp parallel for collapse(3)
    for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++) {
            for (int k = 1; k < Nz - 1; k++) {
                u_next[get_index(i, j, k, Ny, Nz)] =
                    2 * u_curr[get_index(i, j, k, Ny, Nz)] - u_prev[get_index(i, j, k, Ny, Nz)] +
                    tau * tau * (
                        (u_curr[get_index(i + 1, j, k, Ny, Nz)] - 2 * u_curr[get_index(i, j, k, Ny, Nz)] + u_curr[get_index(i - 1, j, k, Ny, Nz)]) / (hx * hx) +
                        (u_curr[get_index(i, j + 1, k, Ny, Nz)] - 2 * u_curr[get_index(i, j, k, Ny, Nz)] + u_curr[get_index(i, j - 1, k, Ny, Nz)]) / (hy * hy) +
                        (u_curr[get_index(i, j, k + 1, Ny, Nz)] - 2 * u_curr[get_index(i, j, k, Ny, Nz)] + u_curr[get_index(i, j, k - 1, Ny, Nz)]) / (hz * hz)
                    );
            }
        }
    }
    
    // exchange_boundaries_3d(u_next, Nx, Ny, Nz, rank, size, cart_comm, coords, dims);
}

void compute_errors(double* u_curr, int local_Nx, int local_Ny, int local_Nz,
                    double hx, double hy, double hz, double Lx, double Ly, double Lz,
                    double tau, int n,
                    int x_offset, int y_offset, int z_offset,
                    MPI_Comm comm, int rank, int Nx, int Ny, int Nz,
                    double* global_max_error, double* global_avg_error, int* max_error_coords) {
    double local_max_error = 0.0, local_avg_error = 0.0;
    int local_max_coords[3] = {0, 0, 0};
    
    // Подсчет локальных ошибок
    for (int i = 1; i < local_Nx - 1; i++) {
        for (int j = 1; j < local_Ny - 1; j++) {
            for (int k = 1; k < local_Nz - 1; k++) {
                double x = (i - 1 + x_offset) * hx;
                double y = (j - 1 + y_offset) * hy;
                double z = (k - 1 + z_offset) * hz;
                double analytical_value = analytical_func(x, y, z, n * tau, Lx, Ly, Lz);
                double error = fabs(u_curr[get_index(i, j, k, local_Ny, local_Nz)] - analytical_value);
                if (error > local_max_error) {
                    local_max_error = error;
                    local_max_coords[0] = i + x_offset;
                    local_max_coords[1] = j + y_offset;
                    local_max_coords[2] = k + z_offset;
                }
                local_avg_error += error;
            }
        }
    }

    // Сбор глобальных
    struct {
        double error;
        int rank;
    } local_error_rank = {local_max_error, rank}, global_error_rank;

    MPI_Allreduce(&local_error_rank, &global_error_rank, 1, MPI_DOUBLE_INT, MPI_MAXLOC, comm);

    int global_coords[3] = {0, 0, 0};
    if (rank == global_error_rank.rank) {
        for (int i = 0; i < 3; i++) {
            global_coords[i] = local_max_coords[i];
        }
    }
    MPI_Bcast(global_coords, 3, MPI_INT, global_error_rank.rank, comm);

    MPI_Allreduce(&local_avg_error, global_avg_error, 1, MPI_DOUBLE, MPI_SUM, comm);

    if (rank == 0) {
        *global_max_error = global_error_rank.error;
        *global_avg_error /= (Nx * Ny * Nz);
        for (int i = 0; i < 3; i++) {
            max_error_coords[i] = global_coords[i];
        }
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 5) {
        if (rank == 0) {
            fprintf(stderr, "Usage: %s Lx Ly Lz N\n", argv[0]);
        }
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    // Параметры задачи
    double Lx = atof(argv[1]);
    double Ly = atof(argv[2]);
    double Lz = atof(argv[3]);
    int N = atoi(argv[4]); // Глобальное количество узлов вдоль одной стороны

    int Nx = N, Ny = N, Nz = N; // Размеры сетки
    double hx = Lx / (Nx - 1);
    double hy = Ly / (Ny - 1);
    double hz = Lz / (Nz - 1);
    double tau = 0.0002; // Шаг по времени
    int time_steps = 20;

    // Настройка декартовой топологии MPI
    int dims[3] = {0, 0, 0};
    MPI_Dims_create(size, 3, dims); // Автоматическое определение размеров
    int periods[3] = {0, 0, 1};     // Периодическая граница только по Z
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods, 1, &cart_comm);

    int coords[3];
    MPI_Cart_coords(cart_comm, rank, 3, coords);

    int neighbors[6];
    MPI_Cart_shift(cart_comm, 0, 1, &neighbors[0], &neighbors[1]); // Соседи по X
    MPI_Cart_shift(cart_comm, 1, 1, &neighbors[2], &neighbors[3]); // Соседи по Y
    MPI_Cart_shift(cart_comm, 2, 1, &neighbors[4], &neighbors[5]); // Соседи по Z

    if (rank == 0) {
        printf("Cartesian topology initialized with dimensions: %d x %d x %d\n", dims[0], dims[1], dims[2]);
    }

    // Локальные размеры сетки
    int local_Nx = Nx / dims[0] + 2;
    int local_Ny = Ny / dims[1] + 2;
    int local_Nz = Nz / dims[2] + 2;

    // Если остаток от деления есть, последний процесс берет больше
    if (coords[0] == dims[0] - 1) 
    {
        local_Nx += Nx % dims[0];
    }
    if (coords[1] == dims[1] - 1) 
    {
        local_Ny += Ny % dims[1];
    }
    if (coords[2] == dims[2] - 1) 
    {
        local_Nz += Nz % dims[2];
    }

    printf("Rank %d: local_Nx=%d, local_Ny=%d, local_Nz=%d\n", rank, local_Nx, local_Ny, local_Nz);

    int x_offset = coords[0] * (Nx / dims[0]);
    int y_offset = coords[1] * (Ny / dims[1]);
    int z_offset = coords[2] * (Nz / dims[2]);

    printf("Rank %d, offset %d %d %d\n", rank, x_offset, y_offset, z_offset);

    // Выделяем память для локальных массивов
    double* u_curr = (double*)malloc(local_Nx * local_Ny * local_Nz * sizeof(double));
    double* u_prev = (double*)malloc(local_Nx * local_Ny * local_Nz * sizeof(double));
    double* u_next = (double*)malloc(local_Nx * local_Ny * local_Nz * sizeof(double));

    // Инициализация начальных условий
    initialize(u_curr, local_Nx, local_Ny, local_Nz, hx, hy, hz, Lx, Ly, Lz, 0.0, x_offset, y_offset, z_offset);
    initialize(u_prev, local_Nx, local_Ny, local_Nz, hx, hy, hz, Lx, Ly, Lz, tau, x_offset, y_offset, z_offset);

    apply_boundary_conditions(u_curr, local_Nx, local_Ny, local_Nz, neighbors);
    apply_boundary_conditions(u_prev, local_Nx, local_Ny, local_Nz, neighbors);

    double global_max_error = 0.0, global_avg_error = 0.0;
    int max_error_coords[3] = {0, 0, 0};


    double start_time = MPI_Wtime();
    for (int n = 1; n <= time_steps; n++) {
       // Вычисление следующего временного слоя
       compute_function(local_Nx, local_Ny, local_Nz, u_next, u_curr, u_prev, tau, hx, hy, hz, cart_comm, neighbors);
       // Подсчет ошибки
       double global_max_error, global_avg_error;
       int max_error_coords[3];
       compute_errors(u_curr, local_Nx, local_Ny, local_Nz, hx, hy, hz, Lx, Ly, Lz,
                      tau, n, x_offset, y_offset, z_offset, cart_comm, rank, Nx, Ny, Nz,
                      &global_max_error, &global_avg_error, max_error_coords);

       if (rank == 0) {
           printf("Time step %d: max_error = %e, avg_error = %e, max_error_coords = (%d, %d, %d)\n",
                   n, global_max_error, global_avg_error, 
                   max_error_coords[0], max_error_coords[1], max_error_coords[2]);
       }

        // Обновление массивов
        double* temp = u_prev;
        u_prev = u_curr;
        u_curr = u_next;
        u_next = temp;
        
        apply_boundary_conditions(u_curr, local_Nx, local_Ny, local_Nz, neighbors);
    }

    double end_time = MPI_Wtime();

    if (rank == 0) {
        printf("Total simulation time: %f seconds\n ", end_time - start_time);
    }
    
    int inner_Nx = local_Nx - 2; // Убираем граничные элементы
    int inner_Ny = local_Ny - 2;
    int inner_Nz = local_Nz - 2;

    double* u_inner = (double*)malloc(inner_Nx * inner_Ny * inner_Nz * sizeof(double));

    // Копируем данные из u_curr в u_inner, исключая границы
    #pragma omp parallel for collapse(3)
    for (int i = 1; i < local_Nx - 1; ++i) {
        for (int j = 1; j < local_Ny - 1; ++j) {
            for (int k = 1; k < local_Nz - 1; ++k) {
                int inner_index = get_index(i - 1, j - 1, k - 1, inner_Ny, inner_Nz);
                int curr_index = get_index(i, j, k, local_Ny, local_Nz);
                u_inner[inner_index] = u_curr[curr_index];
            }
        }
    }
    
    int inner_size = inner_Nx * inner_Ny * inner_Nz;
    double* global_u = NULL;

    if (rank == 0) {
        global_u = (double*)malloc(Nx * Ny * Nz * sizeof(double)); // Глобальный массив без перекрытий
    }

    MPI_Gather(u_inner, inner_size, MPI_DOUBLE, 
               global_u, inner_size, MPI_DOUBLE, 
               0, MPI_COMM_WORLD);

    if (rank == 0) {
        double* sorted_u = (double*)malloc(Nx * Ny * Nz * sizeof(double));

        for (int r = 0; r < size; ++r) {
            int coords_r[3];
            MPI_Cart_coords(cart_comm, r, 3, coords_r);

            int x_offset_r = coords_r[0] * (Nx / dims[0]);
            int y_offset_r = coords_r[1] * (Ny / dims[1]);
            int z_offset_r = coords_r[2] * (Nz / dims[2]);
            
            #pragma omp parallel for collapse(3)
            for (int i = 0; i < inner_Nx; ++i) {
                for (int j = 0; j < inner_Ny; ++j) {
                    for (int k = 0; k < inner_Nz; ++k) {
                        int global_index = get_index(x_offset_r + i, 
                                                     y_offset_r + j, 
                                                     z_offset_r + k, 
                                                     Ny, Nz);
                        int local_index = get_index(i, j, k, inner_Ny, inner_Nz);
                        sorted_u[global_index] = global_u[r * inner_size + local_index];
                    }
                }
            }
        }
        
        FILE* output_file = fopen("final_grid.bin", "wb");
        if (output_file) {
            fwrite(sorted_u, sizeof(double), Nx * Ny * Nz, output_file);
            fclose(output_file);
            printf("File saved\n");
        } else {
            printf("Error opening the file for saving data.\n");
        }
        free(sorted_u);
    } 
    
  // Очистка памяти
    free(u_curr);
    free(u_prev);
    free(u_next);
    free(global_u);
    MPI_Finalize();
    return EXIT_SUCCESS;
}