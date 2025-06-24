#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>
#include <math.h>

double P[2][2] = {{1.0,3.0},{-1.5, 1.5}};
double eps = 10.0;

int k_func(double x, double y){ //функция фиктивной области
    if ((x*x - 4*y*y > 1.0)&((x >= 1)&(x <= 3.0))){
        return 1;
    }
    return 0;
}

double intersection_x(double x, double y1, double y2){
    double y_inter = sqrt((pow(x, 2)-1)/4);
    if (y_inter > y2){
        return -y_inter;
    }
    return y_inter;
}

double intersection_y(double y){
    return sqrt(pow(y, 2)*4+1);
}

void print_2d(double** arr, int n, int m){
    printf("\n");
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            printf("%2.1lf ", arr[i][j]);
        }
        printf("\n");
    }
}

double operator(double** A, double** B, double** T, int i, int j, double h1, double h2){ //разностная схема
    double first = (A[i+1][j]*(T[i+1][j]-T[i][j])-A[i][j]*(T[i][j]-T[i-1][j]))*(-1);
    double second = (B[i][j+1]*(T[i][j+1]-T[i][j])-B[i][j]*(T[i][j]-T[i][j-1]))/(-1);
    return first + second;
}

double F_value(double x, double y, double h1, double h2){//вычисление правой части
    double sum = 0;
    for (int i = -1; i < 2; i+=2){
        for (int j = -1; j < 2; j+=2){
            if (k_func(x + i*h1/2, y +j*h2/2) == 1){
                sum += 1;
            }
        }
    }
    if (sum == 0){
        return 0;
    } else if (sum == 4) {
        return 1;
    } else if (sum == 1) {
        double l_x = 0, l_y = 0;
        
        if (x == 3){
            if (y > 0){
                double a = intersection_x(x-h1/2,y-h2/2, y+h2/2)-(y-h2/2);
                double b = intersection_x(x,y-h2/2, y+h2/2)-(y-h2/2);
                return ((a+b)/2)/h2;
            } else if (y < 0) {
                double a = (y+h2/2) - intersection_x(x-h1/2,y-h2/2, y+h2/2);
                double b = (y+h2/2) - intersection_x(x,y-h2/2, y+h2/2);
                return ((a+b)/2)/h2;
            }            
        }
        
        if (y > 0){
            l_x = (x + h1/2 - intersection_y(y - h2/2))/h1;
            l_y = (intersection_x(x + h1/2, y - h2/2, y + h2/2) - (y - h2/2))/h2;
        } else if (y < 0) {
            l_x = (x + h1/2 - intersection_y(y + h2/2))/h1;
            l_y = ((intersection_x(x + h1/2, y - h2/2, y + h2/2) - (y + h2/2)))*(-1)/h2;
        }
        
        return l_x*l_y/2;
        
    } else if (sum == 2) {
        if (x == 3){
            return 0.5;
        } else if (x == 1) {
            return ((x + h1/2) - intersection_y(y + h2/2))/h1;
        } else {
            if (y > 0){
                return (intersection_x(x-h1/2,y-h2/2, y+h2/2)-(y-h2/2)+
                    intersection_x(x+h1/2,y-h2/2, y+h2/2)-(y-h2/2))/(2*h2);
            } else if (y < 0) {
                return ((y+h2/2) - intersection_x(x-h1/2,y-h2/2, y+h2/2)+
                    (y+h2/2) - intersection_x(x+h1/2,y-h2/2, y+h2/2))/(2*h2);
            }
        }
    } else if (sum == 3) {
        if (y > 0){
            double a = (y+h2/2) - intersection_x(x-h1/2,y-h2/2, y+h2/2);
            double b = intersection_y(y+h2/2)-(x-h1/2);
            return 1 - (a/h2)*(b/h1)*0.5;
        } else if (y < 0) {
            double a = (y-h2/2)*(-1) + intersection_x(x-h1/2,y-h2/2, y+h2/2);
            double b = intersection_y(y-h2/2)-(x-h1/2);
            return 1 - (a/h2)*(b/h1)*0.5;
        }
    }
    return 0;
}

void A_B_value(double x, double y, double h1, double h2, double *res){
    if (k_func(x - h1/2, y - h2/2) == 1){
        if (k_func(x - h1/2, y + h2/2) == 1){
            res[0] = 1;
        } else {
            double l = intersection_x(x - h1/2, y - h2/2, y + h2/2) - (y - h1/2);
            res[0] = l/h2 + (1 - l/h2)/eps;
        }
        if (k_func(x + h1/2, y - h2/2) == 1){
            res[1] = 1;
        } else {
            //impossible
        }
    } else {
        if (k_func(x - h1/2, y + h2/2) == 1){
            double l =  (y + h1/2) - intersection_x(x - h1/2, y - h2/2, y + h2/2);
            res[0] = l/h2 + (1 - l/h2)/eps;
        } else {
            res[0] = 1/eps;
        }
        if (k_func(x + h1/2, y - h2/2) == 1){
            double l = intersection_y(y - h2/2) - (x - h1/2);
            res[1] = l/h2 + (1 - l/h2)/eps;
        } else {
            res[1] = 1/eps;
        }
    }
    res[0]/=(h1*h2);
    res[1]/=(h1*h2);
}

void initial_values(double** A, double** B, double** F, int n, int m, double h1, double h2, int proc_x, int proc_y){//вычисление коэффициентов уравнений
    double add_x = proc_x * m;
    double add_y = proc_y * n;
    for (int i = 0; i <= n+1; i++){
        for (int j = 0; j <= m+1; j++){
            double x = P[0][0] + (j + add_x)*h1;
            double y = P[1][0] + (i + add_y)*h2;
            F[i][j] =  F_value(x, y, h1, h2);
            double res[] = {0,0};
            A_B_value(x,y,h1,h2, res);
            A[i][j] = res[0];
            B[i][j] = res[1];
        }
    }
}

void update_borders(int X, int Y, int proc_x, int proc_y, int m, int n,\
                  double** borders, double** arr, int tag) {
    int i, j;
    
    MPI_Status status;
    MPI_Request r1, r2, r3, r4;

    //пересылка налево
    if (proc_x > 0) {
        for (i = 1; i <= n; i++)
            borders[0][i - 1] = arr[i][1];
        int dest = proc_y * X  + proc_x - 1;
        MPI_Isend(borders[0], n, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &r1);
    }

    //пересылка вниз
    if (proc_y > 0) {
        for (j = 1; j <= m; j++)
            borders[1][j - 1] = arr[1][j];
        int dest = (proc_y - 1) * X  + proc_x;
        MPI_Isend(borders[1], m, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &r2);
    }

    //пересылка направо
    if (proc_x < X - 1) {
        for (i = 1; i <= n; i++)
            borders[2][i - 1] = arr[i][m];
        int dest = proc_y * X  + proc_x + 1;
        MPI_Isend(borders[2], n, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &r3);
    }


    //пересылка вверх
    if (proc_y < Y - 1) {
        for (j = 1; j <= m; j++)
            borders[3][j - 1] = arr[n][j];
        int dest = (proc_y + 1) * X  + proc_x;
        MPI_Isend(borders[3], m, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD, &r4);
    }

    //получение слева
    if (proc_x > 0) {
        int sender = proc_y * X  + proc_x - 1;
        MPI_Recv(borders[4], n, MPI_DOUBLE, sender, tag, MPI_COMM_WORLD, &status);

        for (i = 1; i <= n; i++)
            arr[i][0] = borders[4][i - 1];
    }

    //получение сверху
    if (proc_y < Y - 1) {
        int sender = (proc_y + 1) * X  + proc_x;
        MPI_Recv(borders[5], m, MPI_DOUBLE, sender, tag, MPI_COMM_WORLD, &status);

        for (j = 1; j <= m; j++)
            arr[n+1][j] = borders[5][j - 1];
    }

    //получение справа
    if (proc_x < X - 1) {
        int sender = proc_y * X  + proc_x + 1;
        MPI_Recv(borders[6], n, MPI_DOUBLE, sender, tag, MPI_COMM_WORLD, &status);

        for (i = 1; i <= n; i++)
            arr[i][m+1] = borders[6][i - 1];
    }


    //получение снизу
    if (proc_y > 0) {
        int sender = (proc_y - 1) * X  + proc_x;
        MPI_Recv(borders[7], m, MPI_DOUBLE, sender, tag, MPI_COMM_WORLD, &status);

        for (j = 1; j <= m; j++)
            arr[0][j] = borders[7][j - 1];
    }

    if (proc_x > 0) MPI_Wait(&r1, &status);
    if (proc_y > 0) MPI_Wait(&r2, &status);
    if (proc_x < X - 1) MPI_Wait(&r3, &status);
    if (proc_y < Y - 1) MPI_Wait(&r4, &status);
}



int main (int argc, char** argv)
{   
    double start, end;
    int rank, size;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int N, M;
    //int N = 10, M = 10;
    N = atoi(argv[2]);
    M = atoi(argv[1]);
    int x_grid = atoi(argv[3]), y_grid = atoi(argv[4]);
    //int x_grid = 1, y_grid = 1;
    int n_proc = x_grid*y_grid;
    int threads = atoi(argv[5]);
    
    
    
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    
    //scanf("%d %d", &M, &N);
    printf("M = %d, N = %d, proc = %d, rank = %d\n", M, N, n_proc, rank);
    int m = (M-2)/x_grid, n = (N-2)/y_grid;
    int proc_x = rank % x_grid;
    int proc_y = rank / x_grid;
    
    
    double h1 = (P[0][1]-P[0][0])/(M-1), h2 = (P[1][1]-P[1][0])/(N-1);//ширина разбиения
    
    //выделение памяти под массивы коэффициентов а b и значений правой части
    
    double **A = calloc(n + 2, sizeof(double*));
    for (int i = 0; i <= n + 1; i++){
        A[i] = calloc(m + 2, sizeof(double));
    }
    
    double **B = calloc(n + 2, sizeof(double*));
    for (int i = 0; i <= n + 1; i++){
        B[i] = calloc(m + 2, sizeof(double));
    }
    
    double **F = calloc(n + 2, sizeof(double*));
    for (int i = 0; i <= n + 1; i++){
        F[i] = calloc(m + 2, sizeof(double));
    }
    
    //заполнение начальный данных
    //initial_values(A, B, F, n, m, h1, h2, proc_x, proc_y);
    double add_x = proc_x * m;
    double add_y = proc_y * n;
    for (int i = 0; i <= n + 1; i++){//вычисление коэффициентов уравнений
        for (int j = 0; j <= m + 1; j++){
            double x = P[0][0] + (j+add_x)*h1;
            double y = P[1][0] + (i+add_y)*h2;
            //printf("x %4.2lf y %4.2lf\n", x, y);
            if (k_func(x - h1/2, y - h2/2) == 1){
                if (k_func(x - h1/2, y + h2/2) == 1){
                    A[i][j] = 1;
                } else {
                    double l = intersection_x(x - h1/2, y - h2/2, y + h2/2) - (y - h1/2);
                    A[i][j] = l/h2 + (1 - l/h2)/eps;
                }
                if (k_func(x + h1/2, y - h2/2) == 1){
                    B[i][j] = 1;
                } else {
                    //impossible
                }
            } else {
                if (k_func(x - h1/2, y + h2/2) == 1){
                    double l =  (y + h1/2) - intersection_x(x - h1/2, y - h2/2, y + h2/2);
                    A[i][j] = l/h2 + (1 - l/h2)/eps;
                } else {
                    A[i][j] = 1/eps;
                }
                if (k_func(x + h1/2, y - h2/2) == 1){
                    double l = intersection_y(y - h2/2) - (x - h1/2);
                    B[i][j] = l/h2 + (1 - l/h2)/eps;
                } else {
                    B[i][j] = 1/eps;
                }
            }
            F[i][j] = F_value(x, y, h1, h2);
            A[i][j] /= h1*h1;
            B[i][j] /= h2*h2;
        }
    }
    
    //printf("%d made init\n", rank);
    //printf("F:\n");
    //print_2d(F, n+2, m+2);
    //printf("A:\n");
    //print_2d(A, n+2, m+2);
    //printf("B:\n");
    //print_2d(B, n+2, m+2);
    
    
    
    double delta = 0.000001;
    double difference = 1;
    //выделение памяти под массивы для хранения функции и невязки
    
    double **W = calloc(n + 2, sizeof(double*));
    for (int i = 0; i <= n + 1; i++){
        W[i] = calloc(m + 2, sizeof(double));
    }
    
    double **R = calloc(n + 2, sizeof(double*));
    for (int i = 0; i <= n + 1; i++){
        R[i] = calloc(m + 1, sizeof(double));
    }
    
    int k = 0;
    double rr = 0, Arr = 0, scal = 0;
    double loc_rr = 0, loc_Arr = 0, loc_scal = 0;
    double step;
    
    double **borders = calloc(8, sizeof(double*));
    for (int i = 0; i < 8; i++){
        if (i%2 == 0){
            borders[i] = calloc(n, sizeof(double));
        }else{
            borders[i] = calloc(m, sizeof(double));
        }
    }
    
    omp_set_num_threads(threads);
    
    //printf("%d start descent\n", rank);
    for (; (difference > delta); k++){  //метод спуска
        
        int i, j;
        
        
        #pragma omp parallel for private(i,j) reduction(+:loc_rr)
        for (i = 1; i <= n; i++){
            for (j = 1; j <= m; j++){
                R[i][j] = operator(A,B,W,i,j,h1,h2) - F[i][j];
                loc_rr += R[i][j] * R[i][j];
            }
        }
        
        
        //printf("%d after first omp\n", rank);
        MPI_Barrier(MPI_COMM_WORLD);

        update_borders(x_grid, y_grid, proc_x, proc_y, m, n, borders, R, 1);
        
        //printf("proc %d updated borders r\n", rank);
        #pragma omp parallel for reduction(+: loc_Arr) private(i,j) shared(R)
        for (i = 1; i <= n; i++){
            for (j = 1; j <= m; j++){
                loc_Arr += R[i][j] * operator(A,B,R,i,j,h1,h2);
            }
        }
        //printf("proc %d before reduce", rank);
        //printf("%d after second omp\n", rank);
        //printf("iter %d, rr %lf, Arr %lf\n", k,loc_rr,loc_Arr);
        MPI_Allreduce(&loc_rr, &rr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&loc_Arr, &Arr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        //printf("iter %d after reduce, rr %lf, Arr %lf\n", k,rr,Arr);

        
        if (rank == 0){
            step = rr/Arr;
            //printf("counted step = %lf, rr = %lf, Arr = %lf \n", step, rr, Arr);
        }
        MPI_Bcast(&step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //printf("step %lf rr %lf arr %lf", step, rr, Arr);
        #pragma omp parallel for reduction(+:loc_scal) private(i,j) shared(W)
        for (i = 1; i <= n; i++){
            for (j = 1; j <= m; j++){
                double tmp = step*R[i][j];
                loc_scal += tmp*tmp;
                W[i][j] -= tmp;
            }
        }
        //printf("%d after third omp\n", rank);
        MPI_Reduce(&loc_scal, &scal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        
        

        update_borders(x_grid, y_grid, proc_x, proc_y, m, n, borders, W, 2);
        //printf("proc %d updated borders r", rank);
        //printf("scal = %lf\n", scal);
        if (rank == 0){ 
            difference = sqrt(scal);
            //printf("iter = %d, diffference = %lf, step = %lf\n", k, difference, step);
        }
        MPI_Bcast(&difference, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        loc_rr = 0;
        loc_Arr = 0;
        loc_scal = 0;
        rr = 0;
        Arr = 0;
        scal = 0;
    }



    
    //print_2d(R, N+1, M+1);
    free(A);
    free(B);
    free(F);
    free(W);
    free(R);
    free(borders);
    MPI_Barrier(MPI_COMM_WORLD);
    end = MPI_Wtime();

    

    if (rank == 0){
        printf("%d iterations, difference = %lf\n", k, difference);
        printf("time = %fs\n", end - start);
    }
    MPI_Finalize();
    return 0;
}

