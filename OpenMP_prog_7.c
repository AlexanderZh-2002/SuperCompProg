#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

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


int main (int argc, char** argv)
{
    double start, end;
    start = omp_get_wtime();
    int N, M;
    //int N = 180, M = 160;
    N = atoi(argv[2]);
    M = atoi(argv[1]);
    
    //scanf("%d %d", &M, &N);
    printf("M = %d, N = %d\n", M, N);
    
    double h1 = (P[0][1]-P[0][0])/M, h2 = (P[1][1]-P[1][0])/N;//ширина разбиения
    
    //выделение памяти под массивы коэффициентов а b и значений правой части
    
    double **A = calloc(N + 1, sizeof(double*));
    for (int i = 0; i <= N; i++){
        A[i] = calloc(M + 1, sizeof(double));
    }
    
    double **B = calloc(N + 1, sizeof(double*));
    for (int i = 0; i <= N; i++){
        B[i] = calloc(M + 1, sizeof(double));
    }
    
    double **F = calloc(N + 1, sizeof(double*));
    for (int i = 0; i <= N; i++){
        F[i] = calloc(M + 1, sizeof(double));
    }
    
    
    for (int i = 0; i <= N; i++){//вычисление коэффициентов уравнений
        for (int j = 0; j <= M; j++){
            double x = P[0][0] + j*h1;
            double y = P[1][0] + i*h2;
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
    //print_2d(F, N+1, M+1);
    
    double delta = 0.000001;
    double difference = 1;
    
    //выделение памяти под массивы для хранения функции и невязки
    
    double **W = calloc(N + 1, sizeof(double*));
    for (int i = 0; i <= N; i++){
        W[i] = calloc(M + 1, sizeof(double));
    }
    
    double **R = calloc(N + 1, sizeof(double*));
    for (int i = 0; i <= N; i++){
        R[i] = calloc(M + 1, sizeof(double));
    }
    
    int k = 0;
    double rr = 0, Arr = 0, scal = 0;
    double step;
        
    for (; (difference > delta); k++){  //метод спуска
        //printf("k = %d  dif = %10.6lf\n", k, difference);
        #pragma omp parallel for collapse(2) shared(R)
        for (int i = 1; i < N; i++){
            for (int j = 1; j < M; j++){
                R[i][j] = operator(A,B,W,i,j,h1,h2) - F[i][j];
            }
        }
                
        #pragma omp parallel for collapse(2) reduction(+:rr, Arr)
        for (int i = 1; i < N; i++){
            for (int j = 1; j < M; j++){
                rr += R[i][j] * R[i][j];
                Arr += R[i][j] * operator(A,B,R,i,j,h1,h2);
            }
        }
        
        step = rr/Arr;
        //printf("step %lf rr %lf arr %lf", step, rr, Arr);
        #pragma omp parallel for collapse(2) shared(W, R) reduction(+:scal)
        for (int i = 1; i < N; i++){
            for (int j = 1; j < M; j++){
                double tmp = step*R[i][j];
                scal += tmp*tmp;
                W[i][j] -= tmp;
            }
        }
            
        //printf("scal = %lf\n", scal);
        difference = sqrt(scal);
        rr = 0;
        Arr = 0;
        scal = 0;
    }
    
    end = omp_get_wtime();     
    printf("%d iterations, difference = %lf\ntime: %lf s", k, difference, end - start);
    //print_2d(R, N+1, M+1);
    free(A);
    free(B);
    free(F);
    free(W);
    free(R);
}
