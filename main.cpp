#define __rose_lt(x,y) ((x)<(y)?(x):(y))
#define __rose_gt(x,y) ((x)>(y)?(x):(y))
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <omp.h>

#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#define SIZE 15000

void f1(int **A, int T, int N);
void f2(int **A, int T, int N);
void f3(int **A, int T, int N);
void f4(int **A, int T, int N);
int **A = new int *[SIZE], **B = new int *[SIZE], **C = new int *[SIZE];

int main(int argc, char* argv[]){

    for(int i = 0 ; i < SIZE ; i++) {
        A[i] = new int[SIZE];
        B[i] = new int[SIZE];
        C[i] = new int[SIZE];
        for (int j = 0; j < SIZE; j++) {

            A[i][j] = (int) (rand() % 1000000);
            B[i][j] = A[i][j];
            C[i][j] = A[i][j];

        }
    }
    auto start = std::chrono::high_resolution_clock::now();
    f1( A, SIZE, SIZE);
    auto finish = std::chrono::high_resolution_clock::now();
    auto microseconds = std::chrono::duration_cast<std::chrono::microseconds>(finish-start);
    printf("Time taken finish in milliseconds : %f - original code\n", microseconds.count()/1000.0);

    auto start3 = std::chrono::high_resolution_clock::now();
    f4( C, SIZE, SIZE);
    auto finish3  = std::chrono::high_resolution_clock::now();
    auto microseconds3 = std::chrono::duration_cast<std::chrono::microseconds>(finish3-start3);
    printf("Time taken finish in milliseconds : %f - chill without lb & ub code\n", microseconds3.count()/1000.0);

    auto start1 = std::chrono::high_resolution_clock::now();
    f2( B, SIZE, SIZE);
    auto finish1 = std::chrono::high_resolution_clock::now();
    auto microseconds1 = std::chrono::duration_cast<std::chrono::microseconds>(finish1-start1);
    printf("Time taken finish in milliseconds : %f - chill transformation\n", microseconds1.count()/1000.0);


    auto start5 = std::chrono::high_resolution_clock::now();
    f4( C, SIZE, SIZE);
    auto finish5  = std::chrono::high_resolution_clock::now();
    auto microseconds5 = std::chrono::duration_cast<std::chrono::microseconds>(finish5-start5);
    printf("Time taken finish in milliseconds : %f - chill without lb & ub code\n", microseconds5.count()/1000.0);
    auto start2 = std::chrono::high_resolution_clock::now();
    f3( C, SIZE, SIZE);
    auto finish2 = std::chrono::high_resolution_clock::now();
    auto microseconds2 = std::chrono::duration_cast<std::chrono::microseconds>(finish2-start2);
    printf("Time taken finish in milliseconds : %f - pluto code\n", microseconds2.count()/1000.0);




    return 0;
}
void f1( int **A , int T , int N)
{
//    int t,i ;
//    for( t = 0 ; t < T-1 ; t++) {
//        for ( i = 2 ; i < N-3 ; i++)
//            A[t+1][i] = (A[t][i-2] + A[t][i]+A[t][i+2])/3 ;   // S
//
//    }

    int t,i ;
    for( t = 0 ; t < T-1 ; t++) {
        for ( i = 2 ; i < N-3 ; i++)
            A[t+1][i] = (A[t][i-2] + A[t][i] + A[t][i+2])/3 ;   // S

    }


}

void f2(int **A,int T,int N)
{
//    int t8;
//    int t6;
//    int t4;
//    int t2;
//    int t;
//    int i;
//    if (2 <= T && 6 <= N)
//
//        for (t2 = 2; t2 <= 2 * T + N - 8; t2 += 1)
//            for (t4 = __rose_gt((__rose_gt(-t2 - 22,-N + 4)),t2 - 2 * N - 2); t4 <= __rose_lt((__rose_lt((__rose_lt(2 * T - 6,t2 + 6)),4 * T - t2 - 8)),2 * T + N - 12); t4 += 1)
//                for (t6 = __rose_gt((__rose_gt(-t2 - 11,t4)),t2 - 2 * N + 8); t6 <= __rose_lt((__rose_lt((__rose_lt(2 * T - 6,-t2 + 4 * T - 8)),t4 + 11)),t2 + 7); t6 += 1)
//                    if ((t4 + N - 4) % 12 == 0 && (t2 - 2) % 12 == 0)
//                        #pragma omp parallel
//                        for (t8 = __rose_gt((__rose_gt(t2 + (-t6 - t2) % 4,t6 + 4 + (-t6 - (t6 + 4)) % 4)),-t6); t8 <= __rose_lt((__rose_lt(t6 + 2 * N - 8,t2 + 11)),-t6 + 4 * T - 8); t8 += 4)
//// S
//                            A[(t8 + t6) / 4 + 1][(t8 - t6) / 2] = (A[(t8 + t6) / 4][(t8 - t6) / 2 - 2] + A[(t8 + t6) / 4][(t8 - t6) / 2] + A[(t8 + t6) / 4][(t8 - t6) / 2 + 2]) / 3;

//    int t4;
//    int t2;
//    int t;
//    int i;
//    if (6 <= N)
//        //#pragma omp parallel for
//        for (t2 = 0; t2 <= 4 * T - 8; t2 += 4)
//            for (t4 = -(-(t2 + 4) / 2); t4 <= (t2 + 2 * N - 8) / 2; t4 += 1)
//// S
//                A[t2 / 4 + 1][(2 * t4 - t2) / 2] = (A[t2 / 4][(2 * t4 - t2) / 2 - 2] + A[t2 / 4][(2 * t4 - t2) / 2] + A[t2 / 4][(2 * t4 - t2) / 2 + 2]) / 3;


    int t8;
    int t6;
    int t4;
    int t2;
    int t;
    int i;

    int lb, ub, lbp, ubp, lb2, ub2;
    register int lbv, ubv;
    if (2 <= T && 6 <= N) {
//#pragma omp parallel for

        for (t2 = __rose_gt(0,-(-(-N + 36) / 32)); t2 <= __rose_lt((T + 13) / 8,(N + 18 * T + 215) / 144); t2 += 1){
            //lbp = __rose_gt((__rose_gt(1,-(-(-T + 10 * t2 - 12) / 4))),-(-t2 / 2));
            //ubp = __rose_lt((__rose_lt((2 * T + N + 23) / 32,(N + 32 * t2 - 4) / 32)),(N + 16 * t2 + 11) / 32);
#pragma omp parallel for private(t6, t8)
                for (t4 = __rose_gt((__rose_gt(1,-(-(-T + 10 * t2 - 12) / 4))),-(-t2 / 2)); t4 <= __rose_lt((__rose_lt((2 * T + N + 23) / 32,(N + 32 * t2 - 4) / 32)),(N + 16 * t2 + 11) / 32); t4 += 1){
//#pragma omp for ordered schedule(static)
//#pragma omp parallel for
                    for (t6 = __rose_gt((__rose_gt((__rose_gt(16 * t2 - 16 * t4 - 14,0)),-(-(32 * t4 - N - 27) / 2))),8 * t2 - 15); t6 <= __rose_lt((__rose_lt((__rose_lt(16 * t4 - 1,T - 2)),8 * t2)),(32 * t2 - 32 * t4 + N - 4) / 2); t6 += 1) {
#pragma ivdep
#pragma vector always
                        //lbv=__rose_gt((__rose_gt(32 * t4 - 2 * t6 - 31,32 * t4 + 2 * t6 - 32 * t2)),2);
                        //ubv=__rose_lt((__rose_lt(32 * t4 - 2 * t6,32 * t4 + 2 * t6 - 32 * t2 + 31)),N - 4);

                        for (t8 =__rose_gt((__rose_gt(32 * t4 - 2 * t6 - 31,32 * t4 + 2 * t6 - 32 * t2)),2) ; t8 <= __rose_lt((__rose_lt(32 * t4 - 2 * t6,32 * t4 + 2 * t6 - 32 * t2 + 31)),N - 4); t8 += 1)
                            A[t6 + 1][t8] = (A[t6][t8 - 2] + A[t6][t8] + A[t6][t8 + 2]) / 3;
                    }
                }
        }
    }


}

void f4(int **A,int T,int N) {
    int t8;
    int t6;
    int t4;
    int t2;
    int t;
    int i;

    int lb, ub, lbp, ubp, lb2, ub2;
    register int lbv, ubv;
    if (2 <= T && 6 <= N) {
//#pragma omp parallel for
        for (t2 = __rose_gt(0, -(-(-N + 36) / 32)); t2 <= __rose_lt((T + 13) / 8, (N + 18 * T + 215) / 144); t2 += 1) {
#pragma omp parallel for private(t6, t8)
            for (t4 = __rose_gt((__rose_gt(1, -(-(-T + 10 * t2 - 12) / 4))), -(-t2 / 2)); t4 <= __rose_lt(
                    (__rose_lt((2 * T + N + 23) / 32, (N + 32 * t2 - 4) / 32)), (N + 16 * t2 + 11) / 32); t4 += 1) {
//#pragma omp for ordered schedule(static)
//#pragma omp parallel for
                for (t6 = __rose_gt((__rose_gt((__rose_gt(16 * t2 - 16 * t4 - 14, 0)), -(-(32 * t4 - N - 27) / 2))),
                                    8 * t2 - 15); t6 <= __rose_lt((__rose_lt((__rose_lt(16 * t4 - 1, T - 2)), 8 * t2)),
                                                                  (32 * t2 - 32 * t4 + N - 4) / 2); t6 += 1) {
#pragma ivdep
#pragma vector always

                    for (t8 = __rose_gt((__rose_gt(32 * t4 - 2 * t6 - 31, 32 * t4 + 2 * t6 - 32 * t2)), 2); t8 <=
                                                                                                            __rose_lt(
                                                                                                                    (__rose_lt(
                                                                                                                            32 *
                                                                                                                            t4 -
                                                                                                                            2 *
                                                                                                                            t6,
                                                                                                                            32 *
                                                                                                                            t4 +
                                                                                                                            2 *
                                                                                                                            t6 -
                                                                                                                            32 *
                                                                                                                            t2 +
                                                                                                                            31)),
                                                                                                                    N -
                                                                                                                    4); t8 += 1)
                        A[t6 + 1][t8] = (A[t6][t8 - 2] + A[t6][t8] + A[t6][t8 + 2]) / 3;
                }
            }
        }
    }

}

    void f3(int **A, int T, int N) {
        int t1, t2, t3, t4;
        int lb, ub, lbp, ubp, lb2, ub2;
        register int lbv, ubv;
/* Start of CLooG code */
        if ((N >= 6) && (T >= 2)) {
            for (t1 = 0; t1 <= floord(3 * T + N - 10, 32); t1++) {
                lbp = max(ceild(2 * t1, 3), ceild(32 * t1 - T + 2, 32));
                ubp = min(min(floord(2 * T + N - 8, 32), floord(64 * t1 + N + 58, 96)), t1);
#pragma omp parallel for private(lbv,ubv,t3,t4)
                for (t2 = lbp; t2 <= ubp; t2++) {
                    for (t3 = max(ceild(32 * t2 - N + 4, 2), 32 * t1 - 32 * t2);
                         t3 <= min(min(T - 2, 16 * t2 + 14), 32 * t1 - 32 * t2 + 31); t3++) {
                        lbv = max(32 * t2, 2 * t3 + 2);
                        ubv = min(32 * t2 + 31, 2 * t3 + N - 4);
#pragma ivdep
#pragma vector always
                        for (t4 = lbv; t4 <= ubv; t4++) {
                            A[t3 + 1][(-2 * t3 + t4)] =
                                    (A[t3][(-2 * t3 + t4) - 2] + A[t3][(-2 * t3 + t4)] + A[t3][(-2 * t3 + t4) + 2]) /
                                    3;;
                        }
                    }
                }
            }
        }
    }

/*void mm_parallel(  double A[][SIZE] , double B[][SIZE] , double C[][SIZE])
{

#pragma omp parallel num_threads(8)
    {
#pragma omp for schedule(dynamic)
        for(int i = 0 ; i < SIZE ; i++) {
            for(int j=0; j < SIZE ; j++) {
                C[i][j] = 0 ;
                for(int k = 0 ; k < SIZE ; k++)
                    C[i][j] += A[i][k]*B[k][j] ;

            }

        }
    }
}

void mm_normal(  double A[][SIZE] , double B[][SIZE] , double C[][SIZE])
{
    for(int i = 0 ; i < SIZE ; i++) {
        for(int j=0; j < SIZE ; j++) {

            C[i][j] = 0 ;

            for(int k = 0 ; k < SIZE ; k++)
                C[i][j] += A[i][k]*B[k][j] ;
        }
    }
}
void mm_tiled_parallel( double A[][SIZE] , double B[][SIZE] , double C[][SIZE] ,int block_size)
{

#pragma omp parallel num_threads(8)
    {
#pragma omp for schedule(dynamic)
        for(int jj = 0; jj < SIZE; jj += block_size){

            for(int kk = 0; kk < SIZE; kk += block_size){

                for(int i = 0; i < SIZE; i++){

                    for(int j = jj; j < (jj + block_size); j++){

                        for(int k = kk; k < (kk + block_size); k++){
                            C[i][j] += A[i][k] * B[k][j];
                        }
                    }
                }
            }
        }  }*/
//}