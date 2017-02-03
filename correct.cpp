//#define __rose_lt(x,y) ((x)<(y)?(x):(y))
//#define __rose_gt(x,y) ((x)>(y)?(x):(y))
//#include<iostream>
//#include<vector>
//#include<ctime>
//#include<cstdlib>
//#include<chrono>
//#include<omp.h>
//
//using namespace std ;
//
//void f( int **A , int T , int N);
//void f1( int **A , int T , int N);
//void f2( int **A , int T , int N);
//
//int main() {
//
//
//    int N = 5000 , T = 5000 ;
//    int** ary = new int*[T];
//    int** ary2 = new int*[T] ;
//    for(int i = 0; i < T; ++i) {
//        ary[i] = new int[N];
//        ary2[i] = new int[N] ;
//    }
//
//
//    srand( time(NULL)) ;
//
//    for(int i = 0 ; i < T ; i++)
//        for(int j = 0 ; j < N ; j++)
//            ary[i][j] = i+j ;
//
//    /*
//    for(int i = 0 ; i < T ; i++) {
//        	  for(int j = 0 ; j < N ; j++)
//        		    cout << ary[i][j] << " " ;
//        	  cout << endl ;
//    }
//    */
//    for(int i = 0 ; i < T ; i++)
//        for(int j = 0 ; j < N ; j++)
//            ary2[i][j] = ary[i][j] ;
//
//
//
//
//    auto t0 = std::chrono::high_resolution_clock::now();
//    f(ary,T,N) ;
//    auto t1 = chrono::high_resolution_clock::now();
//
//    auto microseconds = chrono::duration_cast<std::chrono::microseconds>(t1-t0);
//    cout << "time to execute f " << microseconds.count()/1000.0 << "ms\n";
//
//
//    auto t2 = std::chrono::high_resolution_clock::now();
//    f2(ary2,T,N) ;
//    auto t3 = chrono::high_resolution_clock::now();
//
//    auto microseconds1 = chrono::duration_cast<std::chrono::microseconds>(t3-t2);
//    cout << "time to execute f1 " << microseconds1.count()/1000.0 << "ms\n";
//
//}
//
//
//
//
//void f( int **A , int T , int N)
//{
//
//
//    int t,i ;
//    for( t = 0 ; t < T-1 ; t++)
//        for ( i = 2 ; i < N-3 ; i++)
//            A[t+1][i] = (A[t][i-2] + A[t][i] + A[t][i+2])/3 ;   // S
//
//
//
//
//
//}
//
//void f2(int **A,int T,int N)
//{
//    int t8;
//    int t6;
//    int t4;
//    int t2;
//    int t;
//    int i;
//    if (2 <= T && 6 <= N) {
//#pragma omp parallel for private(t2, t4, t6, t8)
//        for (t2 = __rose_gt(0,-(-(-N + 36) / 32)); t2 <= __rose_lt((T + 13) / 8,(N + 18 * T + 215) / 144); t2 += 1)
////#pragma omp parallel for
//            for (t4 = __rose_gt((__rose_gt(1,-(-(-T + 10 * t2 - 12) / 4))),-(-t2 / 2)); t4 <= __rose_lt((__rose_lt((2 * T + N + 23) / 32,(N + 32 * t2 - 4) / 32)),(N + 16 * t2 + 11) / 32); t4 += 1)
////#pragma omp for ordered schedule(static)
////#pragma omp parallel for
//                    for (t6 = __rose_gt((__rose_gt((__rose_gt(16 * t2 - 16 * t4 - 14,0)),-(-(32 * t4 - N - 27) / 2))),8 * t2 - 15); t6 <= __rose_lt((__rose_lt((__rose_lt(16 * t4 - 1,T - 2)),8 * t2)),(32 * t2 - 32 * t4 + N - 4) / 2); t6 += 1) {
//                        for (t8 = __rose_gt((__rose_gt(32 * t4 - 2 * t6 - 31,32 * t4 + 2 * t6 - 32 * t2)),2); t8 <= __rose_lt((__rose_lt(32 * t4 - 2 * t6,32 * t4 + 2 * t6 - 32 * t2 + 31)),N - 4); t8 += 1)
//                            A[t6 + 1][t8] = (A[t6][t8 - 2] + A[t6][t8] + A[t6][t8 + 2]) / 3;
//                    }
//    }
//}
//
//void f1(int **A,int T,int N)
//{
//    int t8;
//    int t6;
//    int t4;
//    int t2;
//    int t;
//    int i;
//
//    if (2 <= T && 6 <= N) {
//#pragma omp parallel for private(t2,t4,t6,t8) num_threads(3)
//        for (t2 = __rose_gt(0,-(-(-N + 36) / 32)); t2 <= __rose_lt((T + 13) / 8,(N + 18 * T + 215) / 144); t2 += 1) {
//            for (t4 = __rose_gt((__rose_gt(1,-(-(-T + 10 * t2 - 12) / 4))),-(-t2 / 2)); t4 <= __rose_lt((__rose_lt((2 * T + N + 23) / 32,(N + 32 * t2 - 4) / 32)),(N + 16 * t2 + 11) / 32); t4 += 1)
//                for (t6 = __rose_gt((__rose_gt((__rose_gt(16 * t2 - 16 * t4 - 14,0)),-(-(32 * t4 - N - 27) / 2))),8 * t2 - 15); t6 <= __rose_lt((__rose_lt((__rose_lt(16 * t4 - 1,T - 2)),8 * t2)),(32 * t2 - 32 * t4 + N - 4) / 2); t6 += 1)
//                    for (t8 = __rose_gt((__rose_gt(32 * t4 - 2 * t6 - 31,32 * t4 + 2 * t6 - 32 * t2)),2); t8 <= __rose_lt((__rose_lt(32 * t4 - 2 * t6,32 * t4 + 2 * t6 - 32 * t2 + 31)),N - 4); t8 += 1)
//                        A[t6 + 1][t8] = (A[t6][t8 - 2] + A[t6][t8] + A[t6][t8 + 2]) / 3;
//        }
//    }
//}

////
//
//
//
//
//
//
//












#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))
#define __rose_lt(x,y) ((x)<(y)?(x):(y))
#define __rose_gt(x,y) ((x)>(y)?(x):(y))
#include<iostream>
#include<vector>
#include<ctime>
#include<cstdlib>
#include<omp.h>

using namespace std ;

void f1(int **A , int T , int N);
void f3(int **A , int T , int N) ;




int main() {


    int N = 50 , T = 50 ;
//    double ** ary = new double*[T];
//    double ** ary2 = new double*[T] ;
//    for(int i = 0; i < T; ++i) {
//        ary[i] = new double[N];
//        ary2[i] = new double[N] ;
//    }

    int ** ary = new int*[T];
    int ** ary2 = new int*[T] ;
    for(int i = 0; i < T; ++i) {
        ary[i] = new int[N];
        ary2[i] = new int[N] ;
    }


    srand( time(NULL)) ;

    for(int i = 0 ; i < T ; i++)
        for(int j = 0 ; j < N ; j++)
            ary[i][j] = rand()%500 + 1 ;

    /*
    for(int i = 0 ; i < T ; i++) {
        	  for(int j = 0 ; j < N ; j++)
        		    cout << ary[i][j] << " " ;
        	  cout << endl ;
    }
    */
    for(int i = 0 ; i < T ; i++)
        for(int j = 0 ; j < N ; j++)
            ary2[i][j] = ary[i][j] ;



    f1(ary,T,N) ;
    f3(ary2,T,N) ;


    for(int i = 0 ; i < T ; i++) {
        for(int j = 0 ; j < N ; j++)
            cout << ary[i][j]-ary2[i][j] << " " ;
        cout << endl ;
    }



}



void f1( int **A , int T , int N)
{


    int t,i ;
    for( t = 0 ; t < T-1 ; t++)
        for ( i = 2 ; i < N-3 ; i++)
            A[t+1][i] = (A[t][i-2] + A[t][i] + A[t][i+2])/3 ;   // S


}

void f0(int **A , int T , int N) {

    int i, j;
    int t;
    for (t = 1; t <= T; t++) {
        for (i = 1; i < N - 1; i++)
            for (j = 2; j < N - 1; j++)
                A[i][j] = 0.25 * (A[i - 1][j] + A[i + 1][j] + A[i][j - 1] + A[i][j + 1]);


    }
}


void f2(double **A , int T , int N)
{


    int t10;
    int t8;
    int t6;
    int t4;
    int t2;
    int i;
    int j;
    int t;
    if (1 <= T && 4 <= N)
        for (t2 = 2; t2 <= (N + 2 * T + 60) / 32; t2 += 1)
            for (t4 = __rose_gt(-(-t2 / 2),-(-(32 * t2 - T - 31) / 32)); t4 <= __rose_lt((__rose_lt((32 * t2 + N + 29) / 64,(N + T + 29) / 32)),t2 - 1); t4 += 1)
                for (t6 = __rose_gt(32 * t4 - N - 29,32 * t2 - 32 * t4 - 31); t6 <= __rose_lt((__rose_lt(32 * t2 - 32 * t4,32 * t4 - 1)),T); t6 += 1)
                    for (t8 = __rose_gt(32 * t4 - t6 - 31,1); t8 <= __rose_lt(N - 2,32 * t4 - t6); t8 += 1)
                        for (t10 = 2; t10 <= N - 2; t10 += 1)
                            A[t8][t10] = (0.25 * (A[t8 - 1][t10] + A[t8 + 1][t10] + A[t8][t10 - 1] + A[t8][t10 + 1]));


}

void f31( int **A , int T , int N) {


    int t10;
    int t8;
    int t6;
    int t4;
    int t2;
    int i;
    int j;
    int t;
    int lb,ub;
    if (1 <= T && 4 <= N)
        for (t2 = 2; t2 <= (N + 2 * T + 60) / 32; t2 += 1) {
            lb = __rose_gt(-(-t2 / 2), -(-(32 * t2 - T - 31) / 32));
            ub = __rose_lt((__rose_lt((32 * t2 + N + 29) / 64, (N + T + 29) / 32)), t2 - 1);
#pragma omp parallel for private(t4,t6,t8,t10)
                        for (t4 = lb;t4 <= ub; t4 += 1) {

            //for (t4 = __rose_gt(-(-t2 / 2),-(-(32 * t2 - T - 31) / 32)); t4 <= __rose_lt((__rose_lt((32 * t2 + N + 29) / 64,(N + T + 29) / 32)),t2 - 1); t4 += 1){
                for (t6 = __rose_gt(32 * t4 - N - 29, 32 * t2 - 32 * t4 - 31);
                     t6 <= __rose_lt((__rose_lt(32 * t2 - 32 * t4, 32 * t4 - 1)), T); t6 += 1) {
                    for (t8 = __rose_gt(32 * t4 - t6 - 31, 1); t8 <= __rose_lt(N - 2, 32 * t4 - t6); t8 += 1) {
                        for (t10 = 2; t10 <= N - 2; t10 += 1) {
                            A[t8][t10] = (0.25 * (A[t8 - 1][t10] + A[t8 + 1][t10] + A[t8][t10 - 1] + A[t8][t10 + 1]));
                        }
                    }
                }
            }
        }
}

void f6(int **A,int T,int N)
{
    int t8;
    int t6;
    int t4;
    int t2;
    int t;
    int i;
    if (2 <= T && 6 <= N) {
//#pragma omp parallel for
        for (t2 = __rose_gt(0,-(-(-N + 36) / 32)); t2 <= __rose_lt((T + 13) / 8,(N + 18 * T + 215) / 144); t2 += 1)
#pragma omp parallel for private(t6, t8)
                for (t4 = __rose_gt((__rose_gt(1,-(-(-T + 10 * t2 - 12) / 4))),-(-t2 / 2)); t4 <= __rose_lt((__rose_lt((2 * T + N + 23) / 32,(N + 32 * t2 - 4) / 32)),(N + 16 * t2 + 11) / 32); t4 += 1)
//#pragma omp for ordered schedule(static)
//#pragma omp parallel for
                    for (t6 = __rose_gt((__rose_gt((__rose_gt(16 * t2 - 16 * t4 - 14,0)),-(-(32 * t4 - N - 27) / 2))),8 * t2 - 15); t6 <= __rose_lt((__rose_lt((__rose_lt(16 * t4 - 1,T - 2)),8 * t2)),(32 * t2 - 32 * t4 + N - 4) / 2); t6 += 1) {
#pragma ivdep
#pragma vector always
                        for (t8 = __rose_gt((__rose_gt(32 * t4 - 2 * t6 - 31,32 * t4 + 2 * t6 - 32 * t2)),2); t8 <= __rose_lt((__rose_lt(32 * t4 - 2 * t6,32 * t4 + 2 * t6 - 32 * t2 + 31)),N - 4); t8 += 1)
                            A[t6 + 1][t8] = (A[t6][t8 - 2] + A[t6][t8] + A[t6][t8 + 2]) / 3;
                    }
    }
}

void f27(int **A,int T,int N)
{
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
        lbp = __rose_gt((__rose_gt(1,-(-(-T + 10 * t2 - 12) / 4))),-(-t2 / 2));
        ubp = __rose_lt((__rose_lt((2 * T + N + 23) / 32,(N + 32 * t2 - 4) / 32)),(N + 16 * t2 + 11) / 32);
#pragma omp parallel for private(t6, t8,lbv,ubv) schedule(dynamic)
        for (t4 = lbp; t4 <= ubp; t4 += 1)
//#pragma omp for ordered schedule(static)
//#pragma omp parallel for
            for (t6 = __rose_gt((__rose_gt((__rose_gt(16 * t2 - 16 * t4 - 14, 0)), -(-(32 * t4 - N - 27) / 2))),
                                8 * t2 - 15); t6 <= __rose_lt((__rose_lt((__rose_lt(16 * t4 - 1, T - 2)), 8 * t2)),
                                                              (32 * t2 - 32 * t4 + N - 4) / 2); t6 += 1) {
                lbv = __rose_gt((__rose_gt(32 * t4 - 2 * t6 - 31, 32 * t4 + 2 * t6 - 32 * t2)), 2);
                ubv = __rose_lt((__rose_lt(32 * t4 - 2 * t6, 32 * t4 + 2 * t6 - 32 * t2 + 31)), N - 4);
#pragma ivdep
#pragma vector always
                for (t8 =lbv ;
                     t8 <=ubv ; t8 += 1)
                    A[t6 + 1][t8] = (A[t6][t8 - 2] + A[t6][t8] + A[t6][t8 + 2]) / 3;
            }
        }
    }
}

//128x128
void f300(int **A,int T,int N)
{
    int t8;
    int t6;
    int t4;
    int t2;
    int t;
    int i;
    int lb, ub, lbp, ubp, lb2, ub2;
    register int lbv, ubv;
    if (2 <= T && 6 <= N)
        for (t2 = __rose_gt(0,-(-(-N + 132) / 128)); t2 <= __rose_lt((T + 61) / 32,(N + 66 * T + 3959) / 2112); t2 += 1){
            lbp = __rose_gt((__rose_gt(1,-(-(-T + 40 * t2 - 60) / 16))),-(-t2 / 2));
            ubp = __rose_lt((__rose_lt((2 * T + N + 119) / 128,(N + 128 * t2 - 4) / 128)),(N + 64 * t2 + 59) / 128);
#pragma omp parallel for private(t6, t8,lbv,ubv) schedule(dynamic)
            for (t4 = lbp ; t4 <= ubp; t4 += 1){
                for (t6 = __rose_gt((__rose_gt((__rose_gt(64 * t2 - 64 * t4 - 62,0)),-(-(128 * t4 - N - 123) / 2))),32 * t2 - 63); t6 <= __rose_lt((__rose_lt((__rose_lt(64 * t4 - 1,T - 2)),32 * t2)),(128 * t2 - 128 * t4 + N - 4) / 2); t6 += 1){
                    lbv =  __rose_gt((__rose_gt(128 * t4 - 2 * t6 - 127,128 * t4 + 2 * t6 - 128 * t2)),2);
                    ubv = __rose_lt((__rose_lt(128 * t4 - 2 * t6,128 * t4 + 2 * t6 - 128 * t2 + 127)),N - 4);
#pragma ivdep
#pragma vector always
                    for (t8 =lbv; t8 <= ubv ; t8 += 1){
                        A[t6 + 1][t8] = (A[t6][t8 - 2] + A[t6][t8] + A[t6][t8 + 2]) / 3;
                    }

                }

        }
        }
}
//64x64
void f301(int **A,int T,int N)
{
    int t8;
    int t6;
    int t4;
    int t2;
    int t;
    int i;
    int lb, ub, lbp, ubp, lb2, ub2;
    register int lbv, ubv;
    if (2 <= T && 6 <= N)
        for (t2 = __rose_gt(0,-(-(-N + 68) / 64)); t2 <= __rose_lt((T + 29) / 16,(N + 34 * T + 951) / 544); t2 += 1) {
            lbp = __rose_gt((__rose_gt(1, -(-(-T + 20 * t2 - 28) / 8))), -(-t2 / 2));
            ubp = __rose_lt((__rose_lt((2 * T + N + 55) / 64, (N + 64 * t2 - 4) / 64)), (N + 32 * t2 + 27) / 64);
#pragma omp parallel for private(t6, t8,lbv,ubv) schedule(dynamic)
            for (t4 = lbp; t4 <=ubp ; t4 += 1) {
                for (t6 = __rose_gt((__rose_gt((__rose_gt(32 * t2 - 32 * t4 - 30, 0)), -(-(64 * t4 - N - 59) / 2))),
                                    16 * t2 - 31); t6 <=
                                                   __rose_lt((__rose_lt((__rose_lt(32 * t4 - 1, T - 2)), 16 * t2)),
                                                             (64 * t2 - 64 * t4 + N - 4) / 2); t6 += 1) {
                    lbv =  __rose_gt((__rose_gt(64 * t4 - 2 * t6 - 63, 64 * t4 + 2 * t6 - 64 * t2)), 2);
                    ubv = __rose_lt((__rose_lt(
                            64 *
                            t4 -
                            2 *
                            t6,
                            64 *
                            t4 +
                            2 *
                            t6 -
                            64 *
                            t2 +
                            63)),
                                    N -
                                    4);
#pragma ivdep
#pragma vector always
                    for (t8 = lbv; t8 <= ubv; t8 += 1) {
                        A[t6 + 1][t8] = (A[t6][t8 - 2] + A[t6][t8] + A[t6][t8 + 2]) / 3;
                    }
                }
            }
        }


}
//16x16
void f3(int **A,int T,int N)
{
    int t8;
    int t6;
    int t4;
    int t2;
    int t;
    int i;
    int lb, ub, lbp, ubp, lb2, ub2;
    register int lbv, ubv;
    if (2 <= T && 6 <= N)
        for (t2 = __rose_gt(-(-(-N + 20) / 16),0); t2 <= __rose_lt((T + 5) / 4,(N + 10 * T + 39) / 40); t2 += 1) {
            lbp = __rose_gt((__rose_gt(1, -(-(-T + 5 * t2 - 4) / 2))), -(-t2 / 2));
            ubp = __rose_lt((__rose_lt((2 * T + N + 7) / 16, (N + 16 * t2 - 4) / 16)), (N + 8 * t2 + 3) / 16);
#pragma omp parallel for private(t6, t8,lbv,ubv) schedule(dynamic)
            for (t4 = lbp; t4 <= ubp; t4 += 1) {
                for (t6 = __rose_gt(
                        (__rose_gt((__rose_gt(8 * t2 - 8 * t4 - 6, 0)), -(-(16 * t4 - N - 11) / 2))),
                        4 * t2 - 7); t6 <= __rose_lt((__rose_lt((__rose_lt(8 * t4 - 1, T - 2)), 4 * t2)),
                                                     (16 * t2 - 16 * t4 + N - 4) / 2); t6 += 1) {
                    lbv = __rose_gt((__rose_gt(16 * t4 - 2 * t6 - 15, 16 * t4 + 2 * t6 - 16 * t2)), 2);
                    ubv = __rose_lt((__rose_lt(16 * t4 - 2 * t6, 16 * t4 + 2 * t6 - 16 * t2 + 15)), N - 4);
#pragma ivdep
#pragma vector always
                    for (t8 = lbv; t8 <= ubv; t8 += 1) {
                        A[t6 + 1][t8] = (A[t6][t8 - 2] + A[t6][t8] + A[t6][t8 + 2]) / 3;
                    }
                }
            }
        }
}



//64x32
void f303(int **A,int T,int N)
{
    int t8;
    int t6;
    int t4;
    int t2;
    int t;
    int i;
    if (2 <= T && 6 <= N)
        for (t2 = 1; t2 <= __rose_lt((6 * T + N + 107) / 64,(54 * T + 11 * N + 925) / 576); t2 += 1)
            for (t4 = __rose_gt((__rose_gt(-(-(-T + 16 * t2 - 21) / 8),-(-(2 * t2 - 1) / 3))),-(-(-T + 20 * t2 - 22) / 14)); t4 <= __rose_lt((__rose_lt((__rose_lt((64 * t2 + N - 4) / 64,2 * t2)),(32 * t2 + N + 11) / 48)),(2 * T + N + 23) / 32); t4 += 1)
                for (t6 = __rose_gt((__rose_gt((__rose_gt(32 * t2 - 32 * t4 - 30,0)),-(-(32 * t4 - N - 27) / 2))),16 * t2 - 8 * t4 - 23); t6 <= __rose_lt((__rose_lt((__rose_lt(16 * t4 - 1,T - 2)),16 * t2 - 8 * t4)),(64 * t2 - 64 * t4 + N - 4) / 2); t6 += 1)
                    for (t8 = __rose_gt((__rose_gt(32 * t4 - 2 * t6 - 31,64 * t4 + 2 * t6 - 64 * t2)),2); t8 <= __rose_lt((__rose_lt(32 * t4 - 2 * t6,64 * t4 + 2 * t6 - 64 * t2 + 63)),N - 4); t8 += 1)
                        A[t6 + 1][t8] = (A[t6][t8 - 2] + A[t6][t8] + A[t6][t8 + 2]) / 3;
}
//16x32
void f304(int **A,int T,int N)
{
    int t8;
    int t6;
    int t4;
    int t2;
    int t;
    int i;
    if (2 <= T && 6 <= N)
        for (t2 = __rose_gt(-(-(-N + 4) / 32),-(-(-N + 20) / 16)); t2 <= (3 * T + 22) / 16; t2 += 1)
            for (t4 = __rose_gt((__rose_gt((__rose_gt(-(-(-T + 8 * t2 - 4) / 8),-(-t2 / 3))),-t2)),1); t4 <= __rose_lt((__rose_lt((__rose_lt((T - 4 * t2 + 9) / 4,(8 * t2 + N + 11) / 24)),(16 * t2 + N - 4) / 16)),(2 * T + N + 23) / 32); t4 += 1)
                for (t6 = __rose_gt((__rose_gt((__rose_gt(8 * t2 - 8 * t4 - 6,0)),-(-(32 * t4 - N - 27) / 2))),4 * t2 + 4 * t4 - 11); t6 <= __rose_lt((__rose_lt((__rose_lt(16 * t4 - 1,T - 2)),4 * t2 + 4 * t4)),(16 * t2 - 16 * t4 + N - 4) / 2); t6 += 1)
                    for (t8 = __rose_gt((__rose_gt(32 * t4 - 2 * t6 - 31,16 * t4 + 2 * t6 - 16 * t2)),2); t8 <= __rose_lt((__rose_lt(32 * t4 - 2 * t6,N - 4)),16 * t4 + 2 * t6 - 16 * t2 + 15); t8 += 1)
                        A[t6 + 1][t8] = (A[t6][t8 - 2] + A[t6][t8] + A[t6][t8 + 2]) / 3;
}