/*
 * gcc -O2 -fopenmp -std=c99 -o matmul2 matmul2.c -lm
 * uso: ./matmul2 [N]   (N debe ser potencia de 2, defecto 512)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#define BS 64

static double *mnew(int n){ double*p=malloc((size_t)n*n*8); if(!p){perror("");exit(1);} return p; }
static void mzero(double*C,int n){ memset(C,0,(size_t)n*n*8); }
static void mrand(double*A,int n){ for(int i=0;i<n*n;i++)A[i]=(double)rand()/RAND_MAX; }
static void madd(double*C,const double*A,const double*B,int n){ for(int i=0;i<n*n;i++)C[i]=A[i]+B[i]; }
static void msub(double*C,const double*A,const double*B,int n){ for(int i=0;i<n*n;i++)C[i]=A[i]-B[i]; }

static void extr(double*d,const double*s,int n,int r,int c,int h){
    for(int i=0;i<h;i++) memcpy(d+i*h, s+(r+i)*n+c, h*8); }
static void insr(double*d,const double*s,int n,int r,int c,int h){
    for(int i=0;i<h;i++) memcpy(d+(r+i)*n+c, s+i*h, h*8); }

static void mseq(const double*A,const double*B,double*C,int n){
    mzero(C,n);
    for(int bi=0;bi<n;bi+=BS) for(int bk=0;bk<n;bk+=BS) for(int bj=0;bj<n;bj+=BS){
        int il=bi+BS<n?bi+BS:n, kl=bk+BS<n?bk+BS:n, jl=bj+BS<n?bj+BS:n;
        for(int i=bi;i<il;i++) for(int k=bk;k<kl;k++){
            double a=A[i*n+k];
            for(int j=bj;j<jl;j++) C[i*n+j]+=a*B[k*n+j];
        }
    }
}

static void mseq_naive(const double*A,const double*B,double*C,int n){
    mzero(C,n);
    for(int i=0; i<n; i++){
        for(int k=0; k<n; k++){
            double a = A[i*n+k];
            for(int j=0; j<n; j++){
                C[i*n+j] += a * B[k*n+j];
            }
        }
    }
}

/* submatrices / paraleliza los bloques
   collapse(2) reparte pares (bi,bj) entre hilos.
   Cada par posee una region disjunta de C */
static void blk_par(const double *A, const double *B, double *C, int n, int d)
{
    (void)d;
    mzero(C, n);

    #pragma omp parallel for collapse(2) schedule(dynamic) \
        shared(A, B, C) firstprivate(n) default(none)
    for (int bi = 0; bi < n; bi += BS)
        for (int bj = 0; bj < n; bj += BS)
            for (int bk = 0; bk < n; bk += BS) {
                int il = bi+BS < n ? bi+BS : n;
                int kl = bk+BS < n ? bk+BS : n;
                int jl = bj+BS < n ? bj+BS : n;
                for (int i = bi; i < il; i++)
                    for (int k = bk; k < kl; k++) {
                        double a = A[i*n+k];
                        for (int j = bj; j < jl; j++)
                            C[i*n+j] += a * B[k*n+j];
                    }
            }
}

/* Nucleo Strassen
   hybrid=0 → caso base naif  |  hybrid=1 → caso base por bloques
   Los 7 productos Mi son independientes → sin races en parallel sections. */
static void strassen(const double*A,const double*B,double*C,int n,int d,int hyb){
    if(n<=BS){ hyb ? mseq(A,B,C,n) : mseq_naive(A,B,C,n); return; }
    int h=n/2;
    double *a11=mnew(h),*a12=mnew(h),*a21=mnew(h),*a22=mnew(h),
           *b11=mnew(h),*b12=mnew(h),*b21=mnew(h),*b22=mnew(h),
           *m1=mnew(h),*m2=mnew(h),*m3=mnew(h),*m4=mnew(h),
           *m5=mnew(h),*m6=mnew(h),*m7=mnew(h),*t1=mnew(h),*t2=mnew(h);
    extr(a11,A,n,0,0,h); extr(a12,A,n,0,h,h);
    extr(a21,A,n,h,0,h); extr(a22,A,n,h,h,h);
    extr(b11,B,n,0,0,h); extr(b12,B,n,0,h,h);
    extr(b21,B,n,h,0,h); extr(b22,B,n,h,h,h);

    #pragma omp parallel sections if(d<3) \
        shared(m1,m2,m3,m4,m5,m6,m7,a11,a12,a21,a22,b11,b12,b21,b22) \
        firstprivate(h,d,hyb) default(none)
    {
    #pragma omp section
    { double*u=mnew(h),*v=mnew(h); madd(u,a11,a22,h); madd(v,b11,b22,h);
      strassen(u,v,m1,h,d+1,hyb); free(u);free(v); }
    #pragma omp section
    { double*u=mnew(h); madd(u,a21,a22,h); strassen(u,b11,m2,h,d+1,hyb); free(u); }
    #pragma omp section
    { double*u=mnew(h); msub(u,b12,b22,h); strassen(a11,u,m3,h,d+1,hyb); free(u); }
    #pragma omp section
    { double*u=mnew(h); msub(u,b21,b11,h); strassen(a22,u,m4,h,d+1,hyb); free(u); }
    #pragma omp section
    { double*u=mnew(h); madd(u,a11,a12,h); strassen(u,b22,m5,h,d+1,hyb); free(u); }
    #pragma omp section
    { double*u=mnew(h),*v=mnew(h); msub(u,a21,a11,h); madd(v,b11,b12,h);
      strassen(u,v,m6,h,d+1,hyb); free(u);free(v); }
    #pragma omp section
    { double*u=mnew(h),*v=mnew(h); msub(u,a12,a22,h); madd(v,b21,b22,h);
      strassen(u,v,m7,h,d+1,hyb); free(u);free(v); }
    }
    // combina: C11=M1+M4-M5+M7, C12=M3+M5, C21=M2+M4, C22=M1-M2+M3+M6
    madd(t1,m1,m4,h); msub(t2,m7,m5,h); madd(t1,t1,t2,h);
    madd(t2,m3,m5,h);
    mzero(C,n); insr(C,t1,n,0,0,h); insr(C,t2,n,0,h,h);
    madd(t1,m2,m4,h);
    msub(t2,m1,m2,h); madd(t2,t2,m3,h); madd(t2,t2,m6,h);
    insr(C,t1,n,h,0,h); insr(C,t2,n,h,h,h);

    free(a11);free(a12);free(a21);free(a22);
    free(b11);free(b12);free(b21);free(b22);
    free(m1);free(m2);free(m3);free(m4);free(m5);free(m6);free(m7);
    free(t1);free(t2);
}

static void str_par(const double*A,const double*B,double*C,int n,int d){ strassen(A,B,C,n,d,0); }
static void hyb_par(const double*A,const double*B,double*C,int n,int d){ strassen(A,B,C,n,d,1); }

// Benchmark: mide T, Speedup, Eficiencia
static void bench(const char*name,
                  void(*fn)(const double*,const double*,double*,int,int),
                  const double*A,const double*B,double*C,int n){
    static const int ps[]={1,2,4,8,16,32};
    double t1=-1;
    printf("\n%s (n=%d)\n%-4s  %-9s  %-8s  Eficiencia\n",name,n,"p","T(s)","Speedup");
    for(int i=0;i<6;i++){
        omp_set_num_threads(ps[i]);
        double t0=omp_get_wtime(); fn(A,B,C,n,0); double tp=omp_get_wtime()-t0;
        if(t1<0) t1=tp;
        printf("p=%-2d  %.4fs    %.3f     %.3f\n",ps[i],tp,t1/tp,(t1/tp)/ps[i]);
    }
}

int main(int argc,char**argv){
    int n=(argc>1)?atoi(argv[1]):512, n2=1;
    while(n2<n) n2<<=1;
    n=n2;
    srand(42);
    double*A=mnew(n),*B=mnew(n),*C=mnew(n);
    mrand(A,n); mrand(B,n);
    printf("BS=%d  threshold=%d\n",BS,BS);
    bench("3.1 Bloques paralelo",  blk_par, A,B,C,n);
    bench("3.2 Strassen paralelo", str_par, A,B,C,n);
    bench("3.3 Hibrido",           hyb_par, A,B,C,n);
    free(A);free(B);free(C);
}