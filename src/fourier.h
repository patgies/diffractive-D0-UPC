// Last modified on August 16th, 2017 -- F. Gelis


#ifndef FORTRAN
typedef double (myfunc)(double,void *);
#else
typedef double (myfunc)(double*,void *);
#endif


extern int Nzeros;
extern double *workspace;

extern "C"
{

double fourier_j0(double x,myfunc *func,void* param);
double fourier_j0_i(double x,int Ni,myfunc *func,void *param);
double fourier_j0_f(double x,int Nf,myfunc *func,void *param);
double fourier_j0_if(double x,int Ni,int Nf,myfunc *func,void *param);
double fourier_j1(double x,myfunc *func,void* param);
double fourier_j1_i(double x,int Ni,myfunc *func,void *param);
double fourier_j1_f(double x,int Nf,myfunc *func,void *param);
double fourier_j1_if(double x,int Ni,int Nf,myfunc *func,void *param);

void init_workspace_fourier(int N);
void set_fpu_state(void);
void set_fourier_precision(double e,double e1);
}
