/* Reaction-diffusion equation in 1D
 * Solution by Jacobi iteration
 * simple MPI implementation
 *
 * C Michael Hanke 2006-12-12
 */
#define MIN(a,b) ((a) < (b) ? (a) : (b))

/* Use MPI */
#define M_PI 3.14159265358979323846
#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

/* define problem to be solved */
#define N 1000 /* number of inner grid points */
#define SMX 100000 /* number of iterations */

/* implement coefficient functions */
extern double r(const double x);
extern double f(const double x);

/* We assume linear data distribution. The formulae according to the lecture
   are:
      L = N/P;
      R = N%P;
      I = (N+P-p-1)/P;    (number of local elements)
      n = p*L+MIN(p,R)+i; (global index for given (p,i)
   Attention: We use a small trick for introducing the boundary conditions:
      - The first ghost point on p = 0 holds u(0);
      - the last ghost point on p = P-1 holds u(1).
   Hence, all local vectors hold I elements while u has I+2 elements.
*/

int main(int argc, char *argv[])
{
/* local variable */
FILE *fp, *fp1;
int pid, np;
int tag = 0;
double x;
int signal = 1;

/* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Status status;
    int L = N/np;
    int R = N%np;

    if (N < np) {
	fprintf(stdout, "Too few discretization points...\n");
	exit(1);
    }
/* Compute local indices for data distribution */
int I = (N+np-pid-1)/np;
double h = pow(1.0/(double) (N+1), 2);

/* arrays */
    double *unew = (double *) malloc(I*sizeof(double));
    double *fapprox = (double *) malloc(I*sizeof(double));
    double *err = (double *) malloc(I*sizeof(double));
/* Note: The following allocation includes additionally:
   - boundary conditins are set to zero
   - the initial guess is set to zero */
    double *u = (double *) calloc(I+2, sizeof(double));
    for(int i=0; i<I+2;i++){
        u[i] = 0.0;
    }

/* Jacobi iteration */
    for (int step = 0; step < SMX; step++) {
/* RB communication of overlap */
        if(pid%2==0){
            MPI_Send(&u[I], 1, MPI_DOUBLE, pid+1, tag, MPI_COMM_WORLD);
            MPI_Recv(&u[I+1], 1, MPI_DOUBLE, pid+1, tag, MPI_COMM_WORLD, &status);
            if(pid!=0){
                MPI_Send(&u[1], 1, MPI_DOUBLE, pid-1, tag, MPI_COMM_WORLD);
                MPI_Recv(&u[0], 1, MPI_DOUBLE, pid-1, tag, MPI_COMM_WORLD, &status);
            }
        }else{
            MPI_Recv(&u[0], 1, MPI_DOUBLE, pid-1, tag, MPI_COMM_WORLD, &status);
            MPI_Send(&u[1], 1, MPI_DOUBLE, pid-1, tag, MPI_COMM_WORLD);
            if(pid!=np-1){
                MPI_Recv(&u[I+1], 1, MPI_DOUBLE, pid+1, tag, MPI_COMM_WORLD, &status);
                MPI_Send(&u[I], 1, MPI_DOUBLE, pid+1, tag, MPI_COMM_WORLD);
            }
        }
        for(int i=0;i<I;i++){
            x = (double) (pid*L+MIN(pid,R)+i)/(N+1);
            unew[i] = (double)(u[i] + u[i+2] - h*f(x)) / (2.0-h*r(x));
        }
        for(int i=1;i<I+1;i++){
            u[i] = unew[i-1];
        }
/* local iteration step */
/*
	    unew[i] = (u[i]+u[i+2]-h*h*ff[i])/(2.0-h*h*rr[i]);
*/
    }
    if(pid%2==0){
            MPI_Send(&u[I], 1, MPI_DOUBLE, pid+1, tag, MPI_COMM_WORLD);
            MPI_Recv(&u[I+1], 1, MPI_DOUBLE, pid+1, tag, MPI_COMM_WORLD, &status);
            if(pid!=0){
                MPI_Send(&u[1], 1, MPI_DOUBLE, pid-1, tag, MPI_COMM_WORLD);
                MPI_Recv(&u[0], 1, MPI_DOUBLE, pid-1, tag, MPI_COMM_WORLD, &status);
            }
        }else{
            MPI_Recv(&u[0], 1, MPI_DOUBLE, pid-1, tag, MPI_COMM_WORLD, &status);
            MPI_Send(&u[1], 1, MPI_DOUBLE, pid-1, tag, MPI_COMM_WORLD);
            if(pid!=np-1){
                MPI_Recv(&u[I+1], 1, MPI_DOUBLE, pid+1, tag, MPI_COMM_WORLD, &status);
                MPI_Send(&u[I], 1, MPI_DOUBLE, pid+1, tag, MPI_COMM_WORLD);
            }
        }
    for(int i=0;i<I;i++){
        // printf("%f ",u[i]);
        x = (double) (pid*L+MIN(pid,R)+i+1)/(N+1);        
        fapprox[i] = (u[i]-2*u[i+1]+u[i+2])/h + r(x)*u[i+1];      
        err[i] = fapprox[i] - f(x);
    }
/* output for graphical representation */
/* Instead of using gather (which may lead to excessive memory requirements
   on the master process) each process will write its own data portion. This
   introduces a sequentialization: the hard disk can only write (efficiently)
   sequentially. Therefore, we use the following strategy:
   1. The master process writes its portion. (file creation)
   2. The master sends a signal to process 1 to start writing.
   3. Process p waites for the signal from process p-1 to arrive.
   4. Process p writes its portion to disk. (append to file)
   5. process p sends the signal to process p+1 (if it exists).
*/
if(pid==0){
    fp = fopen("hw2.txt","wb");
    fp1 = fopen("hw2_error.txt","wb");
    for(int i=0;i<I;i++){
        fprintf(fp, "%lf ",fapprox[i]);
        fprintf(fp1, "%lf ",err[i]);
    }
    fclose(fp);
    fclose(fp1);
    MPI_Send(&signal, 1, MPI_INT, pid+1, tag, MPI_COMM_WORLD);
}else{
    MPI_Recv(&signal, 1, MPI_INT, pid-1, tag, MPI_COMM_WORLD, &status);
    fp = fopen("hw2.txt","ab");
    fp1 = fopen("hw2_error.txt","ab");
    for(int i=0;i<I;i++){
        fprintf(fp, "%lf ",fapprox[i]);
        fprintf(fp1, "%lf ",err[i]);
    }
    fclose(fp);
    fclose(fp1);
    if(pid!=np-1){
        MPI_Send(&signal, 1, MPI_INT, pid+1, tag, MPI_COMM_WORLD);
    }
}



/* That's it */
    MPI_Finalize();
    exit(0);
}

double r(const double x){
    return -x;
}

double f(const double x){
    // double tmp = 2.0*M_PI*x;
    // return sin(tmp);
    return x*(x-1);
}