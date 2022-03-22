#include<stdio.h>
#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>

int comp(const void *a, const void *b){
    double *x = (double *) a;
    double *y = (double *) b;
    if(*x<*y) return -1;
    else if(*x>*y) return 1;
    return 0;
}

double * merge_first_half(double* a, double* x, int I){
    
    double *tmp_merge = (double *) malloc(I*sizeof(double));
    int i=0, j=0, k=0;
    while(i<I && j<I && k<I){
        if(a[i]<x[j]){
            tmp_merge[k++] = a[i];
            i++;
        }
        else{
            tmp_merge[k++] = x[j];
            j++;
        }
    }
    // if(i>I && k<I) while(j<I) tmp_merge[k++] = x[j++];
    // if(j>I && k<I) while(i<I) tmp_merge[k++] = a[i++];  
    return tmp_merge; 
}

double * merge_second_half(double* a, double* x, int I ){
    double *tmp_merge = (double *) malloc(I*sizeof(double));
    int i=I-1, j=I-1, k=I-1;
    while(i>=0 && j>=0 && k>=0){
        if(a[i]>x[j]){
            tmp_merge[k--] = a[i];
            i--;
        }
        else{
            tmp_merge[k--] = x[j];
            j--;
        }
    }
    // if(i>I && k<I) while(j<I) tmp_merge[k++] = x[j++];
    // if(j>I && k<I) while(i<I) tmp_merge[k++] = a[i++];  
    return tmp_merge; 
}

void error_exit(char* nm){
    printf("%s", nm);
    exit(1);
}

int main(int argc, char **argv)
{  
  int P, myrank, N, I;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &P);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Status status;

  int even_process = myrank%2;
  int even_phase = 1;
  int tag=0;
     
  /* Find problem size N from command line */
  if (argc < 2) error_exit("No size N given");
  N = atoi(argv[1]);
  /* local size. Modify if P does not divide N */
  I = N/P;
  double *x = (double *) malloc(I*sizeof(double));
  double *tmp = (double *) malloc(I*sizeof(double));
  /* random number generator initialization */
  srandom(myrank+1);
  /* data generation */
  for (int i = 0; i < I; i++)
    x[i] = ((double) random())/(RAND_MAX);


  qsort(x, I, sizeof(double), comp);


  for(int i=0;i<P;i++){
      if((even_phase & even_process) | (!even_phase & !even_process)){
          if(myrank!=P-1){
              printf("Sending it first %d -> %d \n", myrank, myrank+1);
          MPI_Send(x, I, MPI_DOUBLE, myrank+1, tag, MPI_COMM_WORLD);
          MPI_Recv(tmp, I, MPI_DOUBLE, myrank+1, tag, MPI_COMM_WORLD, &status);
          x = merge_first_half(x, tmp, I);
      }
      }
      else{
          if(myrank!=0){
              printf("Recieving it first %d -> %d \n", i, myrank);
            MPI_Recv(tmp, I, MPI_DOUBLE, myrank-1, tag, MPI_COMM_WORLD, &status);
            MPI_Send(x, I, MPI_DOUBLE, myrank-1, tag, MPI_COMM_WORLD);
            x =  merge_second_half(x, tmp, I);
          }
      }
      even_phase = !even_phase;
  }
  printf("%d -> ", myrank);
  for(int p=0;p<I;p++){
      printf("%lf ",x[p]);
  } 
  MPI_Finalize();
  exit(0);
}