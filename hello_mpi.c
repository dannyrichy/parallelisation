#include<stdio.h>
#include<complex.h>
#include<stdlib.h>
#include<mpi.h>

unsigned char cal_pixel(double complex d, double b, unsigned char N){
    unsigned char count = 1;
    double complex z = 0;
    while((cabs(z) < b) && (count < N)){
        z = z*z + d;
        count++;
    }
    return count;
}

int main(int argc, char** argv) {
    FILE *fp;
    unsigned char N = 255;
    double b = 2;
    int w = 2048;
    int h = 2048;
    double dx = 2*b/(w-1);
    double dy = 2*b/(h-1);
    double dreal, dimag;
    double complex d;

    int pid, np, wp, n_elements_recieved;
    MPI_Status status;
    unsigned char *color, *color_sub;


    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);
    
    wp = w / np;
    

    if(pid==0){
        color = malloc((w*h)*sizeof(unsigned char));
        color_sub = malloc((wp*h)*sizeof(unsigned char));
    }
    else{
        color = malloc((wp*h)*sizeof(unsigned char));
    }
    int h_iter;
    for(int w_iter=0;w_iter<wp;w_iter++){
        for(h_iter=0;h_iter<h;h_iter++){
            color[w_iter*h+h_iter] = cal_pixel(
                (wp*(pid)+w_iter)*dx - b + (h_iter*dy - b)*I, b, N 
            );
        }
    }
    printf("Calculated values in process %d\n", pid);
    if(pid!=0){
        printf("Sending values from process %d\n", pid);
        MPI_Send(color, wp*h, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
        // free(color);    
    }
    else{
        for(int pro_id=1; pro_id<np;pro_id++){
            printf("receiving values from process %d\n", pro_id);
            MPI_Recv(color_sub, wp*h, MPI_UNSIGNED_CHAR, pro_id, 0, MPI_COMM_WORLD, &status);
            printf("PID:%d -> Starting index: %d Ending index: %d\n",pro_id, (wp*h)*pro_id, (wp*h)*pro_id+wp*h-1);
            for(int it=0; it<wp*h;it++){
                color[(wp*h)*pro_id+it] = color_sub[it];
            }
            printf("\nreceived values from process %d\n", pro_id);
        }
    }


    if(pid==0){
        fp = fopen("color.txt","wb");
        // printf("opened file");
        for(int i = 0; i < h; i++) {
            for (int j = 0; j < w; j++) {
                fprintf(fp, "%hhu ", color[i*w+j]);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }
        
    
    MPI_Finalize();
    
}