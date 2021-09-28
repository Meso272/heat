/**
 *  Copyright (c) 2017 Leonardo A. Bautista-Gomez
 *  All rights reserved
 *
 *  @file   heatdis.c
 *  @author Leonardo A. Bautista Gomez
 *  @date   May, 2014
 *  @brief  Heat distribution code to test FTI.
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
//#include <fti.h>


#define PRECISION   0.001
#define ITER_TIMES  10000
//#define ITER_OUT    1000
#define WORKTAG     50
#define REDUCE      1
typedef union ldouble
{
    double value;
    unsigned long lvalue;
    unsigned char byte[8];
} ldouble;

void writeByteData(unsigned char *bytes, size_t byteLength, char *tgtFilePath, int *status)
{
    FILE *pFile = fopen(tgtFilePath, "wb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 3\n");
        *status = 1;
        return;
    }
    
    fwrite(bytes, 1, byteLength, pFile); //write outSize bytes
    fclose(pFile);
    *status = 0;
}

void writeDoubleData_inBytes(double *data, size_t nbEle, char* tgtFilePath, int *status)
{
    size_t i = 0, index = 0; 
    int state = 0;
    ldouble buf;
    unsigned char* bytes = (unsigned char*)malloc(nbEle*sizeof(double));
    for(i=0;i<nbEle;i++)
    {
        index = i*8;
        buf.value = data[i];
        bytes[index+0] = buf.byte[0];
        bytes[index+1] = buf.byte[1];
        bytes[index+2] = buf.byte[2];
        bytes[index+3] = buf.byte[3];
        bytes[index+4] = buf.byte[4];
        bytes[index+5] = buf.byte[5];
        bytes[index+6] = buf.byte[6];
        bytes[index+7] = buf.byte[7];
    }

    size_t byteLength = nbEle*sizeof(double);
    writeByteData(bytes, byteLength, tgtFilePath, &state);
    free(bytes);
    *status = state;
}

void initData(int nbLines, int M, int rank, double *h) {
    if (rank!=0)
        return;
    int i, j;
    for (i = 0; i < nbLines; i++) {
        for (j = 0; j < M; j++) {
            h[(i*M)+j] = 25.0;
            if(j==0)
                h[(i*M)+j] = 1000.0;
            if(i==0||i==nbLines-1||j==M-1)
                h[(i*M)+j] = 0.0;

        }

    }
    /*
    if (rank == 0) {
        for (j = (M*0.1); j < (M*0.9); j++) {
            h[j] = 100;
        }
    }
    */
}

void print_solution (char *filename, double *u, int size)
{
   int i, j;
   char sep;
   FILE *outfile;
   
   if (!filename) {
      sep = '\t';   /* just for easier view */
      outfile = stdout;
   } else {
      sep = '\n';   /* for gnuplot format */
      outfile = fopen(filename,"w");
      if (outfile == NULL) {
         printf("Can't open output file.");
         exit(-1);
      }
   }

   /* Print the solution array */
   for (i = 0; i < size; i++) {
         fprintf (outfile, "%6.10f%c", u[i], sep);
      fprintf(outfile, "\n"); /* Empty line for gnuplot */
   }
   if (outfile != stdout)
      fclose(outfile);


   
}

double doWork(int numprocs, int rank, int nbLines, int M, double *g,
 double *h) {
    int i, j;
    MPI_Request req1[2], req2[2];
    MPI_Status status1[2], status2[2];
    double localerror;
    localerror = 0;
    
    int total_lines=nbLines-2;
    int start=rank*total_lines/numprocs+1;
    int end=(rank+1)*total_lines/numprocs;
    for (i = start; i <= end; i++) {
        for (j = 1; j < M-1; j++) {
            h[(i*M)+j] = g[(i*M)+j];
        }
    }
    /*
    if (rank > 0) {
        MPI_Isend(g+M, M, MPI_DOUBLE, rank-1, WORKTAG,
         MPI_COMM_WORLD, &req1[0]);
        MPI_Irecv(h,   M, MPI_DOUBLE, rank-1, WORKTAG,
         MPI_COMM_WORLD, &req1[1]);
    }
    if (rank < numprocs-1) {
        MPI_Isend(g+((nbLines-2)*M), M, MPI_DOUBLE, rank+1, WORKTAG,
         MPI_COMM_WORLD, &req2[0]);
        MPI_Irecv(h+((nbLines-1)*M), M, MPI_DOUBLE, rank+1, WORKTAG,
         MPI_COMM_WORLD, &req2[1]);
    }
    if (rank > 0) {
        MPI_Waitall(2, req1, status1);
    }
    if (rank < numprocs-1) {
        MPI_Waitall(2, req2, status2);
    }
    */
    MPI_Barrier(MPI_COMM_WORLD);
    for (i = start; i <=end; i++) {
        for (j = 1; j < M-1; j++) {
            g[(i*M)+j] = 0.25*(h[((i-1)*M)+j]+h[((i+1)*M)+j]+
                h[(i*M)+j-1]+h[(i*M)+j+1]);
            if (localerror < fabs(g[(i*M)+j] - h[(i*M)+j])) {
                localerror = fabs(g[(i*M)+j] - h[(i*M)+j]);
            }
        }
    }
    /*
    if (rank == (numprocs-1)) {
        for (j = 0; j < M; j++) {
            g[((nbLines-1)*M)+j] = g[((nbLines-2)*M)+j];
        }
    }
    */
    return localerror;
}


int main(int argc, char *argv[]) {
    int rank, nbProcs, N, i, M, save_interval;//save_interval optional, -1 means only save last
    char *outfolder;
    double wtime, *h, *g, memSize, localerror, globalerror = 1;
    N = atoi(argv[1]);
    M = atoi(argv[2]);
    outfolder=argv[3];
    if (argc>=5)
        save_interval=atoi(argv[4]);
    else
        save_interval=-1;
    MPI_Init(&argc, &argv);
    /*
    if (FTI_Init(argv[2], MPI_COMM_WORLD) !=  0) {
        printf("FTI could not initialize properly!\n");
        return 1;
    };
    */
    MPI_Comm_size(MPI_COMM_WORLD, &nbProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    
    
    
    h = (double *) malloc(sizeof(double) * N * M);
    g = (double *) malloc(sizeof(double) * N *M);
    initData(N,M,rank, g);
    initData(N,M,rank, h);
    memSize = N * M * 2 * sizeof(double) / (double)(1024 * 1024);

    if (rank == 0) {
        printf("Local data size is %d x %d = %f MB.\n", N,
         M, memSize);
        printf("Target precision : %f \n", PRECISION);
        printf("Maximum number of iterations : %d \n", ITER_TIMES);
    }
    /*
    FTI_Protect(0, &i, 1, FTI_INTG);
    FTI_Protect(1, h, M*nbLines, FTI_DBLE);
    FTI_Protect(2, g, M*nbLines, FTI_DBLE);
    */
    wtime = MPI_Wtime();
    for (i = 0; i < ITER_TIMES; i++) {
        //int checkpointed = FTI_Snapshot();
        localerror = doWork(nbProcs, rank, N, M, g, h);
        
        if ( (save_interval>0)&&((i%save_interval) == 0) && (rank == 0)) {
            printf("Step : %d, error = %f\n", i, globalerror);
            char filename[100];
            sprintf(filename,"%s/%d.dat",outfolder,i);
            int status=-1;
            writeDoubleData_inBytes(g, N*M, filename, &status);
        }
        if ((i%REDUCE) == 0) {
            MPI_Allreduce(&localerror, &globalerror, 1, MPI_DOUBLE, MPI_MAX,
             MPI_COMM_WORLD);
        }
        if (globalerror < PRECISION) {
            break;
        }
    }
    if (rank == 0) {
        printf("Execution finished in %lf seconds.\n", MPI_Wtime() - wtime);
    }
    

    //FTI_Finalize();
    if(rank==0){
        char filename[100];
        sprintf(filename,"%s/%d.dat",outfolder,i);
        int status=-1;
        printf("%f\n",g[200]);
        writeDoubleData_inBytes(g, N*M, filename, &status);
        free(h);
        free(g);
    }
    MPI_Finalize();
    
    return 0;
}