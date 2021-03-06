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
#include <time.h>
//#include <fti.h>


#define PRECISION   0.001
#define ITER_TIMES  100000
//#define ITER_OUT    1000
#define WORKTAG     50
#define REDUCE      1
typedef union ldouble
{
    double value;
    unsigned long lvalue;
    unsigned char byte[8];
} ldouble;
typedef union lfloat
{
    float value;
    unsigned int ivalue;
    unsigned char byte[4];
} lfloat;
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
void writeFloatData_inBytes(float *data, size_t nbEle, char* tgtFilePath, int *status)
{
    size_t i = 0; 
    int state = 0;
    lfloat buf;
    unsigned char* bytes = (unsigned char*)malloc(nbEle*sizeof(float));
    for(i=0;i<nbEle;i++)
    {
        buf.value = data[i];
        bytes[i*4+0] = buf.byte[0];
        bytes[i*4+1] = buf.byte[1];
        bytes[i*4+2] = buf.byte[2];
        bytes[i*4+3] = buf.byte[3];                 
    }

    size_t byteLength = nbEle*sizeof(float);
    writeByteData(bytes, byteLength, tgtFilePath, &state);
    free(bytes);
    *status = state;
}
void initData(int nbProcs,int nbLines, int M, int rank, float initTemp,int random,float *h) {
    
    int i, j;
    for (i = 0; i < nbLines; i++) {
        for (j = 0; j < M; j++) {
            if((i==0&&rank==0)||(i==nbLines-1&&rank==nbProcs-1)||j==M-1)
                h[(i*M)+j] = 0.0;
            else if(j==0)
                h[(i*M)+j] = initTemp;
            else{
                if(random){
                    struct timeval time;
                    gettimeofday(&time,NULL);
                    uint32_t seed=time.tv_sec*1000000+time.tv_usec;
                    srand(seed);
                    h[(i*M)+j]=( ((float)(rand()%10000))/10000.0 ) *initTemp;
                }

                else
                    h[(i*M)+j] = 25.0;
            }
            

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

void print_solution (char *filename, float *u, int size)
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

float doWork(int numprocs, int rank, int nbLines, int M, float *g,
 float *h) {
    int i, j;
    MPI_Request req1[2], req2[2];
    MPI_Status status1[2], status2[2];
    float localerror;
    localerror = 0;
    
    int total_lines=nbLines-2;
    //int start=rank*total_lines/numprocs+1;
    //printf("start:%d\n,rank:%d\n",start,rank);
    //int end=(rank+1)*total_lines/numprocs;
    //printf("end:%d\n,rank:%d\n",end,rank);
    for (i = 0; i < nbLines; i++) {
        for (j = 0; j < M; j++) {
            h[(i*M)+j] = g[(i*M)+j];
        }
    }
    
    if (rank > 0) {
        MPI_Isend(g+M, M, MPI_FLOAT, rank-1, WORKTAG,
         MPI_COMM_WORLD, &req1[0]);
        MPI_Irecv(h,   M, MPI_FLOAT, rank-1, WORKTAG,
         MPI_COMM_WORLD, &req1[1]);
    }
    if (rank < numprocs-1) {
        MPI_Isend(g+((nbLines-2)*M), M, MPI_FLOAT, rank+1, WORKTAG,
         MPI_COMM_WORLD, &req2[0]);
        MPI_Irecv(h+((nbLines-1)*M), M, MPI_FLOAT, rank+1, WORKTAG,
         MPI_COMM_WORLD, &req2[1]);
    }
    if (rank > 0) {
        MPI_Waitall(2, req1, status1);
    }
    if (rank < numprocs-1) {
        MPI_Waitall(2, req2, status2);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    for (i = 1; i <nbLines-1; i++) {
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
    int rank, nbProcs, N, i, M, save_start,save_end,save_interval,random;//save_interval optional, -1 means only save last
    char *outfolder;
    float initTemp,wtime, *h, *g, memSize, localerror, globalerror = 1;
    initTemp=atof(argv[1]);
    N = atoi(argv[2]);
    M = atoi(argv[3]);
    outfolder=argv[4];
    
    
    if (argc>=6){
        save_start=atoi(argv[5]);
        save_end=atoi(argv[6]);
    }
    else
        save_end=-1;
    if (argc>=8)
        save_interval=atoi(argv[7]);
    else
        save_interval=1;
    if (argc>=9)
        random=atoi(argv[8]);
    else
        random=0;
    MPI_Init(&argc, &argv);
    /*
    if (FTI_Init(argv[2], MPI_COMM_WORLD) !=  0) {
        printf("FTI could not initialize properly!\n");
        return 1;
    };
    */
    MPI_Comm_size(MPI_COMM_WORLD, &nbProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    

    int start=rank*N/nbProcs;
    
    int end=(rank+1)*N/nbProcs-1;
    int nbLines=end-start+1;
    
    if (rank>0){
        nbLines++;

    }

    if (rank<nbProcs-1){
        nbLines++;
        
    }


   
    
    
    
   
    h = (float *) malloc(sizeof(float) * nbLines * M);
    g = (float *) malloc(sizeof(float) * nbLines * M);
    initData(nbProcs,nbLines,M,rank, initTemp,random,g);
    //initData(nbProcs,nbLines,M,rank, h);
  
    
    memSize = N * M * 3 * sizeof(float) / (float)(1024 * 1024);
    float * result;
    if (rank == 0) {
        printf("Local data size is %d x %d = %f MB.\n", N,
         M, memSize);
        printf("Target precision : %f \n", PRECISION);
        printf("Maximum number of iterations : %d \n", ITER_TIMES);
        result= (float *) malloc(sizeof(float) * N * M);
    }
    /*
    FTI_Protect(0, &i, 1, FTI_INTG);
    FTI_Protect(1, h, M*nbLines, FTI_DBLE);
    FTI_Protect(2, g, M*nbLines, FTI_DBLE);
    */
    
    //MPI_Barrier(MPI_COMM_WORLD);
    
    wtime = MPI_Wtime();
    
    for (i = 0; i < ITER_TIMES; i++) {
        //int checkpointed = FTI_Snapshot();
        localerror = doWork(nbProcs, rank, nbLines, M, g, h);
        //printf("%f\n",g[900]);
        if ( i<=save_end&&i>=save_start && i%save_interval==0 ) {
            //MPI_Request sreq,rreq[100];
            MPI_Status st;
            
            if(rank>0){
                int linesnum=nbLines-2;
                if (rank==nbProcs-1)
                    linesnum++;
                MPI_Send(g+M,linesnum*M,MPI_FLOAT, 0, WORKTAG,
         MPI_COMM_WORLD);
            }
            if(rank==0){
                int ii;
                for(ii=0;ii<=end;ii++){
                    int j;
                    for(j=0;j<M;j++){
                        result[ii*M+j]=g[ii*M+j];
                    }
                }
                int pid;
                for(pid=1;pid<nbProcs;pid++){
                    int pid_start=pid*N/nbProcs;
                    int pid_end=(pid+1)*N/nbProcs-1;
                    MPI_Recv(result+pid_start*M, (pid_end-pid_start+1)*M, MPI_FLOAT, pid, WORKTAG,
         MPI_COMM_WORLD,&st);

                }
            }
            MPI_Barrier(MPI_COMM_WORLD);

            if(rank==0){
                printf("Step : %d, error = %f\n", i, globalerror);
                char filename[100];
                sprintf(filename,"%s/%d.dat",outfolder,i);

                    
    


                int status=-1;
                writeFloatData_inBytes(result, N*M, filename, &status);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if ((i%REDUCE) == 0) {
            MPI_Allreduce(&localerror, &globalerror, 1, MPI_FLOAT, MPI_MAX,
             MPI_COMM_WORLD);
        }
       
        if (globalerror < PRECISION) {
            break;
        }
       
        
        
    }
    
    if (rank == 0) {
        printf("Execution finished in %lf seconds with %d iterations.\n", MPI_Wtime() - wtime,i);
    }
    
    
    //FTI_Finalize();
    //MPI_Request sreq,rreq[100];
    MPI_Status st;
    
            
    if(rank>0){
        int linesnum=nbLines-2;
        if (rank==nbProcs-1)
            linesnum++;
        //printf("baba%d\n",rank);
        MPI_Send(g+M,linesnum*M,MPI_FLOAT, 0, WORKTAG,
    MPI_COMM_WORLD);
        //printf("send from %d, %f\n",rank,g[M+1]);
        
    }
    if(rank==0){
        int ii;
        for(ii=0;ii<=end;ii++){
            int j;
            for(j=0;j<M;j++){
                result[ii*M+j]=g[ii*M+j];
            }
        }
        //printf("%f\n",g[M]);
        //printf("%f\n",result[M]);
        int pid;
        for(pid=1;pid<nbProcs;pid++){
            int pid_start=pid*N/nbProcs;
            int pid_end=(pid+1)*N/nbProcs-1;
            //printf("mama%d\n",pid_start);
            MPI_Recv(result+pid_start*M, (pid_end-pid_start+1)*M, MPI_FLOAT, pid, WORKTAG,MPI_COMM_WORLD,&st);
            //printf("received from%d %f\n",pid,result[pid_start*M+1]);
            

        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(rank==0){
        char filename[100];
        sprintf(filename,"%s/%d.dat",outfolder,i);
        int status=-1;
        //printf("%f\n",result[M]);
        //printf("%f\n",result[M+M/2]);
        //printf("%f\n",result[(N/2)*M+M/2]);
        writeFloatData_inBytes(result, N*M, filename, &status);
       
       // free(result);
    }
    //free(g);
    //free(h);
    MPI_Finalize();
    
    return 0;
}