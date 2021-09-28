/*
 *   Sequential program that solves the steady-state temperature
 *   distribution problem using the Jacobi method.
 *
 *   Compile on Belief:
 *   export PATH=/opt/SUNWspro/bin:$PATH
 *   cc -xO3 -o jacobi_seq jacobi_seq.c -lm
 *   cc -xO3 -o jacobi_seq jacobi_seq.c -lm -DINTERACTIVE (for interactive output display in gnuplot)
 *
 *   Compile on Gideon:
 *   cc -O3 -o jacobi_seq jacobi_seq.c -lm
 *   cc -O3 -o jacobi_seq jacobi_seq.c -lm -DINTERACTIVE (for interactive output display in gnuplot)
 *
 *   Execute:
 *   ./jacobi_seq [optional: <rows> <cols>]
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <values.h>              /* For timing */
#include <sys/time.h>            /* For timing */

#define MAX(a,b) ((a)>(b)?(a):(b))
#define EPSILON 0.001            /* Termination condition */

#define VIEWABLE_SIZE 20         /* Under this grid size, the array is printed for verification. */
char filename[100];              /* File name of output file */

/* Grid size */
int M = 200;                     /* Number of rows */
int N = 200;                     /* Number of cols */
long max_its = 100000;           /* Maximum iterations */
double final_diff;               /* Temperature difference between iterations at the end */

/* Interactive mode is solely for demonstration purpose on gnuplot */
#ifdef INTERACTIVE
int refresh_its = 1000;          /* Overwrite output file per this batch of iterations */
#endif

int main (int argc, char *argv[])
{
   double** u;                   /* Previous temperatures */
   double** w;                   /* New temperatures */
   int      its;                 /* Iterations to converge */
   double   elapsed;             /* Execution time */
   struct timeval stime, etime;  /* Start and end times */
   
   void allocate_2d_array (int, int, double ***);
   void initialize_array (double ***);
   void print_solution (char *, double **);
   int  find_steady_state (double **, double **);
   
   /* For convenience of other problem size testing */
   if (argc > 1) {
      if (argc != 3) {
         printf("Usage: %s <rows> <cols>\n", argv[0]);
         exit(-1);
      } else {
         M = atoi(argv[1]);
         N = atoi(argv[2]);
      }
   } // Otherwise use default grid size
   printf("Problem size: M=%d, N=%d\n", M, N);
   sprintf(filename, "%s.dat", argv[0]);

   allocate_2d_array (M, N, &u);
   allocate_2d_array (M, N, &w);
   initialize_array (&u);
   initialize_array (&w);

   /* For convenience of correctness checking */
   if (MAX(M, N) <= VIEWABLE_SIZE) {
      printf ("Before:\n");
      print_solution (NULL, u);
   }

   gettimeofday (&stime, NULL);
   its = find_steady_state (u, w);
   gettimeofday (&etime, NULL);
   
   elapsed = ((etime.tv_sec*1000000+etime.tv_usec)-(stime.tv_sec*1000000+stime.tv_usec))/1000000.0;

   /* For convenience of correctness checking */
   if (MAX(M, N) <= VIEWABLE_SIZE) {
      printf ("After:\n");
      print_solution (NULL, w);
   }

   printf("Converged after %d iterations with error: %8.8f.\n", its, final_diff);
   printf("Elapsed time = %8.6f sec.\n", elapsed);
   //print_solution (filename, w);
   double * result=(double *)malloc(M*N*sizeof(double));
   ConvertDoubleArray_2Dto1D(w,result,M,N);
   int status=-1;
   writeDoubleData_inBytes(result, M*N, filename, &status);
   free(result);
   int i;
   /*
   for(i=0;i<M;i++){
       free(u[i]);
       free(w[i]);
   }
   free(u);
   free(w);
   */
   return 0;
}

/* Allocate two-dimensional array. */
void allocate_2d_array (int r, int c, double ***a)
{
   double *storage;
   int     i;
   storage = (double *) malloc (r * c * sizeof(double));
   *a = (double **) malloc (r * sizeof(double *));
   for (i = 0; i < r; i++)
      (*a)[i] = &storage[i * c];
}

/* Set initial and boundary conditions */
void initialize_array (double ***u)
{
   int i, j;
   
   /* Set initial values and boundary conditions */
   for (i = 0; i < M; i++) {
      for (j = 0; j < N; j++)
         (*u)[i][j] = 25.0;      /* Room temperature */
      (*u)[i][0] = 1000.0;       /* Heat source */
      (*u)[i][N-1] = 0.0;
   }
   
   for (j = 0; j < N; j++) {
      (*u)[0][j] = 0.0;
      (*u)[M-1][j] = 0.0;
   }
}

/* Print solution to standard output or a file */
void print_solution (char *filename, double **u)
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
   for (i = 0; i < M; i++) {
      for (j = 0; j < N; j++) 
         fprintf (outfile, "%6.2f%c", u[i][j], sep);
      fprintf(outfile, "\n"); /* Empty line for gnuplot */
   }
   if (outfile != stdout)
      fclose(outfile);
   
#ifdef INTERACTIVE
   /* Copy to another data file for gnuplot to plot */
   if (outfile != stdout) {
      char tempname[100];              /* Temp file name for gnuplot */
      char cmd[1024];                  /* System command buffer */
      sprintf(tempname, "jacobi.dat");
      sprintf(cmd, "cp %s %s", filename, tempname);
      system(cmd);
   }
#endif
}

int find_steady_state (double **u, double **w)
{
   double diff;            /* Maximum temperature difference */
   int    its;             /* Iteration count */
   int    i, j;
#ifdef INTERACTIVE
   int    r = 0;
#endif

   for (its = 0; its < max_its; its++) {
      diff = 0.0;
      for (i = 1; i < M-1; i++) {
         for (j = 1; j < N-1; j++) {
            w[i][j] = 0.25 * (u[i-1][j] + u[i+1][j] + u[i][j-1] + u[i][j+1]);
            if (fabs(w[i][j] - u[i][j]) > diff)
               diff = fabs(w[i][j] - u[i][j]);
         }
      }

      for (i = 1; i < M-1; i++)
         for (j = 1; j < N-1; j++)
            u[i][j] = w[i][j];

      /* Terminate if temperatures have converged */
      if (diff <= EPSILON)
         break;

#ifdef INTERACTIVE
      if (r == refresh_its) {
         print_solution (filename, w);
         r = 0;
      }
      r++;
#endif
   }
   final_diff = diff;
   return its;
}

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

void ConvertDoubleArray_2Dto1D(double **Array_2D,double * Array_1D, size_t M,size_t N){
   int i,j;
   for(i=0;i<M;i++){
      for(j=0;j<N;j++){
         Array_1D[i*M+j]=Array_2D[i][j];
      }
   }
}