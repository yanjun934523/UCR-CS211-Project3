/*
 *   Sieve of Eratosthenes
 *
 *   Programmed by Michael J. Quinn
 *
 *   Last modification: 7 September 2001
 */

#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a,b)  ((a)<(b)?(a):(b))

int main (int argc, char *argv[])
{
   unsigned long int    count;        /* Local prime count */
   double elapsed_time; /* Parallel execution time */
   unsigned long int    first;        /* Index of first multiple */
   int   local_first;
   unsigned long int    global_count = 0; /* Global prime count */
   unsigned long long int    high_value;   /* Highest value on this proc */
   unsigned long int    i;
   int    id;           /* Process ID number */
   unsigned long int    index;        /* Index of current prime */
   unsigned long long int    low_value;    /* Lowest value on this proc */
   char  *marked;       /* Portion of 2,...,'n' */
   char  *local_prime_marked;
   unsigned long long int    n;            /* Sieving from 2, ..., 'n' */
   int    p;            /* Number of processes */
   unsigned long int    proc0_size;   /* Size of proc 0's subarray */
   unsigned long int    prime;
   unsigned long int  local_prime;        /* Current prime */
   unsigned long int    size;         /* Elements in 'marked' */
   unsigned long int  local_prime_size;


   MPI_Init (&argc, &argv);

   /* Start the timer */

   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   MPI_Comm_size (MPI_COMM_WORLD, &p);
   MPI_Barrier(MPI_COMM_WORLD);
   elapsed_time = -MPI_Wtime();

   if (argc != 2) {
      if (!id) printf ("Command line: %s <m>\n", argv[0]);
      MPI_Finalize();
      exit (1);
   }

   n = atoll(argv[1]);

   /* Figure out this process's share of the array, as
      well as the integers represented by the first and
      last array elements */



   /* Add you code here  */
   long long int l_size;
   long long int l_first;
   char* l_marked;
   low_value = 3 + 2 * (long long int)(id * ((n - 1) / 2) / p);
   high_value = 1 + 2 * (long long int)((id + 1) * ((n - 1) / 2) / p);
   size = (high_value - low_value) / 2 + 1;
   proc0_size = (n - 1) / 2 / p;
   l_size = sqrt(n);

   if ((3 + 2 * proc0_size) < (long long int) sqrt((double) n)) {
      if (!id) printf("Too many processes\n");
      MPI_Finalize();
      exit(1);
   }

   marked = (char *) malloc(size);
   l_marked = (char *) malloc(l_size);

   if (marked == NULL || l_marked == NULL) {
      printf("Cannot allocate enough memory\n");
      MPI_Finalize();
      exit(1);
   }

   for (i = 0; i < size; i++) marked[i] = 0;
   for (i = 0; i < l_size; i++) l_marked[i] = 0;

   unsigned long int B = 100000;
   unsigned long int j;
   for (i = 0; i < size; i+= B) {
      index = 0;
      prime = 3; 
      do {
            if (prime * prime > low_value) {
               first = (prime * prime - low_value) / 2;
            } else {
               if (!(low_value % prime)) {
                  first = 0;
               } else {
                  if ((prime - (low_value % prime)) % 2 == 0) {
                     first = (prime - (low_value % prime)) / 2;
                  } else {
                     first = (prime * 2 - (low_value % prime)) / 2;
                  }
               }
            }
            l_first = (prime * prime - 3) / 2;
            for (j = i + first; j < size && j < i + B; j += prime) marked[j] = 1;
            for (j = l_first; j < l_size; j+= prime) l_marked[j] = 1; 
            while (l_marked[++index]);
            prime = index * 2 + 3;
       } while (prime * prime <= n);
           low_value += B * 2;
   }
   count = 0;
   global_count = 0;
   for (i = 0; i < size; i++)
      if (!marked[i]) count++;
   if (p > 1) {
      MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM,
                  0, MPI_COMM_WORLD);
   }

   /* Stop the timer */

   elapsed_time += MPI_Wtime();

 













   /* Print the results */

   if (!id) {
      printf("The total number of prime: %ld, total time: %10.6f, total node %d\n", global_count, elapsed_time, p);

   }
   MPI_Finalize ();
   return 0;
}

