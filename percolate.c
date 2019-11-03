#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "percolate.h"

/*
 * Serial program to test for percolation of a cluster.
 */



int main(int argc, char *argv[])
{
  /*
   *  Define the main arrays for the simulation
   */
  MPI_Init(&argc,&argv);
  int old[M+2][N+2], new[M+2][N+2];

  /*
   *  Additional array WITHOUT halos for initialisation and IO. This
   *  is of size LxL because, even in our parallel program, we do
   *  these two steps in serial
   */

  int map[L][L];

  /*
   *  Variables that define the simulation
   */

  int seed;
  double rho;

  /*
   *  Local variables
   */

  int i, j, nhole, step, maxstep, oldval, newval, nchange, printfreq;
  int itop, ibot, perc;
  double r;

  if (argc != 2)
    {
      printf("Usage: percolate <seed>\n");
      return 1;
    }

  /*
   *  Set most important value: the rock density rho (between 0 and 1)
   */

  rho = 0.411;

  /*
   *  Set the randum number seed and initialise the generator
   */

  seed = atoi(argv[1]);

  printf("percolate: params are L = %d, rho = %f, seed = %d\n", L, rho, seed);

  rinit(seed);


  /*MPI init
  *
  *
  */

  
  int rank,size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  int direction,disp;
  int left,right;
  int dims[ndims];
  int period[ndims];
  int reorder;



  int per_process_M=M/size;
  int per_process_area=per_process_M*N;
  int **local_map=alloc_2d_int(per_process_M,N);

  printf("each M=%d,N=%d\n",per_process_M,N);

  MPI_Barrier(MPI_COMM_WORLD);


  //topology
  dims[0]=0;
  period[0]=TRUE;
  reorder=FALSE;
  direction=0;
  disp=1;

  MPI_Comm topo_comm;

  MPI_Dims_create(size,ndims,dims);
  MPI_Cart_create(MPI_COMM_WORLD,ndims,dims,period,reorder,&topo_comm);
  MPI_Comm_rank(topo_comm,&rank);
  MPI_Cart_shift(topo_comm,direction,disp,&left,&right);



  /*
  for(i=0;i<per_process_M;i++){
                for(j=0;j<N;j++){
                        local_map[i][j]=i*j;
                }
        }
  */

  /*
   *  Initialise map with density rho. Zero indicates rock, a positive
   *  value indicates a hole. For the algorithm to work, all the holes
   *  must be initialised with a unique integer
   */

  //begin to intialize
  if(rank==0){
  nhole = 0;

  for (i=0; i < L; i++)
    {
      for (j=0; j < L; j++)
	{
	  r=uni();
	  
	  if(r < rho)
	    {
	      map[i][j] = 0;
	    }
	  else
	    {
	      nhole++;
	      map[i][j] = nhole;
	    }
	}
    }



   /*for(i=0;i<L;i++){
    for(j=0;j<L;j++){
      map[i][j]=i*L+j;
    }
   }*/ 

  printf("percolate: rho = %f, actual density = %f\n",
	 rho, 1.0 - ((double) nhole)/((double) L*L) );


  //display_matrix(map,M,N);

  }


  MPI_Barrier(MPI_COMM_WORLD);
  //scatter map to each prece
  
  /*
  if(rank==0){
    for(i=per_process_M;i<per_process_M*2;i++){
      for(j=0;j<N;j++){
        printf("%d\t",map[i][j]);
      }
      printf("\n");
    }
  }*/
  

  MPI_Scatter(&(map[0][0]),per_process_area,MPI_INT,&(local_map[0][0]),per_process_area,MPI_INT,0,topo_comm); 



/*
   * Initialise the old array: copy the LxL array map to the centre of
   * old, and set the halo values to zero.
   
*/
  //int local_old[per_process_area+2][N+2];
  //int local_new[per_process_M+2][N+2];

  int ** local_old=alloc_2d_int(per_process_M+2,N+2);
  int ** local_new=alloc_2d_int(per_process_M+2,N+2);


  for (i=1; i <= per_process_M; i++)
  {
      for (j=1; j <= N; j++)
    {
      local_old[i][j] = local_map[i-1][j-1];
    }
  }

   for (i=0; i <= per_process_M+1; i++)  // zero the bottom and top halos
  {
      local_old[i][0]   = 0;
      local_old[i][N+1] = 0;
  }

  for (j=0; j <= N+1; j++)  // zero the left and right halos
    {
      local_old[0][j]   = 0;
      local_old[per_process_M+1][j] = 0;
    }


  /*
   *  Update for a fixed number of iterations
   */ 

  maxstep = 16*L;
  printfreq = 100;

  step = 1;
  nchange = 1;

  while (step <= maxstep)
  {
    nchange = 0;

    for (i=1; i<=per_process_M; i++)
     {
       for (j=1; j<=N; j++)
       {
         oldval = local_old[i][j];
         newval = oldval;

        /*
         * Set local[i][j] to be the maximum value of local_old[i][j]
         * and its four nearest neighbours
         */

         if (oldval != 0)
          {
            if (local_old[i-1][j] > newval) newval = local_old[i-1][j];
            if (local_old[i+1][j] > newval) newval = local_old[i+1][j];
            if (local_old[i][j-1] > newval) newval = local_old[i][j-1];
            if (local_old[i][j+1] > newval) newval = local_old[i][j+1];

            if (newval != oldval)
            {
                ++nchange;
            }
          }

        local_new[i][j] = newval;
      }
    }
      /*
       *  Report progress every now and then
       */

    /*
    if (step % printfreq == 0)
    {
    printf("percolate: number of changes on step %d is %d\n",
     step, nchange);
    }
    */

      /*
       *  Copy back in preparation for next step, omitting halos
       */

      for (i=1; i<=per_process_M; i++)
      {
        for (j=1; j<=N; j++)
        {
        local_old[i][j] = local_new[i][j];
        }
      }

        step++;
    }
  


   /*
    *  Update for a fixed number of iterations
    */

if (nchange != 0)
    {
      printf("percolate: WARNING max steps = %d reached before nchange = 0\n",
       maxstep);
    }

  /*
   *  Copy the centre of old, excluding the halos, into map
   */
  
  for (i=1; i<=per_process_M; i++)
    {
      for (j=1; j<=N; j++)
      {
         local_map[i-1][j-1] = local_old[i][j];
      }
    }











  MPI_Barrier(MPI_COMM_WORLD);

  



  //copy back to rank 0
  
  //MPI_Gather(&(local_map[0][0]),per_process_area,MPI_INT,&(map[0][rank*per_process_M]),per_process_area,MPI_INT,0,MPI_COMM_WORLD);
  
  /*if(rank==1){
    display_matrix(local_map,per_process_M,N);
  }*/
  
  MPI_Gather(&(local_map[0][0]),per_process_area,MPI_INT,&(map[0][rank*per_process_M]),per_process_area,MPI_INT,0,topo_comm);
  

  
  if(rank==0){
    for(i=per_process_M;i<per_process_M*2;i++){
      for(j=0;j<N;j++){
        printf("%d\t",map[i][j]);
      }
      printf("\n");
    }
  }
  




  //
  
  //free(local_map[0]);
  //free(local_map);
  //printf("free successfully\n");
  //MPI_Barrier(MPI_COMM_WORLD);


  //display_matrix(local_map,per_process_M,N);
  //printf("%d",local_map[0][0]);



  

  
  /*
   * Initialise the old array: copy the LxL array map to the centre of
   * old, and set the halo values to zero.
   

   for (i=1; i <= M; i++)
    {
      for (j=1; j <= N; j++)
	{
	  old[i][j] = map[i-1][j-1];
	}
    }

   for (i=0; i <= M+1; i++)  // zero the bottom and top halos
    {
      old[i][0]   = 0;
      old[i][N+1] = 0;
    }

   for (j=0; j <= N+1; j++)  // zero the left and right halos
    {
      old[0][j]   = 0;
      old[M+1][j] = 0;
    }

   /*
    *  Update for a fixed number of iterations
    

  maxstep = 16*L;
  printfreq = 100;

  step = 1;
  nchange = 1;

  while (step <= maxstep)
  {
    nchange = 0;

    for (i=1; i<=M; i++)
	   {
	     for (j=1; j<=N; j++)
	     {
	       oldval = old[i][j];
	       newval = oldval;

	      /*
	       * Set new[i][j] to be the maximum value of old[i][j]
	       * and its four nearest neighbours
	       

	       if (oldval != 0)
		      {
		        if (old[i-1][j] > newval) newval = old[i-1][j];
		        if (old[i+1][j] > newval) newval = old[i+1][j];
		        if (old[i][j-1] > newval) newval = old[i][j-1];
		        if (old[i][j+1] > newval) newval = old[i][j+1];

		        if (newval != oldval)
		        {
		            ++nchange;
		        }
		      }

	      new[i][j] = newval;
	    }
	}
      /*
       *  Report progress every now and then
       

      if (step % printfreq == 0)
	{
	  printf("percolate: number of changes on step %d is %d\n",
		 step, nchange);
	}

      /*
       *  Copy back in preparation for next step, omitting halos
       

      for (i=1; i<=M; i++)
	{
	  for (j=1; j<=N; j++)
	    {
	      old[i][j] = new[i][j];
	    }
	}

      step++;
    }



  /*
   *  We set a maximum number of steps to ensure the algorithm always
   *  terminates. However, if we hit this limit before the algorithm
   *  has finished then there must have been a problem (e.g. maxstep
   *  is too small)
   

  if (nchange != 0)
    {
      printf("percolate: WARNING max steps = %d reached before nchange = 0\n",
	     maxstep);
    }

  /*
   *  Copy the centre of old, excluding the halos, into map
   
  
  for (i=1; i<=M; i++)
    {
      for (j=1; j<=N; j++)
	{
	  map[i-1][j-1] = old[i][j];
	}
    }

  /*
   *  Test to see if percolation occurred by looking for positive numbers
   *  that appear on both the top and bottom edges
   

  perc = 0;

  for (itop=0; itop < L; itop++)
    {
      if (map[itop][L-1] > 0)
	{
	  for (ibot=0; ibot < L; ibot++)
	    {
	      if (map[itop][L-1] == map[ibot][0])
		{
		  perc = 1;
		}
	    }
	}
    }

  if (perc != 0)
    {
      printf("percolate: cluster DOES percolate\n");
    }
  else
    {
      printf("percolate: cluster DOES NOT percolate\n");
    }

  /*
   *  Write the map to the file "map.pgm", displaying only the very
   *  largest cluster (or multiple clusters if exactly the same size).
   *  If the last argument here was 2, it would display the largest 2
   *  clusters etc.
   */

  //percwrite("map.pgm", map, 1);

  
  MPI_Finalize();
  return 0;
}
