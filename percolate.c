#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "percolate.h"
#include <math.h>

/*
 * Serial program to test for percolation of a cluster.
 */



int main(int argc, char *argv[])
{
   /*
    *  MPI init, declare basic variable for message-passing 
    *
    *
    */

  MPI_Init(&argc,&argv);
  MPI_Barrier(MPI_COMM_WORLD);
  double time_init=MPI_Wtime();
  
  
  int rank,size;
  int tag=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Status status;
  MPI_Request request;

  MPI_Status loop_status;
  MPI_Request loop_request;

   /*
   *  Define the main arrays for the simulation
   */
    //int old[L+2][L+2], new[L+2][L+2];

  /*
   *  Additional array WITHOUT halos for initialisation and IO. This
   *  is of size LxL because, even in our parallel program, we do
   *  these two steps in serial
   */ 

  int L=288;

  int ** map=alloc_2d_int(L,L);

  /*
   *  Variables that define the simulation
   */

  int seed=1564;
  double rho;

  

  /*
   *  Local variables
   */

  int i, j, nhole, step, maxstep, oldval, newval, nchange,local_change, printfreq;
  int sum,local_sum;
  double average=0;
  int itop, ibot, perc;
  double r;

  /*if (argc != 2)
    {
      printf("Usage: percolate <seed>\n");
      return 1;
    }
    */

  /*
   *  Set most important value: the rock density rho (between 0 and 1)
   */

  rho = 0.411;

  /*
   *  Set the randum number seed and initialise the generator
   */

  //seed = atoi(argv[1]);


  int k=1;
  while(k<argc){
    if(argv[k][0]!='-'){
      printf("Argument # Error : Arguments must start with -\n");
      exit(0);
    }
  switch(argv[k][1]){
      case 'L':L=atoi(argv[k+1]);break;
      case 'r':rho=atof(argv[k+1]);break;
      case 's':seed=atoi(argv[k+1]);break;
      default : printf("Unrec argumen-t %s\n", argv[k]); break;
    }

    k=k+2;
    
  }  

//inputdata filter
  /*
  if(filter_arguments(&grid)==1){
    exit(1);
  }

  if(filter_output(&output)==1){
    exit(1);
  }  */


  if(rank==0)
    printf("percolate: params are L = %d, rho = %f, seed = %d\n", L, rho, seed);

  rinit(seed);





 
  /*
   * Variables that defines the 2D Cartesian Virtual Topology
   */

  //int direction,disp;
  //int left,right;
  //int dims[ndims];
  //int period[ndims];
  //int reorder;

  int dims_2d[2]={0,0};
  int periods_2d[2]={TRUE,FALSE};
  int reorder_2d=TRUE;
  int ndims_2d=2;
  

  MPI_Comm topo_comm_2d;


   /*
    * Create the 2D Cartesian Virtual Topology
    */
  MPI_Dims_create(size,ndims_2d,dims_2d);
  MPI_Cart_create(MPI_COMM_WORLD,2,dims_2d,periods_2d,reorder_2d,&topo_comm_2d);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank); 

  enum DIRECTIONS {DOWN,UP,LEFT,RIGHT};
  //char* neighbours_names[4]={"down","up","left","right"};
  int neighbours_rank[4];

  MPI_Cart_shift(topo_comm_2d,1,1,&neighbours_rank[LEFT],&neighbours_rank[RIGHT]);
  MPI_Cart_shift(topo_comm_2d,0,1,&neighbours_rank[UP],&neighbours_rank[DOWN]);
  



   /* 
    *   Get the M and N for each process and allocate memory for map in each process
    */
  int M_process=dims_2d[0];
  int N_process=dims_2d[1];
  int M=L/M_process;
  int N=L/N_process;

  int coords[2];
  MPI_Cart_coords(topo_comm_2d,rank,ndims_2d,coords);
  int own_M=M;
  int own_N=N;
  int last_M=L-(M_process-1)*M;
  int last_N=L-(N_process-1)*N;
  if(coords[0]==M_process-1){
    own_M=L-coords[0]*M;

  }
  if(coords[1]==N_process-1){
    own_N=L-coords[1]*own_N;

  }

  printf("each M=%d,N=%d in rank %d\n",own_M,own_N,rank);

  int **local_map=alloc_2d_int(own_M,own_N);

  //printf("own_M in rank %d is %d\n",rank,own_M);



  /*
   *  Initialise map with density rho. Zero indicates rock, a positive
   *  value indicates a hole. For the algorithm to work, all the holes
   *  must be initialised with a unique integer
   */

  
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

  }


  MPI_Barrier(MPI_COMM_WORLD);
  double time_before_scatter=MPI_Wtime();
  /*
   * Scatter the map's data in rank 0 to each process,specify the M and N in each sending
   */
  if(rank==0){
    int index_i;
    int index_j;
    int send_M;//Send different M and N to each process 
    int send_N;
    int target_coords[2];
    for(i=0;i<size;i++){
      
      MPI_Cart_coords(topo_comm_2d,i,ndims_2d,target_coords);
      index_i=target_coords[0]*M;
      index_j=target_coords[1]*N;

      send_M=M;
      send_N=N;
      if(target_coords[0]==M_process-1){
          send_M=last_M;
      }
      if(target_coords[1]==N_process-1){
        send_N=last_N;

      }

      //create vector
      MPI_Datatype batch;
      MPI_Type_vector(send_M,send_N,L, MPI_INT,&batch);
      MPI_Type_commit(&batch); 

      MPI_Issend(&map[index_i][index_j],1,batch,i,tag,topo_comm_2d,&request);
    }
  }

  MPI_Irecv(&local_map[0][0],own_M*own_N,MPI_INT,0,tag,topo_comm_2d,&request);
   
  MPI_Wait(&request,&status);

  MPI_Barrier(MPI_COMM_WORLD);
  double time_finish_scatter=MPI_Wtime();



  /*
   * Initialise the old array: copy the LxL array map to the centre of
   * old, and set the halo values to zero.
   */


  int ** local_old=alloc_2d_int(own_M+2,own_N+2);
  int ** local_new=alloc_2d_int(own_M+2,own_N+2);


  for (i=1; i <= own_M; i++)
  {
      for (j=1; j <= own_N; j++)
    {
      local_old[i][j] = local_map[i-1][j-1];
    }
  }

   for (i=0; i <= own_M+1; i++)  // zero the bottom and top halos
  {
      local_old[i][0]   = 0;
      local_old[i][own_N+1] = 0;
  }

  for (j=0; j <= own_N+1; j++)  // zero the left and right halos
  {
      local_old[0][j]   = 0;
      local_old[own_M+1][j] = 0;
  }




  /*
   *  Update for a fixed number of iterations
   *
   */   





  maxstep = 16*L;
  printfreq = 100;

  step = 1;
  nchange = 1;
 

  /*Create vector to send column at left and right side of local old*/
  MPI_Datatype column_type;
  MPI_Type_vector(own_M,1,own_N+2,MPI_INT,&column_type);
  MPI_Type_commit(&column_type);


  MPI_Barrier(MPI_COMM_WORLD);




  MPI_Barrier(MPI_COMM_WORLD);
  double time_before_loop=MPI_Wtime();

  while (step <= maxstep)
  {
    local_change=0;

     /*
      * halo exchange
      *
      */


    //send up
     MPI_Issend(&(local_old[1][1]),own_N,MPI_INT,neighbours_rank[UP],1,topo_comm_2d,&loop_request);
     MPI_Irecv(&local_old[own_M+1][1],own_N,MPI_INT,neighbours_rank[DOWN],1,topo_comm_2d,&loop_request);

    //send down
     MPI_Issend(&(local_old[own_M][1]),own_N,MPI_INT,neighbours_rank[DOWN],0,topo_comm_2d,&loop_request);
     MPI_Irecv(&local_old[0][1],own_N,MPI_INT,neighbours_rank[UP],0,topo_comm_2d,&loop_request);

     //send left
     MPI_Issend(&(local_old[1][1]),1,column_type,neighbours_rank[LEFT],3,topo_comm_2d,&loop_request);
     MPI_Irecv(&local_old[1][own_N+1],1,column_type,neighbours_rank[RIGHT],3,topo_comm_2d,&loop_request);

     //send right
     MPI_Issend(&(local_old[1][own_N]),1,column_type,neighbours_rank[RIGHT],4,topo_comm_2d,&loop_request);
     MPI_Irecv(&local_old[1][0],1,column_type,neighbours_rank[LEFT],4,topo_comm_2d,&loop_request);



     /*
      * overlap calculation
      *
      */
     for (i=2; i<=own_M-1; i++)
     {
       for (j=2; j<=own_N-1; j++)
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
                ++local_change;
            }
          }

        local_new[i][j] = newval;
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);




     /*
      * Finish overlap calculation and wait the finish of halo exchange
      *
      */

    MPI_Wait(&loop_request,&loop_status);

    MPI_Barrier(MPI_COMM_WORLD);
 

     /*
      * Conduct number exchange on the left and right side of local map
      */

    for(i=1;i<=own_M;i++){
      oldval=local_old[i][1];
      newval=oldval;

      if (oldval != 0)
          {
            if (local_old[i-1][1] > newval) newval =local_old[i-1][1];
            if (local_old[i+1][1] > newval) newval = local_old[i+1][1];
            if (local_old[i][0] > newval) newval = local_old[i][0];
            if (local_old[i][2] > newval) newval = local_old[i][2];

            if (newval != oldval)
            {
                ++local_change;

            }
          }
          local_new[i][1]=newval;


          oldval=local_old[i][own_N];
          newval=oldval;
          if (oldval != 0)
          {
            if (local_old[i-1][own_N] > newval) newval =local_old[i-1][own_N];
            if (local_old[i+1][own_N] > newval) newval = local_old[i+1][own_N];
            if (local_old[i][own_N-1] > newval) newval = local_old[i][own_N-1];
            if (local_old[i][own_N+1] > newval) newval = local_old[i][own_N+1];

            if (newval != oldval)
            {
                ++local_change;
              
            }
          }
          local_new[i][own_N] = newval;

    }


     /*
      * Conduct number exchange on the up and down side of local map
      */

    for(j=1;j<=own_N;j++){
      oldval=local_old[1][j];
      newval=oldval;
  

      if (oldval != 0)
          {
            if (local_old[0][j] > newval) newval =local_old[0][j];
            if (local_old[2][j] > newval) newval = local_old[2][j];
            if (local_old[1][j-1] > newval) newval = local_old[1][j-1];
            if (local_old[1][j+1] > newval) newval = local_old[1][j+1];

            if (newval != oldval)
            {
                ++local_change;
                
            }
          }

          local_new[1][j] = newval;


          oldval=local_old[own_M][j];
          newval=oldval;
          if (oldval != 0)
          {
            if (local_old[own_M-1][j] > newval) newval =local_old[own_M-1][j];
            if (local_old[own_M+1][j] > newval) newval = local_old[own_M+1][j];
            if (local_old[own_M][j-1] > newval) newval = local_old[own_M][j-1];
            if (local_old[own_M][j+1] > newval) newval = local_old[own_M][j+1];

            if (newval != oldval)
            {
                ++local_change;
            }
          }

          local_new[own_M][j]=newval;

    }


    
    MPI_Barrier(topo_comm_2d);//make sure every process processes the same step before geting nchange

    MPI_Reduce(&local_change,&nchange,1,MPI_INT,MPI_SUM,0,topo_comm_2d);


      /*
       *  Copy back in preparation for next step, omitting halos
       */
      local_sum=0;
      for (i=1; i<=own_M; i++)
      {
        for (j=1; j<=own_N; j++)
        {
        local_old[i][j] = local_new[i][j];
        local_sum=local_sum+local_new[i][j];
        }
      }

      /*
       * Calculate the sum of each small matrix and sent back to rank 0 to get average
       */
      

      MPI_Reduce(&local_sum,&sum,1,MPI_INT,MPI_SUM,0,topo_comm_2d);
      if(rank==0){
        average=(double)sum/(L*L);
      }
      



      /*
       *  Report progress every 100 steps, providing 
       */    
    MPI_Barrier(topo_comm_2d);
    MPI_Status state;

    if(rank==0){
      if (step % printfreq == 0 || nchange==0)
      {
        printf("percolate: number of changes on step %d is %d\n",
        step, nchange);
        printf("Average: current average value of grid on step %d is %f\n",step,average);
      }


    /*
     * Send nchange result to each process, if nchange is 0, then end the loop 
     */
    for(i=1;i<size;i++){
            MPI_Ssend(&nchange,1,MPI_INT,i,tag,topo_comm_2d);
      
          }
      
    }

    else{
      MPI_Recv(&nchange,1,MPI_INT,0,tag,topo_comm_2d,&state);
    }

  
    MPI_Barrier(topo_comm_2d);
    if(nchange==0){

      break;
    }
     
  
    MPI_Barrier(topo_comm_2d);
    step++;
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  double time_finish_loop=MPI_Wtime();



   /*
    *  Update for a fixed number of iterations
    */

if (nchange != 0 && rank==0)
    {
      printf("percolate: WARNING max steps = %d reached before nchange = 0\n",
       maxstep);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double time_before_copy_back=MPI_Wtime(); 

  /*
   *  Copy the centre of old, excluding the halos, into map
   */
  
  for (i=1; i<=own_M; i++)
    {
      for (j=1; j<=own_N; j++)
      {
         local_map[i-1][j-1] = local_old[i][j];
      }
    }


    
  /*
   * Send the information of local map in each process back to the map in rank 0
   */
    MPI_Issend(&local_map[0][0],own_M*own_N,MPI_INT,0,tag,topo_comm_2d,&request);
    
    if(rank==0){
      int index_i;
      int index_j;
      int target_coords[2];
      int send_M;
      int send_N;
      for(i=0;i<size;i++){
        MPI_Cart_coords(topo_comm_2d,i,ndims_2d,target_coords);
        index_i=target_coords[0]*M;
        index_j=target_coords[1]*N;


        send_M=M;
        send_N=N;
        if(target_coords[0]==M_process-1){
           send_M=last_M;
        }
        if(target_coords[1]==N_process-1){
          send_N=last_N;
        }
        //create vector
        MPI_Datatype batch;
        MPI_Type_vector(send_M,send_N,L, MPI_INT,&batch);


        MPI_Type_commit(&batch);
        MPI_Recv(&map[index_i][index_j],1,batch,i,tag,topo_comm_2d,&status);



      }


      MPI_Wait(&request,&status);   


    }

    

  MPI_Barrier(MPI_COMM_WORLD);
  double time_finish_copy_back=MPI_Wtime(); 



  if(rank==0){
  /*  
   *  Test to see if percolation occurred by looking for positive numbers
   *  that appear on both the top and bottom edges
   */

  
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

  if(perc != 0)
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

  percwrite("map.pgm", map, L, 1);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  double time_finish_output=MPI_Wtime();

  if(rank==0){
    double init_time=time_before_scatter-time_init;
    double scatter_time=time_finish_scatter - time_before_scatter;
    double loop_time= time_finish_loop-time_before_loop;
    double copy_back_time=time_finish_copy_back - time_before_copy_back;

    printf("Runtime\n");
    printf("Init_time:%f\n",init_time);
    printf("Scatter_time:%f\n",scatter_time);
    printf("Loop_time:%f\n",loop_time);
    printf("Copy_back_time:%f\n",copy_back_time);

  }
  

  
  MPI_Finalize();
  return 0;
}
