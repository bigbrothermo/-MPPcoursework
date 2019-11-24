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
   *  Define the main arrays for the simulation
   */
  MPI_Init(&argc,&argv);


  int old[L+2][L+2], new[L+2][L+2];

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

  int i, j, nhole, step, maxstep, oldval, newval, nchange,local_change, printfreq;
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
  int tag=0;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Status status[4];
  MPI_Request request[4];

  MPI_Status loop_status;
  MPI_Request loop_request;


  int direction,disp;
  int left,right;
  int dims[ndims];
  int period[ndims];
  int reorder;



  


  //topology
  /*
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
  */

  //2-d topology
  int dims_2d[2]={0,0};
  int periods_2d[2]={FALSE,FALSE};
  int reorder_2d=TRUE;
  int ndims_2d=2;
  

  MPI_Comm topo_comm_2d;

  MPI_Dims_create(size,ndims_2d,dims_2d);
  MPI_Cart_create(MPI_COMM_WORLD,2,dims_2d,periods_2d,reorder_2d,&topo_comm_2d);




  MPI_Barrier(MPI_COMM_WORLD);

  int M_process=dims_2d[0];
  int N_process=dims_2d[1];
  int M=L/M_process;
  int N=L/N_process;
  printf("M_process:%d N_process:%d M:%d N:%d\n",M_process,N_process,M,N);


  //int per_process_M=M/size;
  //int per_process_area=per_process_M*N;
  
  int topo_length=sqrt(size);
  int per_process_L=L/topo_length;

  int **local_map=alloc_2d_int(M,N);

  printf("each M=%d,N=%d\n",M,N);

  int per_process_area_2d=M*N;
  

  enum DIRECTIONS {DOWN,UP,LEFT,RIGHT};
  char* neighbours_names[4]={"down","up","left","right"};
  int neighbours_rank[4];

  MPI_Cart_shift(topo_comm_2d,1,1,&neighbours_rank[LEFT],&neighbours_rank[RIGHT]);
  MPI_Cart_shift(topo_comm_2d,0,1,&neighbours_rank[UP],&neighbours_rank[DOWN]);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  //printf(" rank: %d up:%d down:%d left:%d right:%d  \n",rank,neighbours_rank[UP],neighbours_rank[DOWN],neighbours_rank[LEFT],neighbours_rank[RIGHT]);

  int coords[2];
  MPI_Cart_coords(topo_comm_2d,rank,ndims_2d,coords);


  //printf("rank:%d i:%d j:%d \n",rank,coords[0],coords[1]);

  

  //create vector
  MPI_Datatype batch,tmp;
  MPI_Type_vector(N,M,L,MPI_INT,&batch);
  MPI_Type_vector(N,M,L,MPI_INT,&tmp);
  MPI_Type_create_resized(tmp,0,sizeof(int),&batch);

  MPI_Type_commit(&batch); 





  //init new vector
  /*
  MPI_datatype column_type;
  MPI_Type_vector(per_process_area,per_process_M,L,MPI_INT,&column_type);
  MPI_Type_commit(&column_type);
  */





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
  
  
  /*if(rank==0){
    for(i=0;i<M;i++){
      for(j=0;j<N;j++){
        printf("%d\t",map[i][j]);
      }
      printf("\n");
    }
  }*/
  

  

  if(rank==0){
    int index_i;
    int index_j;
    int target_coords[2];
    for(i=0;i<size;i++){
      //int target_coords[2];
      MPI_Cart_coords(topo_comm_2d,i,ndims_2d,target_coords);
      int index_i=target_coords[0]*M;
      int index_j=target_coords[1]*N;

      printf("index_i:%d index_j:%d i:%d\n",index_i,index_j,i);
      //printf("per_process_L:%d topo_length:%d\n",per_process_L,topo_length);
      MPI_Issend(&map[index_i][index_j],1,batch,i,tag,topo_comm_2d,&request[0]);
    }
  }

  MPI_Irecv(&local_map[0][0],M*N,MPI_INT,0,tag,topo_comm_2d,&request[0]);
   
  MPI_Wait(&request[0],&status[0]);

  printf("finish scatter\n");




  

  /*if(rank==2){
    display_matrix(local_map,per_process_L,per_process_L);
  }*/

  


  //MPI_Scatter(&(map[0][0]),per_process_area,MPI_INT,&(local_map[0][0]),per_process_area,MPI_INT,0,topo_comm);




  /*
   * Initialise the old array: copy the LxL array map to the centre of
   * old, and set the halo values to zero.
   */

  //int local_old[per_process_area+2][N+2];
  //int local_new[per_process_M+2][N+2];

  int ** local_old=alloc_2d_int(M+2,N+2);
  int ** local_new=alloc_2d_int(M+2,N+2);


  for (i=1; i <= M; i++)
  {
      for (j=1; j <= N; j++)
    {
      local_old[i][j] = local_map[i-1][j-1];
    }
  }

   for (i=0; i <= M+1; i++)  // zero the bottom and top halos
  {
      local_old[i][0]   = 0;
      local_old[i][N+1] = 0;
  }

  for (j=0; j <= M+1; j++)  // zero the left and right halos
  {
      local_old[0][j]   = 0;
      local_old[N+1][j] = 0;
  }




  /*
   *  Update for a fixed number of iterations
   *
   */   





  maxstep = 16*L;
  printfreq = 100;

  step = 1;
  nchange = 1;


  
  MPI_Datatype column_type;
  MPI_Type_vector(N,1,M+2,MPI_INT,&column_type);
  MPI_Type_commit(&column_type);

  //int *recv_column_left=(int *)malloc(N+2*sizeof(int));
  //int *recv_column_right=(int *)malloc(N+2*sizeof(int));

  /*
  if(rank==2){
    display_matrix(local_old,per_process_L+2,per_process_L+2);
  }
  */


  while (step <= maxstep)
  {
    local_change=0;
    //printf(" rank: %d up:%d down:%d left:%d right:%d  \n",rank,neighbours_rank[UP],neighbours_rank[DOWN],neighbours_rank[LEFT],neighbours_rank[RIGHT]);

  
    /*
    //left
    MPI_Issend(&(local_old[1][0]),1,column_type,neighbours_rank[LEFT],1,topo_comm_2d,&loop_request);
    MPI_Recv(recv_column_left,per_process_L+2,MPI_INT,neighbours_rank[LEFT],2,topo_comm_2d,&loop_request);
    //right
    MPI_Issend(&(local_old[0][per_process_L]),1,column_type,neighbours_rank[RIGHT],2,topo_comm_2d,&loop_request);
    MPI_Recv(recv_column_right,per_process_L+2,MPI_INT,neighbours_rank[RIGHT],tag,topo_comm_2d,&loop_request);
    //up
    MPI_Issend(&(local_old[1][1]),N,MPI_INT,neighbours_rank[UP],tag,topo_comm_2d,&loop_request);
    MPI_Irecv(&local_old[0][1],N,MPI_INT,neighbours_rank[UP],tag,topo_comm_2d,&loop_request);
    //down
    MPI_Issend(&(local_old[per_process_L][1]),N,MPI_INT,neighbours_rank[DOWN],tag,topo_comm_2d,&loop_request);
    MPI_Irecv(&(local_old[per_process_L+1][1]),N,MPI_INT,neighbours_rank[DOWN],tag,topo_comm_2d,&loop_request);
    */

    //send up
     MPI_Issend(&(local_old[1][1]),M,MPI_INT,neighbours_rank[UP],1,topo_comm_2d,&loop_request);
     //MPI_Recv(&local_old[per_process_L+1][1],N,MPI_INT,neighbours_rank[DOWN],1,topo_comm_2d,&loop_status);
     MPI_Irecv(&local_old[M+1][1],M,MPI_INT,neighbours_rank[DOWN],1,topo_comm_2d,&loop_request);

    //send down
     MPI_Issend(&(local_old[M][1]),M,MPI_INT,neighbours_rank[DOWN],0,topo_comm_2d,&loop_request);
     //MPI_Recv(&local_old[0][1],N,MPI_INT,neighbours_rank[UP],2,topo_comm_2d,&loop_status);
     MPI_Irecv(&local_old[0][1],M,MPI_INT,neighbours_rank[UP],0,topo_comm_2d,&loop_request);

     //send left
     MPI_Issend(&(local_old[1][1]),1,column_type,neighbours_rank[LEFT],3,topo_comm_2d,&loop_request);
     //MPI_Recv(&local_old[1][per_process_L+1],1,column_type,neighbours_rank[RIGHT],3,topo_comm_2d,&loop_status);
     MPI_Irecv(&local_old[1][N+1],1,column_type,neighbours_rank[RIGHT],3,topo_comm_2d,&loop_request);

     //send right
     MPI_Issend(&(local_old[1][N]),1,column_type,neighbours_rank[RIGHT],4,topo_comm_2d,&loop_request);
     //MPI_Recv(&local_old[1][0],1,column_type,neighbours_rank[LEFT],4,topo_comm_2d,&loop_status);
     MPI_Irecv(&local_old[1][0],1,column_type,neighbours_rank[LEFT],4,topo_comm_2d,&loop_request);


     //printf("rank %d's right is %d\t",neighbours_rank[RIGHT]);






     /*
      * overlap calculation
      *
      */
     for (i=2; i<=M-1; i++)
     {
       for (j=2; j<=N-1; j++)
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








    MPI_Wait(&loop_request,&loop_status);

    //printf("finish halo exchange in rank %d\n",rank);
    MPI_Barrier(MPI_COMM_WORLD);
 
/*
  if(step==1&&rank==2){
      for(i=0;i<per_process_L+2;i++){
        printf("%d\t",local_old[1][i]);
      }
      printf("\n");

     }

*/

    /*
      *halo calculation
      *
      */

    for(i=1;i<=M;i++){
      oldval=local_old[i][1];
      newval=oldval;
        /*
         * Set local[i][j] to be the maximum value of local_old[i][j]
         * and its four nearest neighbours
         */

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


          oldval=local_old[i][N];
          newval=oldval;
          if (oldval != 0)
          {
            if (local_old[i-1][N] > newval) newval =local_old[i-1][N];
            if (local_old[i+1][N] > newval) newval = local_old[i+1][N];
            if (local_old[i][N-1] > newval) newval = local_old[i][N-1];
            if (local_old[i][N+1] > newval) newval = local_old[i][N+1];

            if (newval != oldval)
            {
                ++local_change;
            }
          }
          local_new[i][N] = newval;

    }

    for(j=1;j<=N;j++){
      oldval=local_old[1][j];
      newval=oldval;
        /*
         * Set local[i][j] to be the maximum value of local_old[i][j]
         * and its four nearest neighbours
         */

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


          oldval=local_old[M][j];
          newval=oldval;
          if (oldval != 0)
          {
            if (local_old[M-1][j] > newval) newval =local_old[M-1][j];
            if (local_old[M+1][j] > newval) newval = local_old[M+1][j];
            if (local_old[M][j-1] > newval) newval = local_old[M][j-1];
            if (local_old[M][j+1] > newval) newval = local_old[M][j+1];

            if (newval != oldval)
            {
                ++local_change;
            }
          }

          local_new[M][j]=newval;



    }




   
 /*
    for (i=1; i<=per_process_L; i++)
     {
       for (j=1; j<=per_process_L; j++)
       {
         oldval = local_old[i][j];
         newval = oldval;

        /*
         * Set local[i][j] to be the maximum value of local_old[i][j]
         * and its four nearest neighbours
        

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

    */
    MPI_Barrier(topo_comm_2d);
    //printf("local_change in rank %d is %d\n",rank,local_change);


    MPI_Reduce(&local_change,&nchange,1,MPI_INT,MPI_SUM,0,topo_comm_2d);

    //if(rank==0){
    //  printf("nchange is %d \n",nchange);
    //}


    

      /*
       *  Copy back in preparation for next step, omitting halos
       */
 
      for (i=1; i<=M; i++)
      {
        for (j=1; j<=N; j++)
        {
        local_old[i][j] = local_new[i][j];
        }
      }

      /*
       *  Report progress every now and then
       */

      
      /*if (step % printfreq == 0)
      {
        printf("percolate: number of changes on step %d is %d\n",
        step, nchange);
      }*/
      

   
  //when nchang ==0 finissssh the loop

    
    MPI_Barrier(topo_comm_2d);
    MPI_Status state;
    //int stop_flag=0;
    if(rank==0){
      //printf("current nchange is %d",nchange);
      if (step % printfreq == 0 || nchange==0)
      {
        printf("percolate: number of changes on step %d is %d\n",
        step, nchange);
      }



      for(i=1;i<size;i++){
            MPI_Ssend(&nchange,1,MPI_INT,i,tag,topo_comm_2d);
            //printf("finish sending nchange in rank %d\n, send to %d",nchange,i);
          }
      
    }

    else{
      MPI_Recv(&nchange,1,MPI_INT,0,tag,topo_comm_2d,&state);
    }

  
    MPI_Barrier(topo_comm_2d);
    if(nchange==0){
      //printf("break from rank %d",rank);

  /*if(rank==2){
      display_matrix(local_old,per_process_L+2,per_process_L+2);
    }*/

      break;
    }




 //}   


  /*for (i=0; i <= per_process_M+1; i++)  // zero the bottom and top halos
  {
      local_old[i][0]   = 0;
      local_old[i][N+1] = 0;
  }

  for (j=0; j <= N+1; j++)  // zero the left and right halos
  {
      local_old[0][j]   = 0;
      local_old[per_process_M+1][j] = 0;
  }  
    */  
  
        MPI_Barrier(topo_comm_2d);
        step++;
  }
  




   /*
    *  Update for a fixed number of iterations
    */

if (nchange != 0 && rank==0)
    {
      printf("percolate: WARNING max steps = %d reached before nchange = 0\n",
       maxstep);
    }

  /*
   *  Copy the centre of old, excluding the halos, into map
   */
  
  for (i=1; i<=M; i++)
    {
      for (j=1; j<=N; j++)
      {
         local_map[i-1][j-1] = local_old[i][j];
      }
    }


    /*if(rank==2){
      display_matrix(local_old,per_process_L+2,per_process_L+2);
    }*/
    
  //copy back to rank 0

    MPI_Issend(&local_map[0][0],M*N,MPI_INT,0,tag,topo_comm_2d,&request[0]);
    
    if(rank==0){
      int index_i;
      int index_j;
      int target_coords[2];
      //display_matrix(map,per_process_L,per_process_L);
      for(i=0;i<size;i++){
        MPI_Cart_coords(topo_comm_2d,i,ndims_2d,target_coords);
        index_i=target_coords[0]*M;
        index_j=target_coords[1]*N;
        //printf("index_i:%d index_j:%d i:%d \n",index_i,index_j,i);
        MPI_Recv(&map[index_i][index_j],1,batch,i,tag,topo_comm_2d,&status[0]);
        //printf("finish recv from rank %d\n",i);
        //MPI_Irecv(&map[index_i][index_j],1,batch,i,tag,topo_comm_2d,&request[0]);


      }
      //MPI_Issend(&local_map[0][0],per_process_L*per_process_L,MPI_INT,0,tag,topo_comm_2d,&request[0]);

      MPI_Wait(&request[0],&status[0]);   
      printf("finish copy back to rank 0\n"); 


    }

    
  //MPI_Gather(&(local_map[0][0]),per_process_area,MPI_INT,&(map[0][rank*per_process_M]),per_process_area,MPI_INT,0,MPI_COMM_WORLD);
  

  /*
  MPI_Gather(&(local_map[0][0]),per_process_area,MPI_INT,&(map[0][rank*per_process_M]),per_process_area,MPI_INT,0,topo_comm);
  */

  /*
  if(rank==0){
    for(i=per_process_M;i<per_process_M*2;i++){
      for(j=0;j<N;j++){
        printf("%d\t",map[i][j]);
      }
      printf("\n");
    }
  }
  */





  //
  
  //free(local_map[0]);
  //free(local_map);
  MPI_Barrier(MPI_COMM_WORLD);
  printf("free successfully\n");

  if(rank==0){
  /*  
   *  Test to see if percolation occurred by looking for positive numbers
   *  that appear on both the top and bottom edges
   */


    
    /*for(i=0;i<M;i++){
      for(j=0;j<N;j++){
        printf("%d\t",map[i][j]);
      }
      printf("\n");
    }*/
    
  
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

  percwrite("map.pgm", map, 1);
  }
  

  

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
