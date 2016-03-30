#include "mpi.h"
#define DEBUG

/*=======================================================================*/
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*=======================================================================*/

 int main (int argc, char *argv[])

/*=======================================================================*/
/*             Begin routine                                             */
{/*begin routine */
/*=======================================================================*/
/*  Local Variables */

  int iii;
  int np,myid;
  int myid_forc;
  int myid_bead_real,myid_bead_junk;
  int np_beads      = 2;
  int np_forc       = 2;

  int i,j,k;
  int ncomm,ngroup;
  int *ranks;

  MPI_Group world_group,incl_group;
  MPI_Comm path_int_comm;
  MPI_Comm path_int_comm_junk;

/*=======================================================================*/
/* Initialize MPI */

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&np);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  myid_forc = (myid % np_forc);

  if(np_beads*np_forc != np){
    printf("Incorrect number of procs %d vs %d\n",np,np_beads*np_forc);
    MPI_Finalize();
    exit();
  }/*endfor*/

/*=======================================================================*/
/* Malloc                                                                */

  ranks  = (int *)malloc(np_beads*sizeof(int));

/*=======================================================================*/
/*             Get rank of processor in new communicator                 */

  MPI_Comm_group(MPI_COMM_WORLD,&world_group);

  for(j=0;j < np_forc;j++){
 /*-----------------------------------------------------------------------*/
 /* i) set the ranks   */
    for(i=0;i<np_beads;i++){
      ranks[i] = np_forc*i+j;
    }/*endfor*/
 /*-----------------------------------------------------------------------*/
 /* ii) Debug the ranks  */
#ifdef DEBUG
     for(k=0;k<np;k++){
       MPI_Barrier(MPI_COMM_WORLD);
       if(myid==k){
         printf("comm number %d myid %d\n",j,myid);
         for(i=0;i<np_beads;i++){ printf("ranks %d ",ranks[i]);}
         printf("\n");
       }/*endif*/
     }/*endfor*/
     if(myid==0){scanf("%d",&iii);}
     MPI_Barrier(MPI_COMM_WORLD);
#endif
 /*-----------------------------------------------------------------------*/
 /* iii) Create the new communicator                                      */
     MPI_Group_incl(world_group,np_beads,ranks,&incl_group);
     if(myid_forc==j){
       MPI_Comm_create(MPI_COMM_WORLD,incl_group,&path_int_comm);
       MPI_Comm_rank(path_int_comm,&myid_bead_real);
       MPI_Comm_size(path_int_comm,&ncomm);
#ifdef DEBUG
       printf("myid's %d %d %d\n",myid,myid_bead_real,ncomm);
#endif
       if(ncomm!=np_beads){
         printf("\n");
       }/*endif*/
     }else{
       MPI_Comm_create(MPI_COMM_WORLD,incl_group,&path_int_comm_junk);
     }/*endif*/
     MPI_Group_free(&incl_group);
#ifdef DEBUG
     if(myid==0){printf("done : \n");scanf("%d",&iii);}
#endif
     MPI_Barrier(MPI_COMM_WORLD);
   }/*endfor*/

   MPI_Group_free(&world_group);

/*==========================================================================*/
/* Exit */

   MPI_Finalize();
   exit(0);

/*------------------------------------------------------------------------*/
  } /*end routine*/ 
/*==========================================================================*/




