/*==========================================================================*/
/*CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC*/
/*==========================================================================*/
/*                                                                          */
/*                         Stochastic DFT:                                  */
/*             The future of density functional theory                      */
/*             ------------------------------------                         */
/*                   Module: init-frag.c                                    */
/*                                                                          */
/* This routine initialize the fragmentation calculation (common part).     */
/*                                                                          */
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/

#include "standard_include.h"
#include "../typ_defs/typedefs_gen.h"
#include "../typ_defs/typedefs_class.h"
#include "../typ_defs/typedefs_bnd.h"
#include "../typ_defs/typedefs_cp.h"
#include "../typ_defs/typedefs_stat.h"
#include "../typ_defs/typedefs_par.h"
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_frag_local.h"
#include "../proto_defs/proto_interface_frag_local.h"

#define TIME_CP_OFF
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void initFrag(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,CP *cp,
	      ANALYSIS *analysis,
	      CLASS **classMiniPoint,BONDED **bondedMiniPoint,
	      GENERAL_DATA **generalDataMiniPoint,
	      ANALYSIS **analysisMiniPoint,CP **cpMiniPoint,int ip_now)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Main driver for fragmentation initialization common part              */
/* This should be done in the stochastic dft initialization		 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"

  STODFTINFO    *stodftInfo       = cp->stodftInfo;
  FRAGINFO      *fragInfo;
  COMMUNICATE   *communicate      = &(cp->communicate);

  int fragOpt           = stodftInfo->fragOpt;
  int numFragProc;
  int myidState		= communicate->myid_state;
  int numProcStates	= communicate->np_states;
  int iFrag;
  MPI_Comm world                = communicate->world;

  if(myidState==0)fragInfo = stodftInfo->fragInfo;
  else{
    stodftInfo->fragInfo = (FRAGINFO*)cmalloc(sizeof(FRAGINFO));
    fragInfo = stodftInfo->fragInfo;
  }
  if(numProcStates>1)Barrier(world);
  //printf("fragOpt %i\n",fragOpt);
  
  switch(fragOpt){ 
    case 1:
      initFragMol(class,bonded,general_data,cp,ip_now);
      break;
    case 3:
      initFragUnitCell(class,bonded,general_data,cp,ip_now); // Use the same initialization
      break;
  }


  numFragProc = fragInfo->numFragProc;
  *classMiniPoint = (CLASS*)cmalloc(numFragProc*sizeof(CLASS));
  *bondedMiniPoint = (BONDED*)cmalloc(numFragProc*sizeof(BONDED));
  *generalDataMiniPoint = (GENERAL_DATA*)cmalloc(numFragProc*sizeof(GENERAL_DATA));
  *analysisMiniPoint = (ANALYSIS*)cmalloc(numFragProc*sizeof(ANALYSIS));
  *cpMiniPoint = (CP*)cmalloc(numFragProc*sizeof(CP));
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    //printf("iFrag %i numFragProc %i\n",iFrag,numFragProc);
    fragInfo->iFrag = iFrag;
    parseFrag(class,bonded,general_data,cp,analysis,&((*classMiniPoint)[iFrag]),
	      &((*bondedMiniPoint)[iFrag]),&((*generalDataMiniPoint)[iFrag]),
	      &((*cpMiniPoint)[iFrag]),&((*analysisMiniPoint)[iFrag]));
  }
  CLASS *classMini = (*classMiniPoint);
  /*
  printf("xyz %lg %lg %lg\n",classMini[0].clatoms_pos[1].x[1],
	classMini[0].clatoms_pos[1].y[1],classMini[0].clatoms_pos[1].z[1]);
  printf("xyz %lg %lg %lg\n",classMini[0].clatoms_pos[1].x[2],
        classMini[0].clatoms_pos[1].y[2],classMini[0].clatoms_pos[1].z[2]);
  printf("xyz %lg %lg %lg\n",classMini[0].clatoms_pos[1].x[3],
        classMini[0].clatoms_pos[1].y[3],classMini[0].clatoms_pos[1].z[3]);
  */
  //exit(0);
 
  if(numFragProc>0){
    if(fragOpt==1){
      initFragEnergy(*cpMiniPoint,*classMiniPoint,class,cp);
    }
  }

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void initFragMol(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,CP *cp,
		 int ip_now)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Fragmentation initialization common part for molecular fragments      */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"

  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  CLATOMS_POS  *clatoms_pos  = &(class->clatoms_pos[ip_now]);
  ATOMMAPS     *atommaps     = &(class->atommaps);
  CELL         *cell         = &(general_data->cell);
  FOR_SCR      *for_scr      = &(class->for_scr);
  EWD_SCR      *ewd_scr      = &(class->ewd_scr);
  PTENS        *ptens        = &(general_data->ptens);

  CPOPTS        *cpopts           = &(cp->cpopts);
  PSEUDO        *pseudo           = &(cp->pseudo);
  CPCOEFFS_INFO *cpcoeffs_info    = &(cp->cpcoeffs_info);
  COMMUNICATE   *communicate      = &(cp->communicate);
  CPCOEFFS_POS  *cpcoeffs_pos     = &(cp->cpcoeffs_pos[ip_now]);
  STODFTINFO    *stodftInfo       = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos    = cp->stodftCoefPos;
  CPSCR         *cpscr            = &(cp->cpscr);
  NEWTONINFO    *newtonInfo;
  FRAGINFO      *fragInfo	  = stodftInfo->fragInfo;
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);

  int iChem,iSamp,iCell,iProc,iCoeff,iMol,iFrag,iType,iAtom;
  int count;
  int cpAtom;
  int numFragTot;
  int numFragProc;
  int div,res;
  int fragOpt           = stodftInfo->fragOpt;
  int fragCellOpt       = stodftInfo->fragCellOpt;
  int rhoRealGridNum	= stodftInfo->rhoRealGridNum;
  int cpLsda		= cpopts->cp_lsda;
  int numMolTot;
  int numMolType        = atommaps->nmol_typ;
  int molIndStart;
  int atomIndTypeStart;
  int countMol,countAtom;
  int molInd,atomInd;
  int numAtomQM		= clatoms_info->nab_initio;
  int numAtomTot	= clatoms_info->natm_tot;
  int myidState		= communicate->myid_state;
  int numProcStates	= communicate->np_states;
  MPI_Comm commStates   = communicate->comm_states;

  int *numMolJmolType	   = atommaps->nmol_jmol_typ;
  int *numAtomJmolType	   = atommaps->natm_1mol_jmol_typ;
  int *jatomJmolTypeStart  = atommaps->jatm_jmol_typ_strt;
  int *numMolFragProc;
  int *numAtomFragProc;
  int *numAtomFragAll;
  int *atomIndStart;
  int *atomNumMol;
  int *cpVlncUp		    = clatoms_info->cp_vlnc_up;
  int *cpVlncDn		    = clatoms_info->cp_vlnc_dn;
  int *cpAtomList	    = atommaps->cp_atm_lst;
  int *atomVlncUpAll,*atomVlncDnAll;
  int *numElecUpFragTot,*numElecDnFragTot,*numElecUpFragProc,*numElecDnFragProc;
  int *molTypeMapAll,*molTypeMapFrag,*numMolTypeFrag;

  int **molTypeFrag,**molNumTypeFrag;
  int **molFragMapProc;
  int **atomFragMapProc;

  char *atomSkinFile;

  double *skinAll;
  FILE *fileSkin;
  FILE *fileNumUC;

/*======================================================================*/
/* I) Get number of fragments per processor                             */

  if(myidState==0){
    printf("**Get number of fragments\n");
    printf("%s\n",fragInfo->atomSkinFile);
  }

  numMolTot = 0;
  for(iMol=1;iMol<=numMolType;iMol++)numMolTot += numMolJmolType[iMol];
  fragInfo->numFragTot = numMolTot;
  numFragTot = fragInfo->numFragTot;
  // I will put in other options later
  div = numFragTot/numProcStates;
  res = numFragTot%numProcStates;
  if(myidState<res)fragInfo->numFragProc = div+1;
  else fragInfo->numFragProc = div;
  numFragProc = fragInfo->numFragProc;

/*======================================================================*/
/* II) Get mol map in each fragments                                    */

  if(myidState==0){
    printf("**Get mol map\n");
  }
  
  // Get Mol Map
  fragInfo->numMolFragProc = (int*)cmalloc(numFragProc*sizeof(int));
  numMolFragProc = fragInfo->numMolFragProc;
  for(iFrag=0;iFrag<numFragProc;iFrag++)numMolFragProc[iFrag] = 1;
  
  fragInfo->molFragMapProc = (int**)cmalloc(numFragProc*sizeof(int*));
  molFragMapProc = fragInfo->molFragMapProc;
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    molFragMapProc[iFrag] = (int*)cmalloc(numMolFragProc[iFrag]*sizeof(int));
  }

  if(myidState<res)molIndStart = myidState*(div+1);
  else molIndStart = res*(div+1)+(myidState-res)*div;
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    molFragMapProc[iFrag][0] = iFrag+molIndStart+1;
  }
  molTypeMapAll = (int*)cmalloc(numMolTot*sizeof(int));
  count = 0;
  // Get Mol Type
  // First get map molTypeMapAll:(molecular index)->(molecular type) for all molecules
  for(iType=1;iType<=numMolType;iType++){
    for(iMol=0;iMol<numMolJmolType[iType];iMol++){
      molTypeMapAll[count] = iType;
      count += 1;
    }
  }
  //for(iMol=0;iMol<numMolTot;iMol++)printf("iMol %i molTypeMapAll %i\n",iMol,molTypeMapAll[iMol]);
  //Second Get number of molecule type in each fragment on this proc
  fragInfo->numMolTypeFrag = (int*)cmalloc(numFragProc*sizeof(int));
  numMolTypeFrag = fragInfo->numMolTypeFrag;
  for(iFrag=0;iFrag<numFragProc;iFrag++)numMolTypeFrag[iFrag] = 1;
  //Get the list of molecule type in each fragment(molTypeFrag) and 
  //the number of molecules in each molecule type in each fragment(molNumTypeFrag)
  fragInfo->molTypeFrag = (int**)cmalloc(numFragProc*sizeof(int*));
  fragInfo->molNumTypeFrag = (int**)cmalloc(numFragProc*sizeof(int*));
  molTypeFrag = fragInfo->molTypeFrag;
  molNumTypeFrag = fragInfo->molNumTypeFrag;
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    molTypeFrag[iFrag] = (int*)cmalloc(numMolTypeFrag[iFrag]*sizeof(int));
    molNumTypeFrag[iFrag] = (int*)cmalloc(numMolTypeFrag[iFrag]*sizeof(int));
  }
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    molTypeMapFrag = (int*)cmalloc(numMolFragProc[iFrag]*sizeof(int));
    for(iMol=0;iMol<numMolFragProc[iFrag];iMol++){
      molTypeMapFrag[iMol] = molTypeMapAll[molFragMapProc[iFrag][iMol]-1];
      //printf("iFrag %i iMol %i molFragMapProc %i molTypeMapFrag %i\n",iFrag,iMol,molFragMapProc[iFrag][iMol],molTypeMapFrag[iMol]);
    }
    molTypeFrag[iFrag][0] = molTypeMapFrag[0];
    molNumTypeFrag[iFrag][0] = 1;
    
    for(iMol=1;iMol<numMolFragProc[iFrag];iMol++){
      if(molTypeFrag[iFrag][numMolTypeFrag[iFrag]-1]==molTypeMapFrag[iMol]){
	// If this molecule is in the same type, this type has one more member
	molNumTypeFrag[iFrag][numMolTypeFrag[iFrag]-1] += 1;
      }
      else{
	// If molecular type does not equal to the old one, add the new mol type
	numMolTypeFrag[iFrag] += 1;
	molTypeFrag[iFrag] = (int*)realloc(molTypeFrag[iFrag],numMolTypeFrag[iFrag]);
	molTypeFrag[iFrag][numMolTypeFrag[iFrag]-1] = molTypeMapFrag[iMol];
	molNumTypeFrag[iFrag] = (int*)realloc(molTypeFrag[iFrag],numMolTypeFrag[iFrag]);
	molNumTypeFrag[iFrag][numMolTypeFrag[iFrag]-1] = 1;
      }//endif
    }//endfor iMol       
    free(molTypeMapFrag);
  }//endfor iFrag
  /*
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    printf("iFrag %i molTypeFrag %i molNumTypeFrag %i\n",iFrag,molTypeFrag[iFrag][0],molNumTypeFrag[iFrag][0]);
  }
  exit(0);
  */
  
/*======================================================================*/
/* 2) Get atoms map in each fragments                                   */

  if(myidState==0){
    printf("**Get atom map\n");
  }
  
  atomIndStart = (int*)cmalloc(numMolTot*sizeof(int));
  atomNumMol = (int*)cmalloc(numMolTot*sizeof(int));
  // Get starting atom index for each molecule
  countMol = 0;
  for(iType=1;iType<=numMolType;iType++){
    atomIndTypeStart = jatomJmolTypeStart[iType];
    for(iMol=0;iMol<numMolJmolType[iType];iMol++){
      atomIndStart[countMol+iMol] = atomIndTypeStart+iMol*numAtomJmolType[iType];
      //printf("countMol+iMol %i,atomIndStart %i\n",countMol+iMol,atomIndStart[countMol+iMol]);
      atomNumMol[countMol+iMol] = numAtomJmolType[iType];
    }
    countMol += numMolJmolType[iType];
  }
  // Get Atom number per fragment
  fragInfo->numAtomFragProc = (int*)cmalloc(numFragProc*sizeof(int));
  numAtomFragProc = fragInfo->numAtomFragProc;
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    numAtomFragProc[iFrag] = 0;
    for(iMol=0;iMol<numMolFragProc[iFrag];iMol++){
      numAtomFragProc[iFrag] += atomNumMol[molFragMapProc[iFrag][iMol]-1];
    }
  }
  // Get Atom Map
  fragInfo->atomFragMapProc = (int**)cmalloc(numFragProc*sizeof(int*));
  atomFragMapProc = fragInfo->atomFragMapProc;
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    atomFragMapProc[iFrag] = (int*)cmalloc(numAtomFragProc[iFrag]*sizeof(int));
  }
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    countAtom = 0;
    for(iMol=0;iMol<numMolFragProc[iFrag];iMol++){
      for(iAtom=0;iAtom<atomNumMol[molFragMapProc[iFrag][iMol]-1];iAtom++){
	atomFragMapProc[iFrag][countAtom+iAtom] = atomIndStart[molFragMapProc[iFrag][iMol]-1]+iAtom;
      }//endfor iAtom
      countAtom += atomNumMol[molFragMapProc[iFrag][iMol]-1];
    }//endfor iMol
  }//endfor iFrag

  /*
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    //printf("iFrag %i numAtomFragProc %i\n",iFrag,numAtomFragProc[iFrag]);
    printf("iFrag %i atomFragMap St %i atomFragMap Ed %i\n",iFrag,atomFragMapProc[iFrag][0],atomFragMapProc[iFrag][numAtomFragProc[iFrag]-1]);
  }
  */

/*======================================================================*/
/* 3) Get electron number in each fragments                             */

  if(myidState==0){
    printf("**Get electron number\n");
  }

  if(numAtomQM!=numAtomTot){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("I am not interested in calculating 'QMMM'. Let all\n");
    printf("atoms be ab initial!\n");
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(0);
  }
  atomVlncUpAll = (int*)cmalloc(numAtomTot*sizeof(int));
  atomVlncDnAll = (int*)cmalloc(numAtomTot*sizeof(int));
  for(iAtom=1;iAtom<=numAtomQM;iAtom++){
   cpAtom = cpAtomList[iAtom];
   atomVlncUpAll[cpAtom-1] = cpVlncUp[cpAtom];
   atomVlncDnAll[cpAtom-1] = cpVlncDn[cpAtom];
  }/*endfor*/
  
  fragInfo->numElecUpFragTot = (int*)cmalloc(numFragTot*sizeof(int));
  fragInfo->numElecDnFragTot = (int*)cmalloc(numFragTot*sizeof(int));
  fragInfo->numElecUpFragProc = (int*)cmalloc(numFragProc*sizeof(int));
  fragInfo->numElecDnFragProc = (int*)cmalloc(numFragProc*sizeof(int));
  numElecUpFragTot = fragInfo->numElecUpFragTot;
  numElecDnFragTot = fragInfo->numElecDnFragTot;
  numElecUpFragProc = fragInfo->numElecUpFragProc;
  numElecDnFragProc = fragInfo->numElecDnFragProc;
  // The following part is only correct for molecule fragment
  countAtom = 0;
  for(iMol=0;iMol<numFragTot;iMol++){
    numElecUpFragTot[iMol] = 0;
    numElecDnFragTot[iMol] = 0;
    for(iAtom=0;iAtom<atomNumMol[iMol];iAtom++){
      atomInd = countAtom+iAtom;
      numElecUpFragTot[iMol] += atomVlncUpAll[atomInd];
      numElecDnFragTot[iMol] += atomVlncDnAll[atomInd];
    }
    countAtom += atomNumMol[iMol];
  }
  //end molecule fragment only part
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    numElecUpFragProc[iFrag] = 0;
    numElecDnFragProc[iFrag] = 0;
    for(iAtom=0;iAtom<numAtomFragProc[iFrag];iAtom++){
      atomInd = atomFragMapProc[iFrag][iAtom];
      numElecUpFragProc[iFrag] += atomVlncUpAll[atomInd-1];
      numElecDnFragProc[iFrag] += atomVlncDnAll[atomInd-1];
    }//endfor iAtom
  }//endfor iFrag
  
  /*
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    printf("iFrag %i numElecUpFragProc %i numElecDnFragProc %i\n",iFrag,numElecUpFragProc[iFrag],numElecDnFragProc[iFrag]);
  }
  */

/*======================================================================*/
/* 3) Malloc fragment wave functions	                                */

/*--------------------------------------------------------------------------*/
/*  Partially Malloc Grid Mapping					    */

  if(myidState==0){
    printf("**Partially allocate grid mapping\n");
  }

  fragInfo->numGridFragProc = (int*)cmalloc(numFragProc*sizeof(int));  
  fragInfo->numGridFragTot = (int*)cmalloc(numFragTot*sizeof(int));
  fragInfo->gridMapProc = (int**)cmalloc(numFragProc*sizeof(int*));
  fragInfo->gridMapTot  = (int**)cmalloc(numFragTot*sizeof(int*));
  fragInfo->numGridFragDim = (int**)cmalloc(numFragProc*sizeof(int*));
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    fragInfo->numGridFragDim[iFrag] = (int*)cmalloc(3*sizeof(int));
  }
  // Only used for unit cell fragment

/*======================================================================*/
/* 3) Initialize other things	                                        */

  if(myidState==0){
    printf("**Finish other initializtion\n");
  }

  fragInfo->molSetName = (char *)cmalloc(MAXWORD*sizeof(char));
  if(numProcStates==1)strcpy(fragInfo->molSetName,general_data->filenames.molsetname);
  else{
    if(myidState==0)strcpy(fragInfo->molSetName,general_data->filenames.molsetname);
    Barrier(commStates);
    Bcast(fragInfo->molSetName,MAXWORD,MPI_CHAR,0,commStates);
  }
  fragInfo->cellHmat = (double**)cmalloc(numFragProc*sizeof(double*));
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    fragInfo->cellHmat[iFrag] = (double*)cmalloc(9*sizeof(double));
  }

/*======================================================================*/
/* 4) Initialize skin	                                                */

  if(myidState==0){
    printf("**Initialize skin\n");
  }
  fragInfo->skinAll = (double*)cmalloc(numAtomTot*sizeof(double));
  fragInfo->skinFragBox = (double**)cmalloc(numFragProc*sizeof(double*));
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    fragInfo->skinFragBox[iFrag] = (double*)cmalloc(numAtomFragProc[iFrag]*sizeof(double));
  }
  skinAll = fragInfo->skinAll;
  if(myidState==0){
    atomSkinFile = fragInfo->atomSkinFile;
    printf("%s\n",fragInfo->atomSkinFile);
    fileSkin = NULL;
    fileSkin = fopen(atomSkinFile,"r");
    if(fileSkin==NULL){
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      printf("Skin file doesn't exist!\n");
      printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(0);
    }
    for(iAtom=0;iAtom<numAtomTot;iAtom++){
      fscanf(fileSkin,"%lg",&skinAll[iAtom]);
    }
    fclose(fileSkin);
  }
  if(numProcStates>1){
    Barrier(commStates);
    Bcast(skinAll,numAtomTot,MPI_DOUBLE,0,commStates);
  }
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    for(iAtom=0;iAtom<numAtomFragProc[iFrag];iAtom++){
      fragInfo->skinFragBox[iFrag][iAtom] = skinAll[atomFragMapProc[iFrag][iAtom]-1];
    }//endfor iAtom
  }//endfor iFrag

/*======================================================================*/
/* 4) Initialize skin                                                   */
  
  if(myidState==0){
    printf("**Pass control to parser\n");
  }
  
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void initFragUnitCell(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,CP *cp,
		 int ip_now)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Fragmentation initialization common part for molecular fragments      */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"

  CLATOMS_INFO *clatoms_info = &(class->clatoms_info);
  CLATOMS_POS  *clatoms_pos  = &(class->clatoms_pos[ip_now]);
  ATOMMAPS     *atommaps     = &(class->atommaps);
  CELL         *cell         = &(general_data->cell);
  FOR_SCR      *for_scr      = &(class->for_scr);
  EWD_SCR      *ewd_scr      = &(class->ewd_scr);
  PTENS        *ptens        = &(general_data->ptens);

  CPOPTS        *cpopts           = &(cp->cpopts);
  PSEUDO        *pseudo           = &(cp->pseudo);
  CPCOEFFS_INFO *cpcoeffs_info    = &(cp->cpcoeffs_info);
  COMMUNICATE   *communicate      = &(cp->communicate);
  CPCOEFFS_POS  *cpcoeffs_pos     = &(cp->cpcoeffs_pos[ip_now]);
  STODFTINFO    *stodftInfo       = cp->stodftInfo;
  STODFTCOEFPOS *stodftCoefPos    = cp->stodftCoefPos;
  CPSCR         *cpscr            = &(cp->cpscr);
  NEWTONINFO    *newtonInfo;
  FRAGINFO      *fragInfo	  = stodftInfo->fragInfo;
  PARA_FFT_PKG3D *cp_para_fft_pkg3d_lg = &(cp->cp_para_fft_pkg3d_lg);

  int iChem,iSamp,iCell,iProc,iCoeff,iMol,iFrag,iType,iAtom;
  int count;
  int cpAtom;
  int numFragTot;
  int numFragProc;
  int div,res;
  int fragOpt           = stodftInfo->fragOpt;
  int fragCellOpt       = stodftInfo->fragCellOpt;
  int rhoRealGridNum	= stodftInfo->rhoRealGridNum;
  int skinUCNum;
  int cpLsda		= cpopts->cp_lsda;
  int numMolTot;
  int numMolType        = atommaps->nmol_typ;
  int molIndStart;
  int atomIndTypeStart;
  int countMol,countAtom;
  int molInd,atomInd;
  int numAtomQM		= clatoms_info->nab_initio;
  int numAtomTot	= clatoms_info->natm_tot;
  int myidState		= communicate->myid_state;
  int numProcStates	= communicate->np_states;
  int numFragA,numFragB,numFragC;
  int indx,indy,indz;

  MPI_Comm commStates   = communicate->comm_states;

  int *sysRootInd;
  int *numMolJmolType	   = atommaps->nmol_jmol_typ;
  int *numAtomJmolType	   = atommaps->natm_1mol_jmol_typ;
  int *jatomJmolTypeStart  = atommaps->jatm_jmol_typ_strt;
  int *numMolFragProc;
  int *numAtomFragProc;
  int *numAtomFragAll;
  int *atomIndStart;
  int *atomNumMol;
  int *cpVlncUp		    = clatoms_info->cp_vlnc_up;
  int *cpVlncDn		    = clatoms_info->cp_vlnc_dn;
  int *cpAtomList	    = atommaps->cp_atm_lst;
  int *atomVlncUpAll,*atomVlncDnAll;
  int *numElecUpFragTot,*numElecDnFragTot,*numElecUpFragProc,*numElecDnFragProc;
  int *molTypeMapAll,*molTypeMapFrag,*numMolTypeFrag;
  int *molType;
  int *numGridBox;

  int **molTypeFrag,**molNumTypeFrag;
  int **molFragMapProc;
  int **atomFragMapProc;

  char *atomSkinFile;

  double xRoot,yRoot,zRoot,x,y,z;
  double xDiff,yDiff,zDiff,xTemp,yTemp,zTemp;
  double comSys[3],boxCnt[3];
  double aBig[3],bBig[3],cBig[3];
  double *sysRoot;

  double *xTot = clatoms_pos->x;
  double *yTot = clatoms_pos->y;
  double *zTot = clatoms_pos->z;
  double *comMol,*comMolReduce;
  double *relativeCoord;
  double *hmat  = cell->hmat;
  double *hmati = cell->hmati;
  double *xTotTemp,*yTotTemp,*zTotTemp;

  double *skinAll;
  FILE *fileNumUC;


/*======================================================================*/
/* II) Get mol map in each fragments                                    */

  // Get COM of each molecule, I assume all molecules(atoms) are in the 
  // supper cell
  numMolTot = 0;
  for(iType=1;iType<=numMolType;iType++)numMolTot += numMolJmolType[iType];
  
  molType = (int*)cmalloc(numMolTot*sizeof(int));
  comMolReduce = (double*)cmalloc(numMolTot*3*sizeof(double));
  atomIndStart = (int*)cmalloc(numMolTot*sizeof(int));
  atomNumMol = (int*)cmalloc(numMolTot*sizeof(int));
  fragInfo->sysRootInd = (int*)cmalloc(3*sizeof(int));
  sysRootInd = fragInfo->sysRootInd;
  sysRoot = (double*)cmalloc(3*sizeof(double));
  numGridBox = (int*)cmalloc(3*sizeof(int));
  
  numGridBox[2] = cp_para_fft_pkg3d_lg->nkf3;
  numGridBox[1] = cp_para_fft_pkg3d_lg->nkf2;
  numGridBox[0] = cp_para_fft_pkg3d_lg->nkf1;

  shiftSystem(numMolTot,numAtomTot,numMolType,&molType[0],&atomIndStart[0],&atomNumMol[0],
              &comMolReduce[0],&sysRootInd[0],&sysRoot[0],&numGridBox[0],
              numAtomJmolType,jatomJmolTypeStart,numMolJmolType,xTot,yTot,zTot,cell);
  //for(iMol=0;iMol<100;iMol++)printf("11111111111111111 comMolReduce %lg\n",comMolReduce[iMol]);

/*======================================================================*/
/* II) Read the fragment file and partition molecule to unit cells      */

  partMolUC(&comMolReduce[0],numMolTot,&numGridBox[0],communicate,fragInfo);
  numFragTot = fragInfo->numFragTot;
  skinUCNum = fragInfo->skinUCNum;

  
/*======================================================================*/
/* II) Combine the unit cell to fragment cell and partition molecules   */

  if(myidState==0)printf("**Get mol map\n");
  
  if(skinUCNum>=0){
    mapFragMol(fragInfo,communicate,numMolType,molType,&numGridBox[0]);
  }
  else{//skin=0.5*UC_siez
    mapFragMolHalf(fragInfo,communicate,numMolType,molType,
		    &numGridBox[0],comMolReduce);
  }
  molFragMapProc = fragInfo->molFragMapProc;
  numMolFragProc = fragInfo->numMolFragProc;
  numFragProc = fragInfo->numFragProc;


/*======================================================================*/
/* 2) Get atoms map in each fragments                                   */

  if(myidState==0)printf("**Get atom map\n");
  
  // Get Atom number per fragment
  fragInfo->numAtomFragProc = (int*)cmalloc(numFragProc*sizeof(int));
  numAtomFragProc = fragInfo->numAtomFragProc;
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    numAtomFragProc[iFrag] = 0;
    for(iMol=0;iMol<numMolFragProc[iFrag];iMol++){
      numAtomFragProc[iFrag] += atomNumMol[molFragMapProc[iFrag][iMol]-1];
    }
  }
  // Get Atom Map
  fragInfo->atomFragMapProc = (int**)cmalloc(numFragProc*sizeof(int*));
  atomFragMapProc = fragInfo->atomFragMapProc;
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    atomFragMapProc[iFrag] = (int*)cmalloc(numAtomFragProc[iFrag]*sizeof(int));
  }
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    countAtom = 0;
    for(iMol=0;iMol<numMolFragProc[iFrag];iMol++){
      for(iAtom=0;iAtom<atomNumMol[molFragMapProc[iFrag][iMol]-1];iAtom++){
	atomFragMapProc[iFrag][countAtom+iAtom] = atomIndStart[molFragMapProc[iFrag][iMol]-1]+iAtom;
      }//endfor iAtom
      countAtom += atomNumMol[molFragMapProc[iFrag][iMol]-1];
    }//endfor iMol
  }//endfor iFrag

  /*
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    //printf("iFrag %i numAtomFragProc %i\n",iFrag,numAtomFragProc[iFrag]);
    printf("iFrag %i atomFragMap St %i atomFragMap Ed %i\n",iFrag,atomFragMapProc[iFrag][0],atomFragMapProc[iFrag][numAtomFragProc[iFrag]-1]);
  }
  */

/*======================================================================*/
/* 3) Get electron number in each fragments                             */

  if(myidState==0){
    printf("**Get electron number\n");
  }

  if(numAtomQM!=numAtomTot){
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    printf("I am not interested in calculating 'QMMM'. Let all\n");
    printf("atoms be ab initial!\n");
    printf("@@@@@@@@@@@@@@@@@@@@_ERROR_@@@@@@@@@@@@@@@@@@@@\n");
    fflush(stdout);
    exit(0);
  }
  atomVlncUpAll = (int*)cmalloc(numAtomTot*sizeof(int));
  atomVlncDnAll = (int*)cmalloc(numAtomTot*sizeof(int));
  for(iAtom=1;iAtom<=numAtomQM;iAtom++){
   cpAtom = cpAtomList[iAtom];
   atomVlncUpAll[cpAtom-1] = cpVlncUp[cpAtom];
   atomVlncDnAll[cpAtom-1] = cpVlncDn[cpAtom];
  }/*endfor*/
  
  fragInfo->numElecUpFragTot = (int*)cmalloc(numFragTot*sizeof(int));
  fragInfo->numElecDnFragTot = (int*)cmalloc(numFragTot*sizeof(int));
  fragInfo->numElecUpFragProc = (int*)cmalloc(numFragProc*sizeof(int));
  fragInfo->numElecDnFragProc = (int*)cmalloc(numFragProc*sizeof(int));
  numElecUpFragTot = fragInfo->numElecUpFragTot;
  numElecDnFragTot = fragInfo->numElecDnFragTot;
  numElecUpFragProc = fragInfo->numElecUpFragProc;
  numElecDnFragProc = fragInfo->numElecDnFragProc;
  // The following part is only correct for molecule fragment
  countAtom = 0;
  for(iMol=0;iMol<numFragTot;iMol++){
    numElecUpFragTot[iMol] = 0;
    numElecDnFragTot[iMol] = 0;
    for(iAtom=0;iAtom<atomNumMol[iMol];iAtom++){
      atomInd = countAtom+iAtom;
      numElecUpFragTot[iMol] += atomVlncUpAll[atomInd];
      numElecDnFragTot[iMol] += atomVlncDnAll[atomInd];
    }
    countAtom += atomNumMol[iMol];
  }
  //end molecule fragment only part
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    numElecUpFragProc[iFrag] = 0;
    numElecDnFragProc[iFrag] = 0;
    for(iAtom=0;iAtom<numAtomFragProc[iFrag];iAtom++){
      atomInd = atomFragMapProc[iFrag][iAtom];
      numElecUpFragProc[iFrag] += atomVlncUpAll[atomInd-1];
      numElecDnFragProc[iFrag] += atomVlncDnAll[atomInd-1];
    }//endfor iAtom
  }//endfor iFrag
  
  
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    printf("iFrag %i numElecUpFragProc %i numElecDnFragProc %i\n",iFrag,numElecUpFragProc[iFrag],numElecDnFragProc[iFrag]);
  }
    

/*======================================================================*/
/* 3) Malloc fragment wave functions	                                */

/*--------------------------------------------------------------------------*/
/*  Partially Malloc Grid Mapping					    */

  if(myidState==0){
    printf("**Partially allocate grid mapping\n");
  }

  fragInfo->numGridFragProc = (int*)cmalloc(numFragProc*sizeof(int));  
  fragInfo->numGridFragTot = (int*)cmalloc(numFragTot*sizeof(int));
  fragInfo->gridMapProc = (int**)cmalloc(numFragProc*sizeof(int*));
  fragInfo->gridMapTot  = (int**)cmalloc(numFragTot*sizeof(int*));
  fragInfo->gridMapProcSmall = (int**)cmalloc(numFragProc*sizeof(int*));
  fragInfo->numGridFragDim = (int**)cmalloc(numFragProc*sizeof(int*));
  fragInfo->numGridFragDimSmall = (int**)cmalloc(numFragProc*sizeof(int*));
  fragInfo->numGridFragProcSmall = (int*)cmalloc(numFragProc*sizeof(int));
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    fragInfo->numGridFragDim[iFrag] = (int*)cmalloc(3*sizeof(int));
    fragInfo->numGridFragDimSmall[iFrag] = (int*)cmalloc(3*sizeof(int));
  }
  // Only used for unit cell fragment

/*======================================================================*/
/* 3) Initialize other things	                                        */

  if(myidState==0){
    printf("**Finish other initializtion\n");
  }

  fragInfo->molSetName = (char *)cmalloc(MAXWORD*sizeof(char));
  if(numProcStates==1)strcpy(fragInfo->molSetName,general_data->filenames.molsetname);
  else{
    if(myidState==0)strcpy(fragInfo->molSetName,general_data->filenames.molsetname);
    Barrier(commStates);
    Bcast(fragInfo->molSetName,MAXWORD,MPI_CHAR,0,commStates);
  }
  fragInfo->cellHmat = (double**)cmalloc(numFragProc*sizeof(double*));
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    fragInfo->cellHmat[iFrag] = (double*)cmalloc(9*sizeof(double));
  }

/*======================================================================*/
/* 4) Initialize skin                                                   */
  
  if(myidState==0){
    printf("**Pass control to parser\n");
  }
  
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void initFragEnergy(CP *cpMini,CLASS *classMini,CLASS *class,CP *cp)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Initialize fragment ke/PNL part				         */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  STODFTINFO *stodftInfo	= cp->stodftInfo;
  FRAGINFO *fragInfo		= stodftInfo->fragInfo;
  CLATOMS_INFO *clatomsInfo	= &(class->clatoms_info);
  CLATOMS_INFO *clatomsInfoMini = &(classMini->clatoms_info);
  CPOPTS *cpOpts		= &(cp->cpopts);
  
  int iFrag;
  int cpLsda = cpOpts->cp_lsda;
  int numAtomTot = clatomsInfo->natm_tot;
  int numAtomFrag = clatomsInfoMini->natm_tot;
  int numFragProc = fragInfo->numFragProc;
  int numStateUpMini,numStateDnMini;
  
  fragInfo->vnlFxCor = (double*)cmalloc(numAtomTot*sizeof(double));
  fragInfo->vnlFyCor = (double*)cmalloc(numAtomTot*sizeof(double));
  fragInfo->vnlFzCor = (double*)cmalloc(numAtomTot*sizeof(double));

  fragInfo->wfProjUp = (double**)cmalloc(numFragProc*sizeof(double*));
  fragInfo->wfProjDn = (double**)cmalloc(numFragProc*sizeof(double*));
  fragInfo->keMatrixUp = (double**)cmalloc(numFragProc*sizeof(double*));
  fragInfo->keMatrixDn = (double**)cmalloc(numFragProc*sizeof(double*));
  fragInfo->vnlMatrixUp = (double**)cmalloc(numFragProc*sizeof(double*));
  fragInfo->vnlFxMatrixUp = (double**)cmalloc(numFragProc*sizeof(double*));
  fragInfo->vnlFyMatrixUp = (double**)cmalloc(numFragProc*sizeof(double*));
  fragInfo->vnlFzMatrixUp = (double**)cmalloc(numFragProc*sizeof(double*));
  fragInfo->vnlMatrixDn = (double**)cmalloc(numFragProc*sizeof(double*));
  fragInfo->vnlFxMatrixDn = (double**)cmalloc(numFragProc*sizeof(double*));
  fragInfo->vnlFyMatrixDn = (double**)cmalloc(numFragProc*sizeof(double*));
  fragInfo->vnlFzMatrixDn = (double**)cmalloc(numFragProc*sizeof(double*));
  fragInfo->Fx = (double**)cmalloc(numFragProc*sizeof(double*));
  fragInfo->Fy = (double**)cmalloc(numFragProc*sizeof(double*));
  fragInfo->Fz = (double**)cmalloc(numFragProc*sizeof(double*));

  for(iFrag=0;iFrag<numFragProc;iFrag++){
    numStateUpMini = cpMini[iFrag].cpcoeffs_info.nstate_up_proc;
    fragInfo->wfProjUp[iFrag] = (double*)cmalloc(numStateUpMini*sizeof(double));
    fragInfo->keMatrixUp[iFrag] = (double*)cmalloc(numStateUpMini*numStateUpMini*sizeof(double));
    fragInfo->vnlMatrixUp[iFrag] = (double*)cmalloc(numStateUpMini*numStateUpMini*sizeof(double));
    fragInfo->vnlFxMatrixUp[iFrag] = (double*)cmalloc(numAtomFrag*numStateUpMini*numStateUpMini*sizeof(double));
    fragInfo->vnlFyMatrixUp[iFrag] = (double*)cmalloc(numAtomFrag*numStateUpMini*numStateUpMini*sizeof(double));
    fragInfo->vnlFzMatrixUp[iFrag] = (double*)cmalloc(numAtomFrag*numStateUpMini*numStateUpMini*sizeof(double));
    fragInfo->Fx[iFrag] = (double*)cmalloc(numAtomFrag*sizeof(double));
    fragInfo->Fy[iFrag] = (double*)cmalloc(numAtomFrag*sizeof(double));
    fragInfo->Fz[iFrag] = (double*)cmalloc(numAtomFrag*sizeof(double));

  }
  if(cpLsda==1){
    for(iFrag=0;iFrag<numFragProc;iFrag++){
      numStateDnMini = cpMini[iFrag].cpcoeffs_info.nstate_dn_proc;
      fragInfo->wfProjDn[iFrag] = (double*)cmalloc(numStateDnMini*sizeof(double));
      fragInfo->keMatrixDn[iFrag] = (double*)cmalloc(numStateDnMini*numStateDnMini*sizeof(double));
      fragInfo->vnlMatrixDn[iFrag] = (double*)cmalloc(numStateDnMini*numStateDnMini*sizeof(double));
      fragInfo->vnlFxMatrixDn[iFrag] = (double*)cmalloc(numAtomFrag*numStateDnMini*numStateDnMini*sizeof(double));
      fragInfo->vnlFyMatrixDn[iFrag] = (double*)cmalloc(numAtomFrag*numStateDnMini*numStateDnMini*sizeof(double));
      fragInfo->vnlFzMatrixDn[iFrag] = (double*)cmalloc(numAtomFrag*numStateDnMini*numStateDnMini*sizeof(double));

    }
  }

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void reInitFrag(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,CP *cp,
	       CLASS *classMini,BONDED *bondedMini,GENERAL_DATA *generalDataMini,
	       CP *cpMini)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Every time a new set of coordinates passed into fragment calculation  */
/* I need ot redo */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"

  double geoCnt[3];

  //cleanFFT(generalData,class,cp,generalDataMini,classMini,cpMini,1);

/*======================================================================*/
/* 1) First pass new coords to fragment calculation                     */

  passAtomCoord(general_data,class,cp,generalDataMini,classMini,cpMini,1,geoCnt);

/*======================================================================*/
/* 2) Recalculate fragment box size		                        */

  //initFFTMap(general_data,class,cp,generalDataMini,classMini,cpMini,1,geoCnt);

/*======================================================================*/
/* 2) Reinitialize FFT for fragments box size				*/

  reInitFFT(general_data,class,cp,generalDataMini,classMini,cpMini,1);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void shiftSystem(int numMolTot,int numAtomTot,int numMolType,int *molType,
		 int *atomIndStart,int *atomNumMol,
		 double *comMolReduce,int *sysRootInd,double *sysRoot,int *numGridBox,
		 int *numAtomJmolType,int *jatomJmolTypeStart,int *numMolJmolType,
		 double *xTot,double *yTot,double *zTot,CELL *cell)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Shift the system so that all COM of molecules are in the big box and  */
/* the system is located in box center.					 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  int indx,indy,indz;
  int iType,iMol,iAtom;
  int countMol;
  int atomIndTypeStart;

  double xRoot,yRoot,zRoot,x,y,z;
  double xDiff,yDiff,zDiff,xTemp,yTemp,zTemp;
  double comSys[3] = {0.0};
  double aBig[3],bBig[3],cBig[3];
  double rate;

  double *comMol;
  double *relativeCoord;
  double *hmat  = cell->hmat;
  double *hmati = cell->hmati;
  double *xTotTemp,*yTotTemp,*zTotTemp;

  comMol = (double*)cmalloc(numMolTot*3*sizeof(double));
  xTotTemp = (double*)cmalloc(numAtomTot*sizeof(double));
  yTotTemp = (double*)cmalloc(numAtomTot*sizeof(double));
  zTotTemp = (double*)cmalloc(numAtomTot*sizeof(double));
  for(iMol=0;iMol<3*numMolTot;iMol++)comMol[iMol] = 0.0;

  countMol = 0;
  for(iType=1;iType<=numMolType;iType++){
    atomIndTypeStart = jatomJmolTypeStart[iType];
    relativeCoord = (double*)cmalloc(numAtomJmolType[iType]*3*sizeof(double));
    for(iMol=0;iMol<numMolJmolType[iType];iMol++){
      atomIndStart[countMol+iMol] = atomIndTypeStart+iMol*numAtomJmolType[iType];
      //printf("countMol+iMol %i,atomIndStart %i\n",countMol+iMol,atomIndStart[countMol+iMol]);
      atomNumMol[countMol+iMol] = numAtomJmolType[iType];
      molType[countMol+iMol] = iType;
      // 1. Select a root atom
      xRoot = xTot[atomIndStart[countMol+iMol]];
      yRoot = yTot[atomIndStart[countMol+iMol]];
      zRoot = zTot[atomIndStart[countMol+iMol]];
      for(iAtom=0;iAtom<numAtomJmolType[iType];iAtom++){
	// 2. If the molecule is cut into two pieces due to PBC
	// connect the atoms w.r.t. the root atom
        x = xTot[atomIndStart[countMol+iMol]+iAtom];
        y = yTot[atomIndStart[countMol+iMol]+iAtom];
        z = zTot[atomIndStart[countMol+iMol]+iAtom];
	xDiff = x-xRoot;
	yDiff = y-yRoot;
	zDiff = z-zRoot;
	xTemp = xDiff*hmati[1]+yDiff*hmati[4]+zDiff*hmati[7];
	yTemp = xDiff*hmati[2]+yDiff*hmati[5]+zDiff*hmati[8];
	zTemp = xDiff*hmati[3]+yDiff*hmati[6]+zDiff*hmati[9];
	if(xTemp>0.5)xTemp -= 1.0;
	if(yTemp>0.5)yTemp -= 1.0;
	if(zTemp>0.5)zTemp -= 1.0;
	if(xTemp<-0.5)xTemp += 1.0;
	if(yTemp<-0.5)yTemp += 1.0;
	if(zTemp<-0.5)zTemp += 1.0;
	xDiff = xTemp*hmat[1]+yTemp*hmat[4]+zTemp*hmat[7];
	yDiff = xTemp*hmat[2]+yTemp*hmat[5]+zTemp*hmat[8];
	zDiff = xTemp*hmat[3]+yTemp*hmat[6]+zTemp*hmat[9];
	x = xDiff+xRoot;
	y = yDiff+yRoot;	      
	z = zDiff+zRoot;
	// 3. Calculate the mol COM w.r.t.
        comMol[(countMol+iMol)*3] += x;
        comMol[(countMol+iMol)*3+1] += y;
        comMol[(countMol+iMol)*3+2] += z;
	relativeCoord[iAtom*3] = x;
        relativeCoord[iAtom*3+1] = y;
        relativeCoord[iAtom*3+2] = z;
      }
      comMol[(countMol+iMol)*3] /= (double)numAtomJmolType[iType];
      comMol[(countMol+iMol)*3+1] /= (double)numAtomJmolType[iType];
      comMol[(countMol+iMol)*3+2] /= (double)numAtomJmolType[iType]; 
      // 4. Calculate the relative atom positions w.r.t. mol COM
      for(iAtom=0;iAtom<numAtomJmolType[iType];iAtom++){
	relativeCoord[iAtom*3] -= comMol[(countMol+iMol)*3];
        relativeCoord[iAtom*3+1] -= comMol[(countMol+iMol)*3+1];
        relativeCoord[iAtom*3+2] -= comMol[(countMol+iMol)*3+2];
      }//endfor iAtom
      // 5. Shift the com back to the box
      x = comMol[(countMol+iMol)*3];
      y = comMol[(countMol+iMol)*3+1];
      z = comMol[(countMol+iMol)*3+2];
      xTemp = x*hmati[1]+y*hmati[4]+z*hmati[7];
      yTemp = x*hmati[2]+y*hmati[5]+z*hmati[8];
      zTemp = x*hmati[3]+y*hmati[6]+z*hmati[9];
      if(xTemp>1.0)xTemp -= 1.0;
      if(yTemp>1.0)yTemp -= 1.0;
      if(zTemp>1.0)zTemp -= 1.0;
      if(xTemp<-1.0)xTemp += 1.0;
      if(yTemp<-1.0)yTemp += 1.0;
      if(zTemp<-1.0)zTemp += 1.0;
      comMolReduce[(countMol+iMol)*3] = xTemp;
      comMolReduce[(countMol+iMol)*3+1] = yTemp;
      comMolReduce[(countMol+iMol)*3+2] = zTemp;
      x = xTemp*hmat[1]+yTemp*hmat[4]+zTemp*hmat[7];
      y = xTemp*hmat[2]+yTemp*hmat[5]+zTemp*hmat[8];
      z = xTemp*hmat[3]+yTemp*hmat[6]+zTemp*hmat[9];
      comMol[(countMol+iMol)*3] = x;
      comMol[(countMol+iMol)*3+1] = y;
      comMol[(countMol+iMol)*3+2] = z;
      // 6. Shift the whole molecule w.r.t. updated COM      
      for(iAtom=0;iAtom<numAtomJmolType[iType];iAtom++){
	xTotTemp[atomIndStart[countMol+iMol]-1] = relativeCoord[iAtom*3]+x;
        yTotTemp[atomIndStart[countMol+iMol]-1] = relativeCoord[iAtom*3+1]+y;
        zTotTemp[atomIndStart[countMol+iMol]-1] = relativeCoord[iAtom*3+2]+z;
      }      
    }//endfor iMol
    countMol += numMolJmolType[iType];
    free(relativeCoord);
  }//endfor iType
  // 7. Generate the system COM
  for(iAtom=0;iAtom<numAtomTot;iAtom++){
    comSys[0] += xTotTemp[iAtom];
    comSys[1] += yTotTemp[iAtom];
    comSys[2] += zTotTemp[iAtom];    
  }
  comSys[0] /= (double)numAtomTot;
  comSys[1] /= (double)numAtomTot;
  comSys[2] /= (double)numAtomTot;
  // 8. Transfrom the system COM to reduced box
  x = comSys[0]*hmati[1]+comSys[1]*hmati[4]+comSys[2]*hmati[7];
  y = comSys[0]*hmati[2]+comSys[1]*hmati[5]+comSys[2]*hmati[8];
  z = comSys[0]*hmati[3]+comSys[1]*hmati[6]+comSys[2]*hmati[9];
  //printf("x %lg\n",x);
  // 9. Calculate the difference between com and box center
  // These values are prototype of root point of the fragmentation
  sysRoot[0] = x-0.5;sysRoot[1] = y-0.5;sysRoot[2] = z-0.5;
  // 10. Round the root point to the grid point
  //printf("sysRoot %lg\n",sysRoot[0]);
  rate = sysRoot[0]*numGridBox[0];
  sysRootInd[0] = NINT(rate);
  rate = sysRoot[1]*numGridBox[1];
  sysRootInd[1] = NINT(rate);
  rate = sysRoot[2]*numGridBox[2];
  sysRootInd[2] = NINT(rate);
  sysRoot[0] = ((double)sysRootInd[0])/numGridBox[0];
  sysRoot[1] = ((double)sysRootInd[1])/numGridBox[1];
  sysRoot[2] = ((double)sysRootInd[2])/numGridBox[2];
  // 11. Shift all COMs so that they are in the middle of the box

  for(iMol=0;iMol<numMolTot;iMol++){
    comMolReduce[iMol*3] -= sysRoot[0];
    comMolReduce[iMol*3+1] -= sysRoot[1];
    comMolReduce[iMol*3+2] -= sysRoot[2];
  }
  free(comMol);
  free(xTotTemp);
  free(yTotTemp);
  free(zTotTemp);
  printf("comMolReduce %lg\n",comMolReduce[0]);
/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void partMolUC(double *comMolReduce,int numMolTot,int *numGridBox,
		COMMUNICATE *communicate,FRAGINFO *fragInfo)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Partition all molecules to unit cells				 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
#include "../typ_defs/typ_mask.h"

  FILE *fileNumUC;
  int iFrag,iMol,iCell;
  int indUC;
  int numUCA,numUCB,numUCC;
  int indx,indy,indz;
  int numUCTot;
  int numFragTot;
  int myidState         = communicate->myid_state;
  int numProcStates     = communicate->np_states;
  int *molNumUC;
  int **molIndexUC;

  MPI_Comm commStates   = communicate->comm_states;
    
  double xBin,yBin,zBin;
  double x,y,z;
  double readTemp;

  fragInfo->numUnitCellDim = (int*)cmalloc(3*sizeof(int));
  if(myidState==0){
    // File format: frag=non-overlap part
    // (#Frag) (# UC in skin) (# uc per a) (# uc per b) (# uc per c)
    // (frag 1 start UC index along a) ('' b) ('' c) (frag length/uc length in a) ('' b) ('' c)
    fileNumUC = fopen("NUC-dim","r");
    fscanf(fileNumUC,"%i",&(fragInfo->numFragTot));
    //fscanf(fileNumUC,"%i",&(fragInfo->skinUCNum));
    fscanf(fileNumUC,"%lg",&readTemp);
    if(readTemp!=0.5&&(double)((int)(readTemp))!=readTemp){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("We can only accept skin as integer multiplication of unit cell\n");
      printf("or half of the unit cell. Please set it to integer of 0.5!\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(0);
    }
    if(readTemp==0.5)fragInfo->skinUCNum = -1;
    else fragInfo->skinUCNum = (int)readTemp;
    fscanf(fileNumUC,"%i",&(fragInfo->numUnitCellDim[0]));
    fscanf(fileNumUC,"%i",&(fragInfo->numUnitCellDim[1]));
    fscanf(fileNumUC,"%i",&(fragInfo->numUnitCellDim[2]));
    numFragTot = fragInfo->numFragTot;
    fragInfo->fragStInd = (int*)cmalloc(numFragTot*3*sizeof(int));
    fragInfo->fragLengthInd = (int*)cmalloc(numFragTot*3*sizeof(int));
    for(iFrag=0;iFrag<numFragTot;iFrag++){
      fscanf(fileNumUC,"%i",&(fragInfo->fragStInd[3*iFrag]));
      fscanf(fileNumUC,"%i",&(fragInfo->fragStInd[3*iFrag+1]));
      fscanf(fileNumUC,"%i",&(fragInfo->fragStInd[3*iFrag+2]));
      fscanf(fileNumUC,"%i",&(fragInfo->fragLengthInd[3*iFrag]));
      fscanf(fileNumUC,"%i",&(fragInfo->fragLengthInd[3*iFrag+1]));
      fscanf(fileNumUC,"%i",&(fragInfo->fragLengthInd[3*iFrag+2]));
    }
    fclose(fileNumUC);
  }
  if(numProcStates>1){
    Bcast(&(fragInfo->numUnitCellDim[0]),3,MPI_INT,0,commStates);
    Bcast(&(fragInfo->numFragTot),1,MPI_INT,0,commStates);
    Bcast(&(fragInfo->skinUCNum),1,MPI_INT,0,commStates);
  }
  numFragTot = fragInfo->numFragTot;
  if(myidState!=0){
    fragInfo->fragStInd = (int*)cmalloc(numFragTot*3*sizeof(int));
    fragInfo->fragLengthInd = (int*)cmalloc(numFragTot*3*sizeof(int));
  }
  if(numProcStates>1){
    Bcast(&(fragInfo->fragStInd[0]),3*numFragTot,MPI_INT,0,commStates);
    Bcast(&(fragInfo->fragLengthInd[0]),3*numFragTot,MPI_INT,0,commStates);
  }

  numUCA = fragInfo->numUnitCellDim[0];
  numUCB = fragInfo->numUnitCellDim[1];
  numUCC = fragInfo->numUnitCellDim[2];
  numUCTot = numUCC*numUCB*numUCA;
  xBin = 1.0/(double)(numUCA);
  yBin = 1.0/(double)(numUCB);
  zBin = 1.0/(double)(numUCC);

  if(myidState==0){
    if(numGridBox[0]%numUCA!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("Number of grid point along a direction is not an integer\n");
      printf("multiplication of the number of fragment. Please change\n");
      printf("multiplication of the number of fragment. Please change\n");
      printf("energy cutoff or number of fragment along this direction.\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(0);
    }
    if(numGridBox[1]%numUCB!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("Number of grid point along b direction is not an integer\n");
      printf("multiplication of the number of fragment. Please change\n");
      printf("multiplication of the number of fragment. Please change\n");
      printf("energy cutoff or number of fragment along this direction.\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(0);
    }
    if(numGridBox[2]%numUCC!=0){
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      printf("Number of grid point along c direction is not an integer\n");
      printf("multiplication of the number of fragment. Please change\n");
      printf("multiplication of the number of fragment. Please change\n");
      printf("energy cutoff or number of fragment along this direction.\n");
      printf("@@@@@@@@@@@@@@@@@@@@_error_@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
      fflush(stdout);
      exit(0);
    }
  }

  fragInfo->molNumUC = (int*)cmalloc(numUCTot*sizeof(int));
  fragInfo->molIndexUC = (int**)cmalloc(numUCTot*sizeof(int*));
  molNumUC = fragInfo->molNumUC;
  molIndexUC = fragInfo->molIndexUC;
  for(iCell=0;iCell<numUCTot;iCell++){
    molNumUC[iCell] = 0;
    molIndexUC[iCell] = (int*)cmalloc(100*sizeof(int));
  }
  for(iMol=0;iMol<numMolTot;iMol++){
    x = comMolReduce[iMol*3];
    y = comMolReduce[iMol*3+1];
    z = comMolReduce[iMol*3+2];
    indx = (int)(x/xBin);
    indy = (int)(y/yBin);
    indz = (int)(z/zBin);
    printf("indx %i indy %i indz %i\n",indx,indy,indz);
    indUC = indz*numUCB*numUCA+indy*numUCA+indx;
    //printf("indx %i indy %i indz %i indUC %i\n",indx,indy,indz,indUC);
    molNumUC[indUC] += 1;
    if(molNumUC[indUC]%100==0){
      molIndexUC[indUC] = (int*)realloc(molIndexUC[indUC],
				            (molNumUC[indUC]+100)*sizeof(double));
    }
    molIndexUC[indUC][molNumUC[indUC]-1] = iMol+1;
  }
  if(numProcStates>1)Barrier(commStates);
  if(myidState==0){
    for(iCell=0;iCell<numUCTot;iCell++){
      printf("iCell %i molNumUC %i\n",iCell,molNumUC[iCell]);
    }
  }
  
  /*
  for(iCell=0;iCell<27;iCell++){
    for(iMol=0;iMol<8;iMol++){
      printf("iCell %i iMol %i molIndexUC %i\n",iCell,iMol,molIndexUC[iCell][iMol]);
    }
  }
  fflush(stdout);
  exit(0);
  */

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void mapFragMol(FRAGINFO *fragInfo,COMMUNICATE *communicate,
		int numMolType,int *molType,int *numGridBox)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Partition all molecules to unit cells                                 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  int iFrag,iuc,juc,kuc,iMol;
  int iGrid,jGrid,kGrid;
  int div,res;
  int numFragTot = fragInfo->numFragTot;
  int numFragProc;
  int numMolFragCent;
  int startFragInd,fragIndNow;
  int ucIndA,ucIndB,ucIndC;
  int ucInd,molInd,molTypeInd;
  int myidState         = communicate->myid_state;
  int numProcStates     = communicate->np_states;
  int skinUCNum		= fragInfo->skinUCNum;

  int *molTypeFragTemp;
  int *molIndFragTemp;
  int molNumFragTemp;
  int molTypeNumFragTemp;
  int inFragFlag;
  int *fragInd;
  int *numUnitCellDim = fragInfo->numUnitCellDim;
  int *molNumUC = fragInfo->molNumUC;
  int *fragStInd = fragInfo->fragStInd;
  int *fragLengthInd = fragInfo->fragLengthInd;
  int *sysRootInd = fragInfo->sysRootInd;
  int *molIndFragCent;
  int numGridUCDim[3];

  int **molIndexUC = fragInfo->molIndexUC;

  div = numFragTot/numProcStates;
  res = numFragTot%numProcStates;
  if(myidState<res){
    fragInfo->numFragProc = div+1;
    startFragInd = (div+1)*myidState;
  }
  else{
    fragInfo->numFragProc = div;
    startFragInd = (div+1)*res+div*(myidState-res);
  }
  numFragProc = fragInfo->numFragProc;
  fragInd = (int*)cmalloc(numFragProc*sizeof(int));
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    fragInd[iFrag] = startFragInd+iFrag;
  }
  fragInfo->numMolTypeFrag = (int*)cmalloc(numFragProc*sizeof(int));
  fragInfo->numMolFragProc = (int*)cmalloc(numFragProc*sizeof(int));
  fragInfo->molTypeFrag = (int**)cmalloc(numFragProc*sizeof(int*));
  fragInfo->molNumTypeFrag = (int**)cmalloc(numFragProc*sizeof(int*));
  fragInfo->molFragMapProc = (int**)cmalloc(numFragProc*sizeof(int*));
  fragInfo->fragRootInd = (int*)cmalloc(3*numFragProc*sizeof(int));
  numGridUCDim[0] = numGridBox[0]/numUnitCellDim[0];
  numGridUCDim[1] = numGridBox[1]/numUnitCellDim[1];
  numGridUCDim[2] = numGridBox[2]/numUnitCellDim[2];

  for(iFrag=0;iFrag<numFragProc;iFrag++){
    fragIndNow = fragInd[iFrag];
    iuc = fragStInd[fragIndNow*3]-skinUCNum;
    juc = fragStInd[fragIndNow*3+1]-skinUCNum;
    kuc = fragStInd[fragIndNow*3+2]-skinUCNum;
    iGrid = iuc*numGridUCDim[0]+sysRootInd[0];
    jGrid = juc*numGridUCDim[1]+sysRootInd[1];
    kGrid = kuc*numGridUCDim[2]+sysRootInd[2];
    if(iGrid<0)iGrid += numGridBox[0];
    if(jGrid<0)jGrid += numGridBox[1];
    if(kGrid<0)kGrid += numGridBox[2];
    fragInfo->fragRootInd[iFrag*3] = iGrid;
    fragInfo->fragRootInd[iFrag*3+1] = jGrid;
    fragInfo->fragRootInd[iFrag*3+2] = kGrid;
  }

  // Start build Fragment
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    fragIndNow = fragInd[iFrag];
    molNumFragTemp = 0;
    molTypeNumFragTemp = 0;
    molTypeFragTemp = (int*)cmalloc(numMolType*sizeof(int));
    molIndFragTemp = NULL;
    molIndFragCent = NULL;
    inFragFlag = 0;
    numMolFragCent = 0;
    for(iuc=-skinUCNum;iuc<fragLengthInd[fragIndNow*3]+skinUCNum;iuc++){
      if(iuc>0&&iuc<fragLengthInd[fragIndNow*3])inFragFlag += 1;
      ucIndA = iuc+fragStInd[fragIndNow*3];
      if(ucIndA<0)ucIndA += numUnitCellDim[0];
      if(ucIndA>=numUnitCellDim[0])ucIndA -= numUnitCellDim[0];
      for(juc=-skinUCNum;juc<fragLengthInd[fragIndNow*3+1]+skinUCNum;juc++){
	if(juc>0&&juc<fragLengthInd[fragIndNow*3+1])inFragFlag += 1;
	ucIndB = juc+fragStInd[fragIndNow*3+1];
	if(ucIndB<0)ucIndB += numUnitCellDim[1];
	if(ucIndB>=numUnitCellDim[1])ucIndB -= numUnitCellDim[1];
	for(kuc=-skinUCNum;kuc<fragLengthInd[fragIndNow*3+2]+skinUCNum;kuc++){
	  if(kuc>0&&kuc<fragLengthInd[fragIndNow*3+2])inFragFlag += 1;
	  if(inFragFlag==3)numMolFragCent += molNumUC[ucInd];
	  ucIndC = kuc+fragStInd[fragIndNow*3+2];
	  if(ucIndC<0)ucIndC += numUnitCellDim[2];
          if(ucIndC>=numUnitCellDim[2])ucIndC -= numUnitCellDim[2];
	  ucInd = ucIndC*numUnitCellDim[1]*numUnitCellDim[0]+ucIndB*numUnitCellDim[0]+ucIndA;
	  molNumFragTemp += molNumUC[ucInd];
	  //printf("iuc %i juc %i kuc %i\n",iuc,juc,kuc);
	  //printf("ucInd %i molNumUC %i\n",ucInd,molNumUC[ucInd]);
	  molIndFragTemp = (int*)crealloc(molIndFragTemp,molNumFragTemp*sizeof(int));
	  molIndFragCent = (int*)crealloc(molIndFragCent,numMolFragCent*sizeof(int));
	  for(iMol=0;iMol<molNumUC[ucInd];iMol++){
	    molInd = molIndexUC[ucInd][iMol];
	    molTypeInd = molType[molInd-1];
	    //printf("molTypeInd %i\n",molTypeInd);
	    if(checkInList(molTypeInd,molTypeFragTemp,molTypeNumFragTemp)==0){
	      molTypeNumFragTemp += 1;
	      molTypeFragTemp[molTypeNumFragTemp-1] = molTypeInd;
	    }
	    //printf("iuc %i juc %i kuc %i ucInd %i iMol %i molInd %i\n",iuc,juc,kuc,ucInd,iMol,molInd);
	    molIndFragTemp[molNumFragTemp-molNumUC[ucInd]+iMol] = molInd;
	    if(inFragFlag==3)molIndFragCent[numMolFragCent-molNumUC[ucInd]+iMol] = molInd;
	  }//endfor iMol
	}//endfor kuc
      }//endfor juc
    }//endfor iuc
    // Reorder the molecule
    fragInfo->numMolFragProc[iFrag] = molNumFragTemp;
    reorderMol(fragInfo,molTypeNumFragTemp,molNumFragTemp,
	       molType,molTypeFragTemp,molIndFragTemp,iFrag);
    /*
    printf("numMolTypeFrag %i\n",fragInfo->numMolTypeFrag[iFrag]);
    int iType;
    for(iType=0;iType<fragInfo->numMolTypeFrag[iFrag];iType++){
      printf("iType %i molTypeFrag %i molNumTypeFrag %i molTypeFragTemp %i\n",iType,fragInfo->molTypeFrag[iFrag][iType],fragInfo->molNumTypeFrag[iFrag][iType],molTypeFragTemp[iType]);
    }
    printf("numMolFragProc %i\n",fragInfo->numMolFragProc[iFrag]);
    for(iMol=0;iMol<fragInfo->numMolFragProc[iFrag];iMol++){
      printf("iMol %i molFragMapProc %i\n",iMol,fragInfo->molFragMapProc[iFrag][iMol]);
    }
    */
    /*
    for(iMol=0;iMol<216;iMol++){
      printf("iMol %i molFragMapProc %i %i\n",iMol,fragInfo->molFragMapProc[iFrag][iMol],molIndFragTemp[iMol]);
    }
    exit(0);
    */
  }//endfor iFrag

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void reorderMol(FRAGINFO *fragInfo,int molTypeNumFragTemp,int molNumFragTemp,
                int *molType,int *molTypeFragTemp,int *molIndFragTemp,int iFrag)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Reorder the molecules so that the order are the same as molecular type*/
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  int *molTypeFrag; 
  int *numMolTypeFrag = fragInfo->numMolTypeFrag;
  int *molNumTypeFrag; 
  int *numMolFragProc = fragInfo->numMolFragProc;
  int *molFragMapProc;
  int iType,iMol;
  int countMol;

  numMolTypeFrag[iFrag] = molTypeNumFragTemp;
  numMolFragProc[iFrag] = molNumFragTemp;

  //printf("iFrag %i molNumFragTemp %i\n",iFrag,molNumFragTemp);
  fragInfo->molTypeFrag[iFrag] = (int*)cmalloc(molTypeNumFragTemp*sizeof(int));
  fragInfo->molNumTypeFrag[iFrag] = (int*)cmalloc(molTypeNumFragTemp*sizeof(int));
  fragInfo->molFragMapProc[iFrag] = (int*)cmalloc(molNumFragTemp*sizeof(int));
  molTypeFrag = fragInfo->molTypeFrag[iFrag];
  molNumTypeFrag = fragInfo->molNumTypeFrag[iFrag];
  molFragMapProc = fragInfo->molFragMapProc[iFrag];
  
  countMol = 0;
  for(iType=0;iType<molTypeNumFragTemp;iType++){
    molTypeFrag[iType] = molTypeFragTemp[iType];
    //printf("iType %i molTypeFrag %i molTypeFragTemp %i\n",iType,molTypeFrag[iType],molTypeFragTemp[iType]);
    molNumTypeFrag[iType] = 0;
    for(iMol=0;iMol<molNumFragTemp;iMol++){
      if(molType[molIndFragTemp[iMol]-1]==molTypeFragTemp[iType]){
	molFragMapProc[countMol] = molIndFragTemp[iMol];	
	molNumTypeFrag[iType] += 1;
	//printf("iType %i iMol %i countMol %i,molNumTypeFrag %i molFragMapProc %i molTypeFrag %i\n",iType,iMol,countMol,molNumTypeFrag[iType],molFragMapProc[countMol],molTypeFrag[0]);
        countMol += 1;
      }
    }//endfor iMol
  }//endfor iType
  //printf("local %i molTypeFrag %i global %i\n",molTypeFragTemp[0],molTypeFrag[0],fragInfo->molTypeFrag[iFrag][0]);
  //exit(0);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
int checkInList(int indNew,int *typeList,int numTypeList)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Check if a mol type is in list					 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  int i,inFlag;
  inFlag = 0;
  for(i=0;i<numTypeList;i++){
    if(typeList[i]==indNew){
      inFlag = 1;
      break;
    }
  }

  return inFlag;

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void mapFragMolHalf(FRAGINFO *fragInfo,COMMUNICATE *communicate,
	int numMolType,int *molType,int *numGridBox,double *comMolReduce)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Partition all molecules to unit cells                                 */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
  int iFrag,iuc,juc,kuc,iMol;
  int iGrid,jGrid,kGrid;
  int iGridSmall,jGridSmall,kGridSmall;
  int div,res;
  int numFragTot = fragInfo->numFragTot;
  int numFragProc;
  int numMolFragCent;
  int startFragInd,fragIndNow;
  int ucIndA,ucIndB,ucIndC;
  int ucInd,molInd,molTypeInd;
  int myidState         = communicate->myid_state;
  int numProcStates     = communicate->np_states;
  int skinUCNum	    = fragInfo->skinUCNum;
  int countMol;

  int *molTypeFragTemp;
  int *molIndFragTemp;
  int molNumFragTemp;
  int molTypeNumFragTemp;
  int inFragFlag;
  int *fragInd;
  int *numUnitCellDim = fragInfo->numUnitCellDim;
  int *molNumUC = fragInfo->molNumUC;
  int *fragStInd = fragInfo->fragStInd;
  int *fragLengthInd = fragInfo->fragLengthInd;
  int *sysRootInd = fragInfo->sysRootInd;
  int *molIndFragCent;
  int *fragRootInd;
  int *gridShift;
  int numGridUCDim[3];

  int **molIndexUC = fragInfo->molIndexUC;

  double iGridF,jGridF,kGridF;
  double xBase,yBase,zBase,x,y,z;
  double aFragLength,bFragLength,cFragLength;

  div = numFragTot/numProcStates;
  res = numFragTot%numProcStates;
  if(myidState<res){
    fragInfo->numFragProc = div+1;
    startFragInd = (div+1)*myidState;
  }
  else{
    fragInfo->numFragProc = div;
    startFragInd = (div+1)*res+div*(myidState-res);
  }
  numFragProc = fragInfo->numFragProc;
  fragInd = (int*)cmalloc(numFragProc*sizeof(int));
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    fragInd[iFrag] = startFragInd+iFrag;
  }
  fragInfo->numMolTypeFrag = (int*)cmalloc(numFragProc*sizeof(int));
  fragInfo->numMolFragProc = (int*)cmalloc(numFragProc*sizeof(int));
  fragInfo->molTypeFrag = (int**)cmalloc(numFragProc*sizeof(int*));
  fragInfo->molNumTypeFrag = (int**)cmalloc(numFragProc*sizeof(int*));
  fragInfo->molFragMapProc = (int**)cmalloc(numFragProc*sizeof(int*));
  fragInfo->fragRootInd = (int*)cmalloc(3*numFragProc*sizeof(int));
  fragInfo->gridShift = (int*)cmalloc(3*numFragProc*sizeof(int));
  fragRootInd = fragInfo->fragRootInd;
  gridShift = fragInfo->gridShift;
  numGridUCDim[0] = numGridBox[0]/numUnitCellDim[0];
  numGridUCDim[1] = numGridBox[1]/numUnitCellDim[1];
  numGridUCDim[2] = numGridBox[2]/numUnitCellDim[2];

  for(iFrag=0;iFrag<numFragProc;iFrag++){
    fragIndNow = fragInd[iFrag];
    iuc = fragStInd[fragIndNow*3];
    juc = fragStInd[fragIndNow*3+1];
    kuc = fragStInd[fragIndNow*3+2];

    iGridSmall = iuc*numGridUCDim[0]+sysRootInd[0];
    jGridSmall = juc*numGridUCDim[1]+sysRootInd[1];
    kGridSmall = kuc*numGridUCDim[2]+sysRootInd[2];

    iGridF = (iuc-0.5)*numGridUCDim[0]+sysRootInd[0];
    jGridF = (juc-0.5)*numGridUCDim[1]+sysRootInd[1];
    kGridF = (kuc-0.5)*numGridUCDim[2]+sysRootInd[2];
    iGrid = NINT(iGrid);
    jGrid = NINT(jGrid);
    kGrid = NINT(kGrid);
    gridShift[3*iFrag] = iGridSmall-iGrid;
    gridShift[3*iFrag+1] = jGridSmall-jGrid;
    gridShift[3*iFrag+2] = kGridSmall-kGrid;
    if(iGrid<0)iGrid += numGridBox[0];
    if(jGrid<0)jGrid += numGridBox[1];
    if(kGrid<0)kGrid += numGridBox[2];
    fragRootInd[iFrag*3] = iGrid;
    fragRootInd[iFrag*3+1] = jGrid;
    fragRootInd[iFrag*3+2] = kGrid;
  }

  // Start build Fragment
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    fragIndNow = fragInd[iFrag];
    molNumFragTemp = 0;
    molTypeNumFragTemp = 0;
    molTypeFragTemp = (int*)cmalloc(numMolType*sizeof(int));
    molIndFragTemp = NULL;
    molIndFragCent = NULL;
    inFragFlag = 0;
    numMolFragCent = 0;
    
    xBase = ((double)fragRootInd[iFrag*3])/numGridBox[0];
    yBase = ((double)fragRootInd[iFrag*3+1])/numGridBox[1];
    zBase = ((double)fragRootInd[iFrag*3+2])/numGridBox[2];
    aFragLength = 1.0/numUnitCellDim[0]*(fragLengthInd[fragInd[iFrag]*3]+1.0);
    bFragLength = 1.0/numUnitCellDim[1]*(fragLengthInd[fragInd[iFrag]*3+1]+1.0);
    cFragLength = 1.0/numUnitCellDim[2]*(fragLengthInd[fragInd[iFrag]*3+2]+1.0);
    countMol = 0;
    for(iuc=-1;iuc<fragLengthInd[fragIndNow*3]+1;iuc++){
      if(iuc>0&&iuc<fragLengthInd[fragIndNow*3])inFragFlag += 1;
      ucIndA = iuc+fragStInd[fragIndNow*3];
      if(ucIndA<0)ucIndA += numUnitCellDim[0];
      if(ucIndA>=numUnitCellDim[0])ucIndA -= numUnitCellDim[0];
      for(juc=-1;juc<fragLengthInd[fragIndNow*3+1]+1;juc++){
	if(juc>0&&juc<fragLengthInd[fragIndNow*3+1])inFragFlag += 1;
	ucIndB = juc+fragStInd[fragIndNow*3+1];
	if(ucIndB<0)ucIndB += numUnitCellDim[1];
	if(ucIndB>=numUnitCellDim[1])ucIndB -= numUnitCellDim[1];
	for(kuc=-1;kuc<fragLengthInd[fragIndNow*3+2]+1;kuc++){
	  if(kuc>0&&kuc<fragLengthInd[fragIndNow*3+2])inFragFlag += 1;
	  if(inFragFlag==3)numMolFragCent += molNumUC[ucInd];
	  ucIndC = kuc+fragStInd[fragIndNow*3+2];
	  if(ucIndC<0)ucIndC += numUnitCellDim[2];
	  if(ucIndC>=numUnitCellDim[2])ucIndC -= numUnitCellDim[2];
	  ucInd = ucIndC*numUnitCellDim[1]*numUnitCellDim[0]+ucIndB*numUnitCellDim[0]+ucIndA;
	  molNumFragTemp += molNumUC[ucInd];
	  //printf("iuc %i juc %i kuc %i\n",iuc,juc,kuc);
	  //printf("ucInd %i molNumUC %i\n",ucInd,molNumUC[ucInd]);
	  molIndFragTemp = (int*)crealloc(molIndFragTemp,molNumFragTemp*sizeof(int));
	  molIndFragCent = (int*)crealloc(molIndFragCent,numMolFragCent*sizeof(int));
	  for(iMol=0;iMol<molNumUC[ucInd];iMol++){
	    molInd = molIndexUC[ucInd][iMol];
	    molTypeInd = molType[molInd-1];
	    x = comMolReduce[3*(molInd-1)]-xBase;
	    y = comMolReduce[3*(molInd-1)+1]-yBase;
	    z = comMolReduce[3*(molInd-1)+2]-zBase;
	    if(x<0.0)x += 1.0;
	    if(y<0.0)y += 1.0;
	    if(z<0.0)z += 1.0;
	    if(x<aFragLength&&y<bFragLength&&z<cFragLength){
	      //printf("molTypeInd %i\n",molTypeInd);
	      if(checkInList(molTypeInd,molTypeFragTemp,molTypeNumFragTemp)==0){
		molTypeNumFragTemp += 1;
		molTypeFragTemp[molTypeNumFragTemp-1] = molTypeInd;
	      }
	      //printf("iuc %i juc %i kuc %i ucInd %i iMol %i molInd %i\n",iuc,juc,kuc,ucInd,iMol,molInd);
	      molIndFragTemp[countMol] = molInd;
	      if(inFragFlag==3)molIndFragCent[countMol] = molInd;
	      countMol += 1;
	    }//endif
	  }//endfor iMol
	}//endfor kuc
      }//endfor juc
    }//endfor iuc
    // Reorder the molecule
    fragInfo->numMolFragProc[iFrag] = countMol;
    reorderMol(fragInfo,molTypeNumFragTemp,molNumFragTemp,
           molType,molTypeFragTemp,molIndFragTemp,iFrag);
    free(molTypeFragTemp);
    free(molIndFragTemp);
  }//endfor iFrag

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/

