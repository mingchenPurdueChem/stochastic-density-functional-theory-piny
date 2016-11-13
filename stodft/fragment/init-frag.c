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
  STODFTINFO    *stodftInfo       = cp->stodftInfo;
  FRAGINFO      *fragInfo;

  int fragOpt           = stodftInfo->fragOpt;
  int numFragProc;
  int iFrag;

  stodftInfo->fragInfo = (FRAGINFO*)cmalloc(sizeof(FRAGINFO));
  fragInfo = stodftInfo->fragInfo;
  
  switch(fragOpt){ 
    case 1:
      initFragMol(class,bonded,general_data,cp,ip_now);
      break;
    case 4:
      initFragUnitCell(class,bonded,general_data,cp,ip_now);
  }


  numFragProc = fragInfo->numFragProc;
  *classMiniPoint = (CLASS*)cmalloc(numFragProc*sizeof(CLASS));
  *bondedMiniPoint = (BONDED*)cmalloc(numFragProc*sizeof(BONDED));
  *generalDataMiniPoint = (GENERAL_DATA*)cmalloc(numFragProc*sizeof(GENERAL_DATA));
  *analysisMiniPoint = (ANALYSIS*)cmalloc(numFragProc*sizeof(ANALYSIS));
  *cpMiniPoint = (CP*)cmalloc(numFragProc*sizeof(CP));
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    parseFrag(class,bonded,general_data,cp,analysis,&(*classMiniPoint[iFrag]),
	      &(*bondedMiniPoint[iFrag]),&(*generalDataMiniPoint[iFrag]),
	      &(*cpMiniPoint[iFrag]),&(*analysisMiniPoint[iFrag]));
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

/*======================================================================*/
/* I) Get number of fragments per processor                             */

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
      countAtom += atomNumMol[molFragMapProc[iFrag][iMol]];
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
   atomVlncUpAll[cpAtom] = cpVlncUp[cpAtom];
   atomVlncDnAll[cpAtom] = cpVlncDn[cpAtom];
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
      numElecUpFragProc[iFrag] += atomVlncUpAll[atomInd];
      numElecDnFragProc[iFrag] += atomVlncDnAll[atomInd];
    }//endfor iAtom
  }//endfor iFrag
  
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    printf("iFrag %i numElecUpFragProc %i numElecDnFragProc %i\n",iFrag,numElecUpFragProc[iFrag],numElecDnFragProc[iFrag]);
  }

/*======================================================================*/
/* 3) Malloc fragment wave functions	                                */
    
  fragInfo->rhoFragSum = (double*)cmalloc(rhoRealGridNum*sizeof(double));
  fragInfo->coefUpFragProc = (double***)cmalloc(numFragProc*sizeof(double**));
  fragInfo->coefUpFragTot = (double***)cmalloc(numFragTot*sizeof(double**));
  
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    fragInfo->coefUpFragProc[iFrag] = (double**)cmalloc(numElecUpFragProc[iFrag]*sizeof(double*));
  }
  for(iFrag=0;iFrag<numFragTot;iFrag++){
    fragInfo->coefUpFragTot[iFrag] = (double**)cmalloc(numElecUpFragTot[iFrag]*sizeof(double*));
  }
  if(cpLsda==1){
    fragInfo->coefDnFragProc = (double***)cmalloc(numFragProc*sizeof(double**));
    fragInfo->coefDnFragTot = (double***)cmalloc(numFragTot*sizeof(double**));
    for(iFrag=0;iFrag<numFragProc;iFrag++){
      fragInfo->coefDnFragProc[iFrag] = (double**)cmalloc(numElecUpFragProc[iFrag]*sizeof(double*));
    }
    for(iFrag=0;iFrag<numFragTot;iFrag++){
      fragInfo->coefDnFragTot[iFrag] = (double**)cmalloc(numElecUpFragTot[iFrag]*sizeof(double*));
    }
  }

/*--------------------------------------------------------------------------*/
/*  Partially Malloc Grid Mapping					    */

  fragInfo->numGridFragProc = (int*)cmalloc(numFragProc*sizeof(int));  
  fragInfo->numGridFragTot = (int*)cmalloc(numFragTot*sizeof(int));
  fragInfo->gridMapProc = (int**)cmalloc(numFragProc*sizeof(int*));
  fragInfo->gridMapTot  = (int**)cmalloc(numFragTot*sizeof(int*));
  fragInfo->numGridFragDim = (int**)cmalloc(numFragProc*sizeof(int*));
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    fragInfo->numGridFragDim[iFrag] = (int*)cmalloc(3*sizeof(int));
  }

/*======================================================================*/
/* 3) Initialize other things	                                        */
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

  FILE *fileSkin;
  fragInfo->skinAll = (double*)cmalloc(numAtomTot*sizeof(double));
  fragInfo->skinFragBox = (double**)cmalloc(numFragProc*sizeof(double*));
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    fragInfo->skinFragBox[iFrag] = (double*)cmalloc(numAtomFragProc[iFrag]*sizeof(double));
  }
  double *skinAll = fragInfo->skinAll;
  if(myidState==0){
    fileSkin = fopen("atomskin","r");
    for(iAtom=0;iAtom<numAtomTot;iAtom++){
      fscanf(fileSkin,"%lg",&skinAll[iAtom]);
    }
  }
  Barrier(commStates);
  Bcast(skinAll,numAtomTot,MPI_DOUBLE,0,commStates);
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    for(iAtom=0;iAtom<numAtomFragProc[iFrag];iAtom++){
      fragInfo->skinFragBox[iFrag][iAtom] = skinAll[atomFragMapProc[iFrag][iAtom]-1];
    }//endfor iAtom
  }//endfor iFrag
  
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

  initFFTMap(general_data,class,cp,generalDataMini,classMini,cpMini,1,geoCnt);

/*======================================================================*/
/* 2) Reinitialize FFT for fragments box size				*/

  reInitFFT(general_data,class,cp,generalDataMini,classMini,cpMini,1);

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


