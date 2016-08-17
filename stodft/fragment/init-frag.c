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
#include "../proto_defs/proto_energy_cpcon_entry.h"
#include "../proto_defs/proto_energy_cpcon_local.h"
#include "../proto_defs/proto_energy_cp_local.h"
#include "../proto_defs/proto_friend_lib_entry.h"
#include "../proto_defs/proto_math.h"
#include "../proto_defs/proto_communicate_wrappers.h"
#include "../proto_defs/proto_stodft_local.h"

#define TIME_CP_OFF
/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void initFrag(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,CP *cp)
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

  stodftInfo->fragInfo = (FRAGINFO*)cmalloc(sizeof(FRAGINFO));
  fraginfo = stodftInfo->fragInfo;
  
  switch(fragOpt){ 
    case 1:
      initFragMol(class,bonded,general_data,cp);
      break;
  }

/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


/*==========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*==========================================================================*/
void initFragMol(CLASS *class,BONDED *bonded,GENERAL_DATA *general_data,CP *cp)
/*==========================================================================*/
/*         Begin Routine                                                    */
   {/*Begin Routine*/
/*************************************************************************/
/* Fragmentation initialization common part for molecular fragments      */
/*************************************************************************/
/*=======================================================================*/
/*         Local Variable declarations                                   */
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
  int numMolTot;
  int numMolType        = atommaps->nmol_typ;
  int molIndStart;
  int atomIndTypeStart;
  int countMol,countAtom;
  int molInd,atomInd;
  int numAtomQM		= clatoms_info->nab_initio;
  int numAtomTot	= clatoms_info->natm_tot;
  MPI_Comm commStates   = communicate->comm_states;

  int *numMolJmolType	   = atommaps->nmol_jmol_typ;
  int *numAtomJmolType	   = atommaps->natm_1mol_jmol_typ;
  int *jatomJmolTypeStart  = atommaps->jatm_jmol_typ_strt
  int *numMolFragProc;
  int *numAtomFragProc;
  int *numAtomFragAll;
  int *atomIndStart;
  int *atomNumMol;
  int *cpVlncUp		    = clatoms_info->cp_vlnc_up;
  int *cpVlncDn		    = clatoms_info->cp_vlnc_dn;
  int *numElecUpFragTot;
  int *cpAtomList	    = atommaps->cp_atm_lst;
  int *atomVlncUpAll,*atomVlncDnAll;
  int *numElecUpFragTot,*numElecDnFragTot,*numElecUpFragProc,*numElecDnFragProc;


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

  fragInfo->numMolFragProc = (int*)cmalloc(numFragProc*sizeof(int));
  numMolFragProc = fragInfo->numMolFragProc;
  for(iFrag=0;iFrag<numFragProc;iFrag++)numMolFragProc[iFrag] = 1;
  
  fragInfo->molFragMapProc = (int**)cmalloc(numFragProc*sizeof(int*));
  molFragMapProc = fragIndo->molFragMapProc;
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    molFragMapProc[iFrag] = (int*)cmalloc(numMolFragProc[iFrag]*sizeof(int));
  }

  if(myidState<res)molIndStart = myidState*(div+1);
  else molIndStart = res*(div+1)+(myidState-res)*div;
  for(iFrag=0;iFrag<numFragProc;iFrag++){
    molFragMapProc[iFrag][0] = iFrag+molIndStart+1;
  }

/*======================================================================*/
/* 2) Get atoms map in each fragments                                   */

  atomIndStart = (int*)cmalloc(numMolTot*sizeof(int));
  atomNumMol = (int*)cmalloc(numMolTot*sizeof(int));
  // Get starting atom index for each molecule
  countMol = 0
  for(iType=1;iType<=numMolType;iType++){
    atomIndTypeStart = jatomJmolTypeStart[iType];
    for(iMol=0;iMol<numMolJmolType[iType];iMol++){
      atomIndStart[countMol+iMol] = atomIndTypeStart+iMol*numAtomJmolType[iType];
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
      numAtomFragProc[iFrag] += atomNumMol[molFragMapProc[iFrag][iMol]];
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
      for(iAtom=0;iAtom<atomNumMol[molFragMapProc[iFrag][iMol]];iAtom++){
	atomFragMapProc[iFrag][countAtom+iAtom] = atomIndStart[molFragMapProc[iFrag][iMol]]+iAtom
      }//endfor iAtom
      countAtom += atomNumMol[molFragMapProc[iFrag][iMol]]
    }//endfor iMol
  }//endfor iFrag

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
  for(iMol=0;iMol<numFragtot;iMol++){
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


/*==========================================================================*/
}/*end Routine*/
/*==========================================================================*/


