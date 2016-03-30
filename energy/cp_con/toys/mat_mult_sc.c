#include <stdio.h>
#include <stdlib.h>

#define DEBUG
#define MULTIPLY_OFF


/*========================================================================*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/*========================================================================*/

main()

/*========================================================================*/
/* Begin Routine                                                          */
  {/* begin routine */
/*========================================================================*/
/*          Local Variables         */
  int ind;
  int nm,nl_tot,nl_tot2;
  int joff,jl,knd,jnd;
  int ioff,il,im;
  int done;
  double frac,whole;
  double *a,*c,*b;

/*========================================================================*/
/* I) Malloc memory                                                       */

   nl_tot = 17;
   nm = 7;
   nl_tot2 = nl_tot*nl_tot;

   a             = (double *) malloc(nm*nl_tot*sizeof(double))-1;
   b             = (double *) malloc(nm*nl_tot*sizeof(double))-1;
   c             = (double *) malloc(nl_tot2*sizeof(double))-1;


/*========================================================================*/
/* II) Create model data                                                */

  ioff=0;
  for(il=1;il<=nl_tot;il++){
    for(im=1;im<=nm;im++){
      frac = (double)im;
      frac /= 10;
      whole = (double)(il);
      a[(im+ioff)] = frac + whole;
/*      b[(im+ioff)] = frac + whole*10;*/
      b[(im+ioff)] = frac + whole;
    }/*endfor*/
    ioff += nm;
  }/*endfor*/

/*========================================================================*/
/* Zero the matrix */

  nl_tot2 = nl_tot*nl_tot;
  for(il=1;il<=nl_tot2;il++){c[il] = 0.0;}

/*========================================================================*/
/* II.5) Create model data                                                */


  ioff=0;
  for(il=1;il<=nl_tot;il++){
   joff=0;
   for(jl=1;jl<=nl_tot;jl++){
    knd = jl + (il-1)*nl_tot;
    for(im=1;im<=nm;im++){
     ind = im + ioff;
     jnd = im + joff;
     c[knd] += a[ind]*b[jnd];
    }
    joff +=nm;
   }
   ioff +=nm;
  }


#ifdef DEBUG
  for(il=1;il<=nl_tot2;il++){
   printf("c[%d] = %.12g\n",il,c[il]);
  }
#endif


/*========================================================================*/
}/* end main */
/*========================================================================*/
