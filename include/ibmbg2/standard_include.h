#define IBM_ESSL
#define PARALLEL
#define NO_PRAGMA
#define SIMP_NINT
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#ifdef PARALLEL
#include "/bgl/BlueLight/ppcfloor/bglsys/include/mpi.h"
#else
#include "../typ_defs/mpi_f.h"
#endif
#include "../typ_defs/defines.h"
typedef struct {
  double re;
  double im;
} zomplex;
