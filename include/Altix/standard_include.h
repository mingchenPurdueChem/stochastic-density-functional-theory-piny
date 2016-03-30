#define ALTIX
#define PARALLEL
#define T3E_SCILIB_OFF
#define NO_PRAGMA
#define SIMP_NINT
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#ifdef PARALLEL
#include "mpi.h"
#else
#include "../typ_defs/mpi_f.h"
#endif
#include "../typ_defs/defines.h"
typedef struct {
  double re;
  double im;
} zomplex;

