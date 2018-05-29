#define DEC_ALPHA
#define PARALLEL
#define NO_PRAGMA
#define SIMP_NINT
#define FFTW3
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <omp.h>
#include <mkl.h>
#ifdef PARALLEL
#include "mpi.h"
#else
#include "../typ_defs/mpi_f.h"
#endif
#include "../typ_defs/defines.h"
#ifdef FFTW3
#include "/home/mingchen/libPippen/include/fftw3.h"
#endif
typedef struct {
  double re;
  double im;
} zomplex;



