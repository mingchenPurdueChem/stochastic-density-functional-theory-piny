#define DEC_ALPHA
#define PARALLEL
#define NO_PRAGMA
#define SIMP_NINT
#define FFTW3
#define MKL_LAPACK

#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#ifdef FFTW3
#include "/home/mingchen/libPippen/include/fftw3.h"
#endif
#include <omp.h>
#include <mkl.h>
#ifdef PARALLEL
#include <mpi.h>
#else
#include "../typ_defs/mpi_f.h"
#endif
#include "../typ_defs/defines.h"
typedef struct {
  double re;
  double im;
} zomplex;



