#define DEC_ALPHA
#define PARALLEL
#define NO_PRAGMA
#define SIMP_NINT
#define FFTW3
#define MKL_LAPACK
#define _GNU_SOURCE
//#define FAST_FILTER

#include <fenv.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <omp.h>
#ifdef FFTW3
#include "fftw3.h"
#endif
#include <mkl.h>
#ifdef PARALLEL
#include <mpi.h>
#else
#include "../typ_defs/mpi_f.h"
#endif
#include "../typ_defs/defines.h"
#include <fcntl.h>
#include <unistd.h>

typedef struct {
  double re;
  double im;
} zomplex;



