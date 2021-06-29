#ifndef BLAS_EXTENDED_H
#define BLAS_EXTENDED_H

#include "blas_enum.h"
#include "blas_malloc.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef x86
#include <fpu_control.h>

/* Declaration of storage space. */
#define x86_FIX_DECL unsigned short __old_cw, __new_cw;
#define x86_FIX_DECL_MAYBE unsigned short __old_cw, __new_cw, __indigenous_flag= 0;

/* Turn off extended bit, turn on double bit. */
#define x86_FIX_START \
  _FPU_GETCW(__old_cw); \
  __new_cw = (__old_cw & ~_FPU_EXTENDED) | _FPU_DOUBLE; \
  _FPU_SETCW(__new_cw);

#define x86_FIX_START_MAYBE \
  _FPU_GETCW(__old_cw); \
  __new_cw = (__old_cw & ~_FPU_EXTENDED) | _FPU_DOUBLE; \
  if(!__indigenous_flag && (__old_cw != __new_cw)) { \
     _FPU_SETCW(__new_cw); \
  }

/* Restore original control word. */
#define x86_FIX_END  _FPU_SETCW(__old_cw);
#define x86_FIX_END_MAYBE  _FPU_SETCW(__old_cw);
/* Set maybe flag for indigenous floating operations */
#define x86_FIX_SET_INDIGENOUS __indigenous_flag = 1;
/* If x86 fix not needed, define as nothing : */
#else
#define x86_FIX_DECL
#define x86_FIX_START
#define x86_FIX_END

#define x86_FIX_DECL_MAYBE
#define x86_FIX_START_MAYBE
#define x86_FIX_END_MAYBE
#define x86_FIX_SET_INDIGENOUS
#endif


/* constants */

#define BITS_S  24
#define BITS_D  53
#define BITS_E  106

/* Split a double into 2 parts with at most 26 bits each. (2^27 + 1) */
#define split 	(134217729.0)


/* macros */

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

/* prototypes */

#include "blas_extended_proto.h"

#endif
   /* BLAS_EXTENDED_H */
