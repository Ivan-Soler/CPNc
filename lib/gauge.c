#ifndef VEC_C
#define VEC_C

#include<complex.h>
#include<math.h>
#include<stdio.h>
#include<string.h>

#include"../include/endianness.h"
#include"../include/macro.h"
#include"../include/random.h"
#include"../include/vec.h"


// A=1
void one_Vec(Gauge * restrict A);


// A=0
void zero_Vec(Vec * restrict A);


// A=B
void equal_Vec(Vec * restrict A, Vec const * const restrict B);


// A=B
void equal_cc_Vec(Vec * restrict A, Vec const * const restrict B);


// A+=B
void plus_equal_Vec(Vec * restrict A, Vec const * const restrict B);


// A-=B
void minus_equal_Vec(Vec * restrict A, Vec const * const restrict B);


// A*=r
void times_equal_real_Vec(Vec * restrict A, double r);


// A*=c
void times_equal_complex_Vec(Vec * restrict A, double complex c);

// A*=c*
void times_equal_complex_conj_Vec(Vec * restrict A, double complex * c);
