#ifndef GAUGE_H
#define GAUGE_H

#include<complex.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"macro.h"

typedef struct Gauge {
   double complex comp[NFLAVOUR] __attribute__((aligned(DOUBLE_ALIGN)));
} GAuge;

