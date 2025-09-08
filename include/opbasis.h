#ifndef OPBASIS_H
#define OPBASIS_H

#include<complex.h>
#include"../include/gparam.h"

typedef struct OPbasis{
      double complex **Poly_ev;
      double complex **Poly_odd;
      //double complex **Glueball_ev;
      //double complex **Glueball_odd;

      int levels;
      int n_polev=2;
      int n_polodd=1;
      //int n_gluev*levels;
      //int n_gluod*levels;


} OPbasis;

void init_op(double complex **operator,
      GParam const * const param);

void poly_line(OPbasis *Basis,
      Conf const * const GC,
      Geometry const * const geo,
      GParam const * const param);

