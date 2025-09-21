#ifndef OPBASIS_H
#define OPBASIS_H

#include"../include/macro.h"
#include<complex.h>
#include <stdlib.h>

#include"../include/gparam.h"
#include"../include/conf.h"
#include"../include/geometry.h"

typedef double complex** OPerator;

typedef struct OPbasis{

      OPerator Poly_ev;
      OPerator Poly_odd;
      OPerator Glue_ev;
      OPerator Glue_odd;

      int polev;
      int polodd;
      int glueev;
      int glueodd;

      int levels;

      int n_polev;
      int n_polodd;
      int n_glueev;
      int n_glueodd;


} OPbasis;

void decl_op(OPerator *op,
      int n_op,
      GParam const * const param);

void init_op(OPerator op,
      int n_op,
      GParam const * const param);

void free_op(OPerator op,
      int n_op);

void decl_opbasis(OPbasis *opbasis,
      GParam const * const param);

void init_opbasis(OPbasis *opbasis,
      GParam const * const param);

void free_basis(OPbasis *opbasis);

double complex poly_line(Conf const * const GC,
      Geometry const * const geo,
      GParam const * const param,
      long r);

void poly_plaq(Conf const * const GC,
      Geometry const * const geo,
      GParam const * const param,
      long r,
      double complex polyline,
      double complex *polyev,
      double complex *polyodd);

void poly_two_plaq(Conf const * const GC,
      Geometry const * const geo,
      GParam const * const param,
      long r,
      double complex polyline,
      double complex *polyev,
      double complex *polyodd);

void poly_averaged(OPbasis *basis,
      Conf const * const GC,
      Geometry const * const geo,
      GParam const * const param,
      int level);

void glueball_averaged(OPbasis *basis,
         Conf const * const GC,
         Geometry const * const geo,
         GParam const * const param,
         int level);

void measure_print_corr(double complex ** operators,
      GParam const * const param,
      FILE * datafilep,
      int nops);

void measure_print_corr_all(OPbasis *opbasis,
                           GParam const * const param,
                           FILE * datafilepev,
                           FILE * datafilepodd,
                           FILE * datafilegev,
                           FILE * datafilegodd);
#endif
