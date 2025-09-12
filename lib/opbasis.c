#ifndef OPBASIS_C
#define OPBASIS_C

#include"../include/macro.h"
#include"../include/opbasis.h"
#include<complex.h>
#include"../include/gparam.h"
#include<string.h>
#include<stdlib.h>

void decl_op(OPerator *operator,
      int n_op,
      GParam const * const param)
   {
   int err;
   err=posix_memalign((void**) operator,(size_t) DOUBLE_ALIGN, (size_t) (n_op) * sizeof(double complex *));
   if(err!=0)
      {
      fprintf(stderr, "Problems in allocating the polyakov correlators! (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   int n;
   for(n=0; n<n_op;n++)
      {
      err=posix_memalign((void**) &((*operator)[n]), (size_t) DOUBLE_ALIGN, (size_t) param->d_size[0] * sizeof(double complex));
      if(err!=0)
         {
         fprintf(stderr, "Problems in allocating the polyakov correlators! (%s, %d)\n", __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
      }
   }

void init_op(OPerator operator,
      int n_op,
      GParam const * const param)
   {
   int n,t;
   for(n=0; n<n_op;n++)
      for(t=0; t<param->d_size[0]; t++)
        {
         operator[n][t]=0;
        }
   }


void free_op(OPerator operator,
         int n_op)
   {
   int n;

   for(n=0; n<n_op; n++)
      {
      free(operator[n]);
      }
   free(operator);
   }

void decl_opbasis(OPbasis *opbasis,
      GParam const * const param)
{
   opbasis->polev=2;
   opbasis->polodd=1;
   opbasis->glueev=1;
   opbasis->glueodd=1;

   opbasis->levels=2*(param->numblock+1); //2 because I measure with and withour smearing
   opbasis->n_polev=opbasis->polev*opbasis->levels;
   opbasis->n_polodd=opbasis->polodd*opbasis->levels;
   opbasis->n_glueev=opbasis->glueev*opbasis->levels;
   opbasis->n_glueodd=opbasis->glueodd*opbasis->levels;

   decl_op(&(opbasis->Poly_odd),opbasis->n_polodd,param);
   decl_op(&opbasis->Poly_ev,opbasis->n_polev,param);
   decl_op(&opbasis->Glue_ev,opbasis->n_glueev,param);
   decl_op(&opbasis->Glue_odd,opbasis->n_glueodd,param);
}

void init_opbasis(OPbasis *opbasis,
      GParam const * const param)
   {
   init_op(opbasis->Poly_odd,opbasis->n_polodd,param);
   init_op(opbasis->Poly_ev,opbasis->n_polev,param);
   init_op(opbasis->Glue_ev,opbasis->n_glueev,param);
   init_op(opbasis->Glue_odd,opbasis->n_glueodd,param);
   }


void free_basis(OPbasis *opbasis)
   {
   free_op(opbasis->Poly_odd,opbasis->n_polodd);
   free_op(opbasis->Poly_ev,opbasis->n_polev);
   free_op(opbasis->Glue_ev,opbasis->n_glueev);
   free_op(opbasis->Glue_odd,opbasis->n_glueodd);
   }

//computes the polyakov line along direction 1 at a single point
double complex poly_line(Conf const * const GC,
      Geometry const * const geo,
      GParam const * const param,
      long r)
   {
   int i;
   double complex poly;
   poly=1.0;
   for(i=0; i<param->d_size[1];i++)
         {
         poly*=GC->lambda[r][1];
         r=nnp(geo,r,1);
         }
   return poly;
   }


//compute the polyakov line with a plaquette inserted in the plane (1,2)
//and orientation (up, down) and summ over it's position to project to zero momentum
//returns the even and odd operator
void poly_plaq(Conf const * const GC,
      Geometry const * const geo,
      GParam const * const param,
      long r,
      double complex polyline,
      double complex *polyev,
      double complex *polyodd)
   {
   int i;
   double complex polypup, polypdown;
   polypup=polypdown=0;

   //then insert a plaquette in each site and project to zero momentum
   for(i=0; i<param->d_size[1];i++)
         {
         polypup+=polyline*conj(plaquette_complex(GC, geo, r, 2, 1));
         polypdown+=polyline*plaquette_complex(GC, geo, nnm(geo,r,2), 2, 1);
         r=nnp(geo,r,1);
         }
   *polyev+=(polypup+polypdown)/(2*param->d_size[1]);
   *polyodd+=(polypup-polypdown)/(2*param->d_size[1]);
   }


void poly_averaged(OPbasis *basis,
      Conf const * const GC,
      Geometry const * const geo,
      GParam const * const param,
      int level)
   {

   double complex tmp_polyline,polyline,polyev,polyodd;


   int i,t;
   long r,r2;
   r=0;
   r2=0;
   for(t=0; t<param->d_size[0]; t++)
      {
      tmp_polyline=polyline=polyev=polyodd=0;
      for(i=0; i<param->d_size[2]; i++)
         {
         tmp_polyline=poly_line(GC,geo,param,r2);
         poly_plaq(GC,geo,param,r2,tmp_polyline,&polyev,&polyodd);
         polyline+=tmp_polyline;
         r2=nnp(geo,r2,2);
         }
      basis->Poly_ev[basis->polev*level][t]=polyline/param->d_size[2];
      basis->Poly_ev[basis->polev*level+1][t]=polyev/param->d_size[2];
      basis->Poly_odd[basis->polodd*level][t]=polyodd/param->d_size[2];
      r=nnp(geo,r,0);
      r2=r;
      }

   }

void glueball_averaged(OPbasis *basis,
         Conf const * const GC,
         Geometry const * const geo,
         GParam const * const param,
         int level)
   {
   double complex tmp_plaq;
   long r,r2,r3;
   int i,j,t;
   r=0;
   r2=0;
   r3=0;

   for(t=0; t<param->d_size[0]; t++)
      {
      tmp_plaq=0;
      r2=r;
      for(i=0; i<param->d_size[1]; i++)
         {
         r3=r2;
         for(j=0; j<param->d_size[2]; j++)
         {
            tmp_plaq+=plaquette_complex(GC, geo, r3, 1, 2);
            r3=nnp(geo,r3,2);
         }
         r2=nnp(geo,r2,1);
         }
      basis->Glue_ev[basis->glueev*level][t]=(tmp_plaq+conj(tmp_plaq))/(param->d_size[1]*param->d_size[2]);
      basis->Glue_odd[basis->glueodd*level][t]=(tmp_plaq-conj(tmp_plaq))/(param->d_size[1]*param->d_size[2]);
      r=nnp(geo,r,0);
      }
   }

void measure_print_corr(double complex ** operators,
      GParam const * const param,
      FILE * datafilep,
      int nops)
   {
   int t1, t2, t ;
   int i,j;
   double complex corr;
   for(t = 0; t<param->d_size[0]/2; t++)
      for(i=0; i<nops; i++)
         for(j=0; j<=i; j++)
        {
        corr = 0.0;
        for(t1 = 0; t1<param->d_size[0]; t1++)
           {
           t2=(t1+t) % param->d_size[0];
           corr += operators[i][t2]*conj(operators[j][t1]);
           }
        corr/=(double) param->d_size[0];
        fprintf(datafilep, "%.12f %.12f ", creal(corr), cimag(corr));
        }
   }

void measure_print_corr_all(OPbasis *opbasis,
                           GParam const * const param,
                           FILE * datafilep)
   {
   measure_print_corr(opbasis->Poly_ev,param,datafilep,opbasis->n_polev);
   measure_print_corr(opbasis->Poly_odd,param,datafilep,opbasis->n_polodd);
   measure_print_corr(opbasis->Glue_ev,param,datafilep,opbasis->n_glueev);
   measure_print_corr(opbasis->Glue_odd,param,datafilep,opbasis->n_glueodd);
   fprintf(datafilep, "\n");
   fflush(datafilep);
   }

#endif










