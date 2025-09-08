#ifndef OPBASIS_C
#define OPBASIS_C

#include"../include/macro.h"
#include"../include/opbasis.h"
#include<complex.h>
#include"../include/gparam.h"
#include<string.h>
#include<stdlib.h>


void decl_opbasis(OPbasis *opbasis,
      GParam const * const param)
{
   opbasis->levels=(param->smearing_steps+1)*(param->numblock+1);
   opbasis->n_polev=2*opbasis->levels;
   opbasis->n_polodd=opbasis->levels;
   int err;

   err=posix_memalign((void**) &(opbasis->Poly_ev),(size_t) DOUBLE_ALIGN, (size_t) (opbasis->n_polev) * sizeof(double complex *));
   if(err!=0)
      {
      fprintf(stderr, "Problems in allocating the polyakov correlators! (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   int n;
   for(n=0; n<opbasis->n_polev;n++)
      {
      err=posix_memalign((void**) &(opbasis->Poly_ev[n]), (size_t) DOUBLE_ALIGN, (size_t) param->d_size[0] * sizeof(double complex));
      if(err!=0)
         {
         fprintf(stderr, "Problems in allocating the polyakov correlators! (%s, %d)\n", __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
      }


   err=posix_memalign((void**) &(opbasis->Poly_odd),(size_t) DOUBLE_ALIGN, (size_t) (opbasis->n_polodd) * sizeof(double complex *));
   if(err!=0)
      {
      fprintf(stderr, "Problems in allocating the polyakov correlators! (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   for(n=0; n<opbasis->n_polodd;n++)
      {
      err=posix_memalign((void**) &(opbasis->Poly_odd[n]), (size_t) DOUBLE_ALIGN, (size_t) param->d_size[0] * sizeof(double complex));
      if(err!=0)
         {
         fprintf(stderr, "Problems in allocating the polyakov correlators! (%s, %d)\n", __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
      }

}

void init_opbasis(OPbasis *opbasis,
      GParam const * const param)
   {

   int n,t;
   for(n=0; n<opbasis->n_polev;n++)
      for(t=0; t<param->d_size[0]; t++)
        {
         opbasis->Poly_ev[n][t]=0;
        }

   for(n=0; n<opbasis->n_polodd;n++)
      for(t=0; t<param->d_size[0]; t++)
        {
         opbasis->Poly_odd[n][t]=0;
        }
   }

void free_basis(OPbasis *opbasis)
   {
   int n;

   for(n=0; n<opbasis->n_polev; n++)
      {
      free(opbasis->Poly_ev[n]);
      }
   free(opbasis->Poly_ev);

   for(n=0; n<opbasis->n_polodd; n++)
      {
      free(opbasis->Poly_odd[n]);
      }
   free(opbasis->Poly_odd);
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
         polypup+=polyline*plaquette_single(GC, geo, r, 1, 2);
         polypdown+=polyline*plaquette_single(GC, geo, nnm(geo,r,2), 2, 1);
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
   tmp_polyline=polyline=polyev=polyodd=0;

   int i,t;
   long r,r2;
   r=0;
   r2=0;
   for(t=0; t<param->d_size[0]; t++)
      {
   for(i=0; i<param->d_size[2]; i++)
      {
      tmp_polyline=poly_line(GC,geo,param,r2);
      poly_plaq(GC,geo,param,r2,tmp_polyline,&polyev,&polyodd);
      polyline+=tmp_polyline;
      r2=nnp(geo,r2,2);
      }
   basis->Poly_ev[2*level][t]=polyline; //2 because two operators, need to be included as a parameter
   basis->Poly_ev[2*level+1][t]=polyev;
   basis->Poly_odd[level][t]=polyodd;
   r=nnp(geo,r,0);
   r2=r;
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
           corr += operators[i][t2]*operators[j][t1];
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
   fprintf(datafilep, "\n");
   fflush(datafilep);
   }

#endif










