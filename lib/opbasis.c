#ifndef OPBASIS_C
#define OPBASIS_C

#include"../include/opbasis.h"
#include<complex.h>
#include"../include/gparam.h"



void init_op(OPbasis *opbasis,
      GParam const * const param,
      int n_ops)
{


   opbasis->levels=param->smearing_steps*param->numblock;
   opbasis->n_polev=2*opbasis->levels;
   opbasis->n_polodd=opbasis->levels;
   int err;
   int t;
   err=posix_memalign((void**) &(opbasis->Poly_ev),(size_t) DOUBLE_ALIGN, (size_t) (opbasis->n_polev) * sizeof(double complex *));
   if(err!=0)
      {
      fprintf(stderr, "Problems in allocating the polyakov correlators! (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   int n;
   for(n=0; n<n_polev;n++)
      {
      err=posix_memalign((void**) &(opbasis->Poly_ev[n]), (size_t) DOUBLE_ALIGN, (size_t) param->d_size[0] * sizeof(double complex));
      if(err!=0)
         {
         fprintf(stderr, "Problems in allocating the polyakov correlators! (%s, %d)\n", __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
      }

   for(n=0; n<opbasis->n_polev;n++)
      for(t=0; t<param->d_size[0]; t++)
        {
         opbasis->Poly_ev[t][n]=0;
        }

   int err;
   int t;
   err=posix_memalign((void**) &(opbasis->Poly_odd),(size_t) DOUBLE_ALIGN, (size_t) (opbasis->n_polodd) * sizeof(double complex *));
   if(err!=0)
      {
      fprintf(stderr, "Problems in allocating the polyakov correlators! (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
   int n;
   for(n=0; n<opbasis->n_polodd;n++)
      {
      err=posix_memalign((void**) &(opbasis->Poly_odd[n]), (size_t) DOUBLE_ALIGN, (size_t) param->d_size[0] * sizeof(double complex));
      if(err!=0)
         {
         fprintf(stderr, "Problems in allocating the polyakov correlators! (%s, %d)\n", __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
      }

   for(n=0; n<opbasis->n_polodd;n++)
      for(t=0; t<param->d_size[0]; t++)
        {
         opbasis->Poly_odd[t][n]=0;
        }

}

//computes the polyakov line along direction 1 at a single point
void poly_line(Conf const * const GC,
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
   polypup=polypdown=1.0

   //then insert a plaquette in each site and project to zero momentum
   for(i=0; i<param->d_size[1];i++)
         {
         polypup+=polyline*plaquette_single(GC, geo, r2, 1, 2)
         polypdown+=polyline*plaquette_single(GC, geo, r2, 1, 2)
         r=nnp(geo,r,1);
         }
   *polyev+=(polypup+polypdown)/(2*param->d_size[1])
   *polyodd+=(polypup-polydown)/(2*param->d_size[1])
   }


void poly_averaged(OPbasis basis,
      Conf const * const GC,
      Geometry const * const geo,
      Gparam const * const param,
      int level)
   {

   double complex polyev,polyodd;
   polyev=polyodd=0;

   int i,r,t,r2;
   r=0;
   r2=0
   for(t=0; t<param->d_size[0]; t++)
      {
   for(i=0; i<param->d_size[2]; i++)
      {
      polyline+=poly_line(GC,geo,param,r2);
      poly_plaq(GC,geo,param,r2,polyline,&polyev,&polyodd);
      r2=nnp(geo,r2,2);
      }
   basis->Poly_ev[t][level+0]=polyline;
   basis->Poly_ev[t][level+1]=polyev;
   basis->Poly_odd[t][level]=polyodd;
   r=nnp(geo,r,0);
   r2=r;

   }











