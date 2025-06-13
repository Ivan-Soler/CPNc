#ifndef CONF_MEAS_C
#define CONF_MEAS_C

#include"../include/macro.h"

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<complex.h>

#include"../include/flavour_matrix.h"
#include"../include/gparam.h"
#include"../include/geometry.h"
#include"../include/conf.h"

// computation of the plaquette in position r and positive directions i,j
double plaquette_single(Conf const * const GC,
                        Geometry const * const geo,
                        long r,
                        int i,
                        int j)
   {

//
//       ^ i
//       |  (3)
//       +---<---+
//       |       |
//   (4) V       ^ (2)
//       |       |
//       +--->---+---> j
//       r  (1)
//

   double complex ris;

   #ifdef CSTAR_BC
     ris = GC->lambda[r][j];  // (1)
     if(bcsitep(geo, r, j)==1){ ris *= GC->lambda[nnp(geo, r, j)][i]; } // (2)
     else { ris *= conj(GC->lambda[nnp(geo, r, j)][i]); }
     if(bcsitep(geo, r, i)==1){ ris *= conj(GC->lambda[nnp(geo, r, i)][j]); } // (3)
     else{ ris *= GC->lambda[nnp(geo, r, i)][j]; }
     ris *= conj(GC->lambda[r][i]);
   #else
     ris = GC->lambda[r][j];  // (1)
     ris *= GC->lambda[nnp(geo, r, j)][i]; // (2)
     ris *= conj(GC->lambda[nnp(geo, r, i)][j]); // (3)
     ris *= conj(GC->lambda[r][i]);
   #endif

   return creal(ris);
   }


double plaquette(Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param)
   {
   long r;
   double ris=0.0;

   for(r=0; r<(param->d_volume); r++)
      {
      double tmp;
      int i, j;

      i=0;
      tmp=0.0;
     
      for(i=0; i<STDIM; i++)
         {
         for(j=i+1; j<STDIM; j++)
            {
            tmp+=plaquette_single(GC, geo, r, i, j);
            }
         }

      ris+=tmp;
      }

   ris*=param->d_inv_vol;
   ris/=((double) STDIM*((double) STDIM-1.0)/2.0);

   return ris;
   }


// compute the average value of Re[ phi_x^{dag} lambda_{x,mu} phi_{x+mu} ]
double higgs_interaction(Conf const * const GC,
                         Geometry const * const geo,
                         GParam const * const param)
  {
  int i;
  long r;
  double aux, ris=0.0;
  Vec v1;

  for(r=0; r<(param->d_volume); r++)
     {
     aux=0.0;

     for(i=0; i<STDIM; i++)
        {
        #ifdef CSTAR_BC
          if(bcsitep(geo, r, i)==1)
            {
            equal_Vec(&v1, &(GC->phi[nnp(geo, r, i)]));
            }
          else
            {
            equal_cc_Vec(&v1, &(GC->phi[nnp(geo, r, i)]));
            }
        #else
          equal_Vec(&v1, &(GC->phi[nnp(geo, r, i)]));
        #endif
        times_equal_complex_Vec(&v1, chargepow(GC->lambda[r][i]) );

        aux+= creal(scal_prod_Vec(&(GC->phi[r]), &v1) );
        }

     ris+=aux;
     }

  ris/=(double) STDIM;
  ris*=param->d_inv_vol;

  return ris;
  }

double gauge_fixing_interaction(Conf const * const GC,
		Geometry const * const geo,
		GParam const * const param)
	{
	int i;
	long r;
	double ris=0.0;

	for(i=0; i<STDIM; i++)
		{
		for(r=0; r<param->d_volume; r++)
			{
				ris+=creal(conj(GC->gauge[r])*GC->lambda[r][i]*GC->gauge[nnp(geo,r,i)]);
			}

		}
	ris*=-param->d_quench_gamma;
	return ris;
	}


// return the average of Re(link)
double realpartlink(Conf const * const GC,
                    GParam const * const param)
  {
  int i;
  long r;
  double ris=0.0;

  for(i=0; i<STDIM; i++)
     {
     for(r=0; r<param->d_volume; r++)
        {
        ris+=creal(GC->lambda[r][i]);
        }
     }

  ris*=param->d_inv_vol;
  ris/=((double) STDIM);

  return ris;
  }

double imagpartlink(Conf const * const GC,
                    GParam const * const param)
  {
  int i;
  long r;
  double ris=0.0;

  for(i=0; i<STDIM; i++)
     {
     for(r=0; r<param->d_volume; r++)
        {
        ris+=cimag(GC->lambda[r][i]);
        }
     }

  ris*=param->d_inv_vol;
  ris/=((double) STDIM);

  return ris;
  }


//Flux and monopole related stuff//
//Flux through a plaquette//
double flux_plaquette(Conf const * const GC,
        Geometry const * const geo,
        long r,
        int j,
        int i)
   {
   double ris;
   ris=0;
   ris+=carg(GC->lambda[r][j]);
   ris+=carg(GC->lambda[nnp(geo,r,j)][i]);
   ris-=carg(GC->lambda[nnp(geo,r,i)][j]);
   ris-=carg((GC->lambda[r][i]));

   return(fmod(ris,2*PI));
   }

//Monopoles inside a cube//

/*
 *   +----------+
    /         / |
   +----<----+  |
   |         |  +
   V         ^ /
   |         |/
   +---->----+




   ^  (2)
   |
   |     (1)
   |  /
   | /
   |/
   +-------------> (0)


               */

long monpoles_cube(Conf const * const GC,
        Geometry const * const geo,
        long r,
        GParam const * const param)
   {
   double flux;
   int ris;

   flux=0;
   flux-=flux_plaquette(GC,geo,r,0,2);  //front (flow points outside)
   flux+=flux_plaquette(GC,geo,r,0,1);  //bottom (flow points inside)
   flux+=flux_plaquette(GC,geo,r,1,2);  //left  (flow points inside)

   flux+=flux_plaquette(GC,geo,nnp(geo,r,1),0,2);  //back (flow points inside)
   flux-=flux_plaquette(GC,geo,nnp(geo,r,2),0,1);  //top (flow points outside)
   flux-=flux_plaquette(GC,geo,nnp(geo,r,0),1,2);  //right (flow points outside)

#ifdef DEBUG
   if(param->d_start==0){
      flux=0;
      flux+=flux_plaquette(GC,geo,r,0,1);  //front (flow points outside)
      flux-=flux_plaquette(GC,geo,nnp(geo,r,2),0,1);  //back (flow points inside)
   }

#endif
   ris=(int)(flux/(2*PI));

#ifdef DEBUG
   if(abs(ris)>0)
   {
   fprintf(stderr, " %d Monopoles at r=%ld \n",ris,r);
   }
#endif

   return(ris);
   }

long measure_monopoles(Conf const * const GC,
        Geometry const * const geo,
        GParam const * const param)
   {
   long r,monopoles;
   monopoles=0;

   for(r=0; r<param->d_volume; r++)
      {
      monopoles+=monpoles_cube(GC,geo,r,param);
      }
   return(monopoles);
   }


// compute flavour related observables in the tensor channel
//
// GC->Qh needs to be initialized before calling this function
//
// tildeG0=Tr[(\sum_x Q_x)(\sum_y Q_y)]/volume
// tildeGminp=ReTr[(\sum_x Q_xe^{ipx})(\sum_y Q_ye^{-ipy)]/volume
//
// tildeG0 is the susceptibility, tildeGminp is used to compute the 2nd momentum correlation function
//
void compute_flavour_observables_tensor(Conf const * const GC,
                                        GParam const * const param,
                                        double *tildeG0,
                                        double *tildeGminp)
  {
  int coord[STDIM];
  long r;
  const double p = 2.0*PI/(double)param->d_size[1];
  FMatrix Q, Qp, Qmp, tmp1, tmp2;

  // Q =sum_x Q_x
  // Qp=sum_x e^{ipx}Q_x
  // Qmp=sum_x e^{-ipx}Q_x

  zero_FMatrix(&Q);
  zero_FMatrix(&Qp);
  zero_FMatrix(&Qmp);
  for(r=0; r<(param->d_volume); r++)
     {
     equal_FMatrix(&tmp1, &(GC->Qh[r]));
     equal_FMatrix(&tmp2, &tmp1);

     plus_equal_FMatrix(&Q, &tmp1);

     si_to_cart(coord, r, param);

     times_equal_complex_FMatrix(&tmp1, cexp(I*((double)coord[1])*p));
     plus_equal_FMatrix(&Qp, &tmp1);

     times_equal_complex_FMatrix(&tmp2, cexp(-I*((double)coord[1])*p));
     plus_equal_FMatrix(&Qmp, &tmp2);
     }

  equal_FMatrix(&tmp1, &Q);
  times_equal_FMatrix(&tmp1, &Q);

  *tildeG0=retr_FMatrix(&tmp1)*param->d_inv_vol;

  equal_FMatrix(&tmp1, &Qp);
  times_equal_FMatrix(&tmp1, &Qmp);
  *tildeGminp=retr_FMatrix(&tmp1)*param->d_inv_vol;
  }

// perform tensor-related measures and save results in a buffer
void perform_tensor_measures_buffer(Conf *GC,
                                 GParam const * const param,
                                 Geometry const * const geo,
                                 double buffer[5])
   {
   (void) geo; // kept just for consistency with other measures
   double tildeG0_v, tildeGminp_v;

   compute_flavour_observables_tensor(GC,
                                      param,
                                      &tildeG0_v,
                                      &tildeGminp_v);

   buffer[0]=tildeG0_v;
   buffer[1]=tildeGminp_v;
   }

// compute flavour related observables in the vector channel
//
// tildeG0=(\sum_x z_x)^{dag}(\sum_y z_y)/volume
// tildeGminp=(\sum_x z_xe^{ipx})^{dag}(\sum_y z_ye^{ipy)]/volume
//
// tildeG0 is the susceptibility, tildeGminp is used to compute the 2nd momentum correlation function
//
void compute_flavour_observables_vector(Conf const * const GC,
                                        GParam const * const param,
                                        double *tildeG0,
                                        double *tildeGminp)
  {
  int coord[STDIM];
  long r;
  const double p = 2.0*PI/(double)param->d_size[1];
  Vec V, Vp, tmp1;

  // V =sum_x Q_x
  // Vp=sum_x e^{ipx}Q_x

  zero_Vec(&V);
  zero_Vec(&Vp);
  for(r=0; r<(param->d_volume); r++)
     {
     equal_Vec(&tmp1, &(GC->phi[r]));

     plus_equal_Vec(&V, &tmp1);

     si_to_cart(coord, r, param);

     times_equal_complex_Vec(&tmp1, cexp(I*((double)coord[1])*p));
     plus_equal_Vec(&Vp, &tmp1);

     }

  equal_Vec(&tmp1, &V);
  *tildeG0=creal(scal_prod_Vec(&tmp1, &V))*param->d_inv_vol;

  equal_Vec(&tmp1, &Vp);
  *tildeGminp=creal(scal_prod_Vec(&tmp1, &Vp))*param->d_inv_vol;

  }

// perform vector-related measures and save results in a buffer
void perform_vec_measures_buffer(Conf *GC,
                                 GParam const * const param,
                                 Geometry const * const geo,
                                 double buffer[5])
   {
   (void) geo; // kept just for consistency with other measures
   double tildeG0_v, tildeGminp_v;

   compute_flavour_observables_vector(GC,
                                      param,
                                      &tildeG0_v,
                                      &tildeGminp_v);

   buffer[0]=tildeG0_v;
   buffer[1]=tildeGminp_v;
   }

// performe gauge measures and save in a buffer
void perform_gauge_measures_buffer(Conf const * const GC,
                                     GParam const * const param,
                                     Geometry const * const geo,
                                     double buffer[5])
  {
	(void) geo;
	int coord[STDIM];
	long r;
	const double p = 2.0*PI/(double)param->d_size[1];
	double complex C, Cp;

	C=0;
	Cp=0;

	// C = sum_x gauge_{x,\mu}
	// Cp= sum_x e^{ipx} gauge_{x,\mu}
	for (r=0; r<(param->d_volume); r++){
		C+=GC->lambda[r][2];

		si_to_cart(coord,r,param);

		Cp+=cexp(I*((double)coord[1])*p)*GC->lambda[r][2];

	}

	C=C*conj(C);
	Cp=Cp*conj(Cp);

	buffer[2]=creal(C)*param->d_inv_vol;
	buffer[3]=creal(Cp)*param->d_inv_vol;

  }

// performe overlap measures and save in a buffer
void perform_overlap_measures_buffer(Conf const * const GC,
                                     Conf const * const GC2,
                                     GParam const * const param,
                                     Geometry const * const geo,
                                     double buffer[5])
  {
  long int r;
  (void) geo; // geo is used only for consistency with other cases
                // so it is possibile to add \xi computation without changes
  double complex overlap;
  overlap=0.0;

  for(r=0; r<param->d_volume; r++)
     {
     overlap+=conj(GC->gauge[r])*GC2->gauge[r];
     }

  buffer[4]=(double) (conj(overlap)*overlap*param->d_inv_vol);

  }

void perform_measures(Conf *GC,
                      GParam const * const param,
                      Geometry const * const geo,
                      FILE *datafilep)
   {
   long r;

   double tildeG0_t, tildeGminp_t;
   double tildeG0_v, tildeGminp_v;
   double scalar_coupling, plaq, relink;
   long monopoles;


   for(r=0; r<(param->d_volume); r++)
      {
      init_FMatrix(&(GC->Qh[r]), &(GC->phi[r]));
      }

   compute_flavour_observables_tensor(GC,
                                      param,
                                      &tildeG0_t,
                                      &tildeGminp_t);

   compute_flavour_observables_vector(GC,
                                      param,
                                      &tildeG0_v,
                                      &tildeGminp_v);
	#ifdef DEBUG_GAUGE_INV
   double buffer[5];
   compute_flavour_observables_vector(GC,
                                      param,
                                      &tildeG0_v,
                                      &tildeGminp_v);

   perform_gauge_measures_buffer(GC,
		   	   	   	   	   	   	 param,
								 geo,
								 buffer);
   tildeG0_t=buffer[2];
   tildeGminp_t=buffer[3];
	#endif
   scalar_coupling=higgs_interaction(GC, geo, param);
   plaq=plaquette(GC, geo, param);
   relink=realpartlink(GC, param);

   monopoles=measure_monopoles(GC,geo,param);

   fprintf(datafilep, "%.12g %.12g ", tildeG0_t, tildeGminp_t);
   fprintf(datafilep, "%.12g %.12g ", tildeG0_v, tildeGminp_v);
   fprintf(datafilep, "%.12g %.12g ", scalar_coupling, plaq);
   fprintf(datafilep, "%.12g %ld", relink,monopoles);
   fprintf(datafilep, "\n");

   fflush(datafilep);
   }


#endif
