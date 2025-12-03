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
     ris *= conj(GC->lambda[r][i]); //4
   #endif

   return creal(ris);
   }

double complex plaquette_complex(Conf const * const GC,
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
     ris *= conj(GC->lambda[r][i]); //4
   #endif

   return ris;
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

double plaquette_spatial(Conf const * const GC,
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

      for(i=1; i<STDIM; i++)
         {
         for(j=i+1; j<STDIM; j++)
            {
            tmp+=plaquette_single(GC, geo, r, i, j);
            }
         }

      ris+=tmp;
      }

   ris*=param->d_inv_vol;
   ris/=((double) (STDIM-1)*((double) STDIM-2.0)/2.0);

   return ris;
   }

double plaquette_temporal(Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param)
   {
   long r;
   double ris=0.0;

   for(r=0; r<(param->d_volume); r++)
      {
      double tmp;
      int j;
      tmp=0.0;

      for(j=1; j<STDIM; j++)
         {
         tmp+=plaquette_single(GC, geo, r, 0, j);
         }
      ris+=tmp;
      }

   ris*=param->d_inv_vol;
   ris/=((double) (STDIM-1));

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
void remove_strings(double * flux)
   {

   double precision=1.0e-12;
   while((fabs(*flux)-PI)>precision)
      {
      if(*flux>0)
         {
         *flux-=2*PI;
         }
      else
         {
         *flux+=2*PI;
         }
      }

   }
void print_links_cube(Conf const * const GC,
        Geometry const * const geo,
        long r)
   {
   int i;
   long rr;
   double angle;
   for(i=0; i<STDIM; i++ )
      {
      angle=carg(GC->lambda[r][i]);
      fprintf(stdout,"Link at r=%ld direction %d angle=%.4g \n", r, i,angle);
      }
   for(i=0; i<STDIM;i++)
      {
   rr=nnp(geo,r,i);
   angle=carg(GC->lambda[rr][(i+1)%STDIM]);
   fprintf(stdout,"Link at r=%ld+%d direction %d angle=%.4g \n", r, i, (i+1)%STDIM, angle);
   angle=carg(GC->lambda[rr][(i+2)%STDIM]);
   fprintf(stdout,"Link at r=%ld+%d direction %d angle=%.4g \n", r, i, (i+2)%STDIM,angle);
      }

   rr=nnp(geo,r,0);
   rr=nnp(geo,rr,1);
   rr=nnp(geo,rr,2);
   for(i=0; i<STDIM; i++)
      {
      angle=carg(GC->lambda[nnm(geo,rr,i)][i]);
      fprintf(stdout,"Link at r=%ld+0+1+2 direction %d angle=%.4g \n", r, i,angle);
      }

   }
double flux_plaquette(Conf const * const GC,
        Geometry const * const geo,
        long r,
        int j,
        int i)
   {
   double flux;
   flux=carg(GC->lambda[r][j]);
   flux+=carg(GC->lambda[nnp(geo,r,j)][i]);
   flux-=carg(GC->lambda[nnp(geo,r,i)][j]);
   flux-=carg((GC->lambda[r][i]));

   remove_strings(&flux);

   return(flux);
   }

double flux_plaquette_new(Conf const * const GC,
        Geometry const * const geo,
        long r,
        int j,
        int i,
        double angle_new)
   {
   double flux;
   flux=angle_new;
   flux+=carg(GC->lambda[nnp(geo,r,j)][i]);
   flux-=carg(GC->lambda[nnp(geo,r,i)][j]);
   flux-=carg((GC->lambda[r][i]));

   remove_strings(&flux);

   return(flux);
   }

double flux_inverted_plaquette_new(Conf const * const GC,
        Geometry const * const geo,
        long r,
        int j,
        int i,
        double angle_new)
   {
   double flux;
   flux=carg(GC->lambda[r][j]);
   flux+=carg(GC->lambda[nnp(geo,r,j)][i]);
   flux-=angle_new;
   flux-=carg((GC->lambda[r][i]));

   remove_strings(&flux);

   return(flux);
   }


int local_action_monopoles(Conf const *GC,
      Geometry const * geo,
      long r,
      int i,
      double angle_new,
      int Q[4])
   {

   int deltaQ;
   double old_flux[4];
   double new_flux[4];
   double delta_flux[4];
   int j, k;
   long rr;

   j=(i+1)%STDIM;
   k=(j+1)%STDIM;

   old_flux[0]=flux_plaquette(GC,geo,r,i,j);
   old_flux[1]=flux_plaquette(GC,geo,nnm(geo,r,k),i,k);
   old_flux[2]=flux_plaquette(GC,geo,nnm(geo,r,j),i,j);
   old_flux[3]=flux_plaquette(GC,geo,r,i,k);

   new_flux[0]=flux_plaquette_new(GC,geo,r,i,j,angle_new);
   new_flux[1]=flux_inverted_plaquette_new(GC,geo,nnm(geo,r,k),i,k,angle_new);
   new_flux[2]=flux_inverted_plaquette_new(GC,geo,nnm(geo,r,j),i,j,angle_new);
   new_flux[3]=flux_plaquette_new(GC,geo,r,i,k,angle_new);

   delta_flux[0]=-old_flux[0]+new_flux[0];
   delta_flux[1]=-old_flux[1]+new_flux[1];
   delta_flux[2]=-old_flux[2]+new_flux[2];
   delta_flux[3]=-old_flux[3]+new_flux[3];

   Q[0]=(int) round((GC->charge[r]*2*PI+delta_flux[0]-delta_flux[3])/(2*PI));
   deltaQ=abs(Q[0])-abs(GC->charge[r]);

   rr=nnm(geo,r,k);
   Q[1]= (int) round((GC->charge[rr]*2*PI-delta_flux[0]-delta_flux[1])/(2*PI));
   deltaQ+=abs(Q[1])-abs(GC->charge[rr]);

   rr=nnm(geo,nnm(geo,r,k),j);
   Q[2]= (int) round((GC->charge[rr]*2*PI-delta_flux[2]+delta_flux[1])/(2*PI));
   deltaQ+=abs(Q[2])-abs(GC->charge[rr]);

   rr=nnm(geo,r,j);
   Q[3]= (int) round((GC->charge[rr]*2*PI+delta_flux[2]+delta_flux[3])/(2*PI));
   deltaQ+=abs(Q[3])-abs(GC->charge[rr]);

   return (deltaQ);

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

int monpoles_cube(Conf const * const GC,
        Geometry const * const geo,
        long r,
        GParam const * const param)
   {
   double flux;
   int ris;

   (void) param;
   flux=-flux_plaquette(GC,geo,r,0,2);  //front (flow points outside)
   flux+=flux_plaquette(GC,geo,r,0,1);  //bottom (flow points inside)
   flux+=flux_plaquette(GC,geo,r,1,2);  //left  (flow points inside)

   flux+=flux_plaquette(GC,geo,nnp(geo,r,1),0,2);  //back (flow points inside)
   flux-=flux_plaquette(GC,geo,nnp(geo,r,2),0,1);  //top (flow points outside)
   flux-=flux_plaquette(GC,geo,nnp(geo,r,0),1,2);  //right (flow points outside)

   ris=(int)round(flux/(2*PI));

#ifdef DEBUG
   if(abs(ris)>0)
   {
   //fprintf(stderr, " %d Monopoles at r=%ld flux= %.8g \n",ris,r,flux);
   //print_links_cube(GC,geo,r);
   //exit(EXIT_FAILURE);
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
      monopoles+=abs(monpoles_cube(GC,geo,r,param));
      }
   return(monopoles);
   }

long action_monopoles(Conf const * const GC,
        Geometry const * const geo,
        GParam const * const param,
        long r,
        int i)
   {
   long rr,monopoles;

   monopoles=0;
   monopoles+=abs(monpoles_cube(GC,geo,r,param));
   rr=nnm(geo,r,(i+1)%STDIM);
   monopoles+=abs(monpoles_cube(GC,geo,rr,param));
   rr=nnm(geo,r,(i+2)%STDIM);
   monopoles+=abs(monpoles_cube(GC,geo,rr,param));
   rr=nnm(geo,r,(i+1)%STDIM);
   rr=nnm(geo,rr,(i+2)%STDIM);
   monopoles+=abs(monpoles_cube(GC,geo,rr,param));

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
   fprintf(datafilep, "%.12g %ld", relink, monopoles);
   fprintf(datafilep, "\n");

   fflush(datafilep);
   }


#endif
