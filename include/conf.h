#ifndef CONF_H
#define CONF_H

#include"macro.h"

#include<complex.h>
#include</opt/openssl/openssl-1.1.1w/include/openssl/md5.h>
#include<stdio.h>

#include"flavour_matrix.h"
#include"gparam.h"
#include"geometry.h"
#include"vec.h"

typedef struct Conf {
  long update_index;

  double complex **lambda;    // [volume][STDIM]
  Vec *phi;                   // [volume]
  double complex *gauge;      // [volume] needed if GAUGE_FIX is defined

  FMatrix *Qh;                // [volume][Nflavors][Nflavors]q
  } Conf;

inline double complex chargepow(double complex x)
  {
  int i;
  double complex ris=1.0+0.0*I;

  if(CHARGE>=0)
    {
    for(i=0;i<CHARGE; i++)
       {
       ris*=x;
       }
    }
  else
    {
    for(i=0;i<(-CHARGE); i++)
       {
       ris*=x;
       }
    }

  return ris;
  }

// in conf_def.c
void init_conf(Conf *GC,
               GParam const * const param);
void restart_gauge_conf(Conf *GC,
               GParam const * const param);
void read_conf(Conf *GC,
               GParam const * const param);
void free_conf(Conf *GC,
               GParam const * const param);
void equal_conf(Conf const * const GC, Conf *GC2,
					GParam const * const param);
void write_conf_on_file_with_name(Conf const * const GC,
                                  GParam const * const param,
                                  char const * const namefile);
void write_conf_on_file(Conf const * const GC,
                        GParam const * const param);
void write_conf_on_file_back(Conf const * const GC,
                             GParam const * const param);
void compute_md5sum_conf(char *res,        // the lenght is 2*MD5_DIGEST_LENGTH
                         Conf const * const GC,
                         GParam const * const param);
void calcstaples_for_gauge(Conf *GC,
		Geometry const * const geo,
		long r,
		double complex *staple);

// in conf_update.c
void calcstaples_for_phi(Conf *GC,
                         Geometry const * const geo,
                         long r,
                         Vec *staple);
int metropolis_for_phi(Conf *GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       long r);
void overrelaxation_for_phi(Conf *GC,
                            Geometry const * const geo,
                            long r);
void overrelaxation_for_gauge(Conf *GC,
                            Geometry const * const geo,
							GParam const * const param,
                            long r);
double complex plaqstaples_for_link(Conf *GC,
                                    Geometry const * const geo,
                                    long r,
                                    int i);
int metropolis_for_link(Conf *GC,
                        Geometry const * const geo,
                        GParam const * const param,
                        long r,
                        int i);
int metropolis_for_link_big(Conf *GC,
                            Geometry const * const geo,
                            GParam const * const param,
                            long r,
                            int i);
int metropolis_for_gauge(Conf *GC,
                         Geometry const * const geo,
                         GParam const * const param,
                         long r);
void update(Conf * GC,
            Geometry const * const geo,
            GParam const * const param,
            double *acc_site,
            double *acc_link,
            double *acc_link_big);
void gauge_apply(Conf *GC,
                 Geometry const * const geo,
                 GParam const * const param);
void inv_gauge_apply(Conf *GC,
                 Geometry const * const geo,
                 GParam const * const param);
double glass_evolution_and_meas(Conf *GC,
                                Conf *GC2,
                                GParam *param,
                                Geometry const * const geo,
                                FILE *datafilep,
								FILE *datafilep2);
// in conf_meas.c
double plaquette_single(Conf const * const GC,
                        Geometry const * const geo,
                        long r,
                        int i,
                        int j);
double plaquette(Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param);
double higgs_interaction(Conf const * const GC,
                         Geometry const * const geo,
                         GParam const * const param);
double gauge_fixing_interaction(Conf const * const GC,
		Geometry const * const geo,
		GParam const * const param);
double realpartlink(Conf const * const GC,
                    GParam const * const param);
double imagpartlink(Conf const * const GC,
                    GParam const * const param);
double flux_plaquette(Conf const * const GC,
        Geometry const * const geo,
        long r,
        int i,
        int j);
long monpoles_cube(Conf const * const GC,
        Geometry const * const geo,
        long r,
        GParam const * const param);
long measure_monopoles(Conf const * const GC,
        Geometry const * const geo,
        GParam const * const param);
void compute_flavour_observables_tensor(Conf const * const GC,
                                        GParam const * const param,
                                        double *tildeG0,
                                        double *tildeGminp);
void compute_flavour_observables_vector(Conf const * const GC,
                                        GParam const * const param,
                                        double *tildeG0,
                                        double *tildeGminp);
void perform_measures(Conf *GC,
                      GParam const * const param,
                      Geometry const * const geo,
                      FILE *datafilep);
void perform_vec_measures_buffer(Conf *GC,
                                 GParam const * const param,
                                 Geometry const * const geo,
                                 double buffer[5]);
void perform_tensor_measures_buffer(Conf *GC,
                                 GParam const * const param,
                                 Geometry const * const geo,
                                 double buffer[5]);
void perform_overlap_measures_buffer(Conf const * const GC,
                                     Conf const * const GC2,
                                     GParam const * const param,
                                     Geometry const * const geo,
                                     double buffer[5]);
void perform_gauge_measures_buffer(Conf const * const GC,
                                     GParam const * const param,
                                     Geometry const * const geo,
                                     double buffer[5]);

#endif
