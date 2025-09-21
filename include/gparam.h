#ifndef GPARAM_H
#define GPARAM_H

#include<stdio.h>
#include<time.h>

#include"macro.h"

typedef struct GParam {
  // lattice dimensions
  int d_size[STDIM];

  // simulation parameters
  double d_J;
  double d_K;
  double d_masssq;
  double d_chemical;


  // simulation details
  int d_sample;
  int d_thermal;
  int d_overrelax;
  int d_measevery;
  // initialization & saving
  int d_start;
  int d_saveconf_back_every;

  // for metropolis
  double d_epsilon_metro_link;
  double d_epsilon_metro_site;

  // for quenched averages if GAUGE_FIX is defined
  double d_quench_gamma;
  double d_quench_epsilon_metro;
  int d_quench_sample;
  int d_quench_thermal;
  int d_quench_measevery;
  int d_quench_overrelax;

  //for smearing and blocking
  int smearing_steps;
  double alpha;
  double blockcoeff;
  int numblock;

  // output file names
  char d_conf_file[STD_STRING_LENGTH];
  char d_data_file[STD_STRING_LENGTH];
  char d_log_file[STD_STRING_LENGTH];

  // random seed
  unsigned int d_randseed;

  // derived constants
  long d_volume;           // total volume
  double d_inv_vol;        // 1 / tot. volume
} GParam;


void remove_white_line_and_comments(FILE *input);
void readinput(char *in_file, GParam *param);
void init_derived_constants(GParam *param);

void init_data_file(FILE **dataf, GParam const * const param, char * fsufix);

void print_parameters(GParam const * const param,
                      time_t time_start,
                      time_t time_end,
                      double acc_site,
                      double acc_link,
                      double acc_link_big,
					  double quench_acc);

#endif
