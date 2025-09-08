#ifndef MONOPOLES_C
#define MONOPOLES_C

#include"../include/macro.h"

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

#include"../include/conf.h"
#include"../include/geometry.h"
#include"../include/gparam.h"
#include"../include/random.h"

void space_polyakov(Conf const * const GC,
                    Geometry const * const geo,
                    GParam const * const param,
                    long r,
                    double complex * polyakov_op_t)
   {
   int i;
   double complex poly;
   double complex poly2_ev,poly2_odd,poly2tmpl,poly2tmpr;
   poly=poly2_ev=poly2_odd=poly2tmpl=poly2tmpr=1.0+0.0*I;
   long int r2;
   r2=r;

   for(i=0; i<param->d_size[1];i++)
      {
      poly*=GC->lambda[r][1];
      r=nnp(geo,r,1);
      }
   for(i=0; i<param->d_size[1];i++)
      {
      poly2tmpl+=poly*plaquette_single(GC, geo, r2, 1, 2);
      poly2tmpr+=poly*plaquette_single(GC, geo, nnm(geo,r2,2), 2, 1);
      r2=nnp(geo,r2,1);
      }
   poly2_ev=(poly2tmpl+poly2tmpr)/(2*param->d_size[1]);
   poly2_odd=(poly2tmpl-poly2tmpr)/(2*param->d_size[1]);

   polyakov_op_t[0]+=poly;
   polyakov_op_t[1]+=poly2_ev;
   polyakov_op_t[2]+=poly2_odd;
   }

void polyakov_averaged(Conf const * const GC,
                       Geometry const * const geo,
                       GParam const * const param,
                       double complex * polyakov_op_t,
                       long r)
   {
   int i;
   long r2;

   r2=r;
   polyakov_op_t[0]=0;
   polyakov_op_t[1]=0;
   polyakov_op_t[2]=0;
   for(i=0; i<param->d_size[2]; i++)
      {
	space_polyakov(GC,geo,param,r2,polyakov_op_t);
        r2=nnp(geo,r2,2);
      }
   polyakov_op_t[0]/=param->d_size[2];
   polyakov_op_t[1]/=param->d_size[2];
   polyakov_op_t[2]/=param->d_size[2];

   }

void polyakov_time_sliced(Conf const * const GC,
                          Geometry const * const geo,
                          GParam const * const param,
                          double complex ** poly,
                          int ind)
   {
   long r;
   int t;

   double complex *polyakov_op;
   int err;
   int ops;
   ops=3;// because we have a basis of three operators for the polyakov
   err=posix_memalign((void**) &(polyakov_op),(size_t) DOUBLE_ALIGN, (size_t) ops * sizeof(double complex ));
   if(err!=0)
      {
      fprintf(stderr, "Problems in allocating the polyakov correlators! (%s, %d)\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }

   //fprintf(stdout, "(%s, %d)\n",  __FILE__, __LINE__);
   r=0;
   polyakov_op[1]=0;
   for(t=0; t<param->d_size[0]; t++)
      {
	 polyakov_averaged(GC,geo,param,polyakov_op,r);
         poly[ind][t]=polyakov_op[0];
         poly[ind+1][t]=polyakov_op[1];
         poly[ind+2][t]=polyakov_op[2];
         r=nnp(geo,r,0);
      }

   free(polyakov_op);
   }

void measure_polyakov_corr(GParam const * const param,
                           double complex ** poly,
                           FILE * datafilep,
                           int ops)
      {
   int t1, t2, t ;
   int i,j;
   double corr_repoly, corr_impoly;
   //fprintf(stdout, "(%s, %d)\n",  __FILE__, __LINE__);
   for(t = 0; t<param->d_size[0]/2; t++)
      for(i=0; i<ops; i++)
         for(j=0; j<=i; j++)
        {
        corr_repoly = 0.0;
        corr_impoly = 0.0;
        for(t1 = 0; t1<param->d_size[0]; t1++)
           {
           t2=(t1+t) % param->d_size[0];
           corr_repoly += creal(poly[i][t2]*poly[j][t1]);
           corr_impoly += cimag(poly[i][t2]*poly[j][t1]);
           }
        corr_repoly/=(double) param->d_size[0];
        corr_impoly/=(double) param->d_size[0];

        fprintf(datafilep, "%.12f %.12f ", corr_repoly, corr_impoly);

        }

   fprintf(datafilep, "\n");
   fflush(datafilep);
   }

void staples_wilson_no_time(Conf const * const GC,
                                Geometry const * const geo,
                                long r,
                                int i,
                                double complex * U)
  {
  int j, l;
  long k;
  double complex link1, link2, link3, link12, stap;

  for(l=i+1; l< i + STDIM; l++)
     {
     j = (l % STDIM);

     if(j!=0)
       {

//
//       i ^
//         |   (1)
//         +----->-----+
//         |           |
//                     |
//         |           V (2)
//                     |
//         |           |
//         +-----<-----+-->   j
//       r     (3)
//

       link1=GC->lambda[nnp(geo, r, i)][j];  // link1 = (1)
       link2= GC->lambda[nnp(geo, r, j)][i];  // link2 = (2)
       link3= GC->lambda[r][j];               // link3 = (3)

       link12=link1*conj(link2);  // link12=link1*link2^{dag}
       stap=link12*conj(link3);   // stap=link12*stap^{dag} (typo to link3)

       *U=stap;

//
//       i ^
//         |   (1)
//         |----<------+
//         |           |
//         |
//     (2) V           |
//         |
//         |           |
//         +------>----+--->j
//        k     (3)    r
//

       k=nnm(geo, r, j);

       link1=GC->lambda[nnp(geo, k, i)][j];  // link1 = (1)
       link2=GC->lambda[k][i];               // link2 = (2)
       link3=GC->lambda[k][j];               // link3 = (3)

       link12=conj(link1*link2); // link12=link1^{dag}*link2^{dag}
       stap=link12*link3;        // stap=link12*link3

       *U+=stap;
       }
     }
   }


// perform smearing on spatial links
void spatial_smearing(Conf const * const GC,
                      Geometry const * const geo,
                      GParam const * const param)
  {
  int i;// step;
  long r;
  double complex U;
  Conf staple_GC;
  U=0;

  //init_conf(&staple_GC,param);
  copy_gauge_conf(&staple_GC,GC,param);
  //fprintf(stdout, "(%s, %d)\n",  __FILE__, __LINE__);
  //fprintf(stdout,"Staple initialized \n");
  //for(step=0; step<1; step++)//step<param->smearing_steps; step++)
  //   {
     //fprintf(stdout, "(%s, %d)\n",  __FILE__, __LINE__);
     for(r = 0; r < param->d_volume; r++)
        {
        for(i = 1; i < STDIM; i++)
           {
           staples_wilson_no_time(GC, geo, r, i, &U);
           staple_GC.lambda[r][i]=U;
           }
        }
     //fprintf(stdout, "(%s, %d)\n",  __FILE__, __LINE__);
     for(r = 0; r < param->d_volume; r++)
        {
        for(i = 1; i < STDIM; i++)
           {
           staple_GC.lambda[r][i]*=param->alpha;
           GC->lambda[r][i]+=conj(staple_GC.lambda[r][i]);
           GC->lambda[r][i]/=sqrt((double) (GC->lambda[r][i]*conj(GC->lambda[r][i]))); //Reunitarize
           }
        //}
     }

  //fprintf(stdout, "(%s, %d)\n",  __FILE__, __LINE__);
  free_conf_gauge(&staple_GC, param);
  }


// compute blocked link for a given site
void spatialblocking_singlesite(Conf const * const GC,
                                Geometry const * const geo,
                                GParam const * const param,
                                long r,
                                int i,
                                double complex * U)
  {
  int j;
  long k, k1;

  double complex link1, link2, link3, link4, stap;

  #ifdef DEBUG
  if(r >= param->d_volume)
    {
    fprintf(stderr, "r too large: %ld >= %ld (%s, %d)\n", r, param->d_volume, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i >= STDIM)
    {
    fprintf(stderr, "i too large: i=%d >= %d (%s, %d)\n", i, STDIM, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  if(i == 0)
    {
    fprintf(stderr, "time direction selected: i=%d (%s, %d)\n", i, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  j = (i==1) ? 2 : 1; //Check the dir perpendicular to time and chosen link
//
//       i ^
//         |    (1)
//         +----->-----+
//         |           |
//         |           V (2)
//         |           |
//         |           |
//       k +-----------+
//         |           |
//         |           |
//         |           V (3)
//         |           |
//         +-----<-----+-->   j
//       r     (4)
//

    k=nnp(geo, r, i);
    link1= GC->lambda[nnp(geo, k, i)][j];  // link1 = (1)
    link2= GC->lambda[nnp(geo, k, j)][i];  // link2 = (2)
    link3= GC->lambda[nnp(geo, r, j)][i];  // link3 = (3)
    link4= GC->lambda[r][j];               // link3 = (4)

    stap=link1*conj(link2)*conj(link3)*conj(link4);

//
//       i ^
//         |   (1)
//         +----<------+
//         |           |
//     (2) V           |
//         |           |
//      k1 +-----------+
//         |           |
//     (3) V           |
//         |           |
//         +------>----+--->j
//        k     (4)    r
//

   k=nnm(geo,r,j);
   k1=nnp(geo,k,i);

   link1= GC->lambda[nnp(geo, k1, i)][j];  // link1 = (1)
   link2= GC->lambda[k1][i];               // link2 = (2)
   link3= GC->lambda[k][i];                // link3 = (3)
   link4= GC->lambda[k][j];                // link4 = (4)

   stap+=conj(link1)*conj(link2)*conj(link3)*link4;

   *U=(1-param->blockcoeff)*GC->lambda[r][i]*GC->lambda[nnp(geo, r, i)][i];

   stap*=param->blockcoeff/6;
   *U+=conj(stap);


   *U/=sqrt((double) (*U*conj(*U)));
   }


// create spatially blocked configurations (i.e. L_spatial->L_spatial/2)
void init_spatial_blocked_conf(Conf *blockGC,
                               GParam * blockparam,
                               Conf const * const GC,
                               Geometry const * const geo,
                               GParam const * const param)
  {

  long rb, r;
  int blockcart[STDIM], cart[STDIM];
  int i, mu, err;
  double complex U;


  long blockvol;

  *blockparam=*param;
  for(i=1; i<STDIM; i++)
     {
     if(param->d_size[i] % 2 != 0)
       {
       fprintf(stderr, "Problem with spatial size not even: %d ! (%s, %d)\n", param->d_size[i], __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  blockvol=param->d_size[0];
  for(i=1; i<STDIM;i++)
         {
         blockparam->d_size[i]=param->d_size[i]/2;
         blockvol*=blockparam->d_size[i];
         }
      blockparam->d_volume=blockvol;
      blockparam->d_inv_vol=1.0/((double) blockparam->d_volume);


  // allocate the lattice
  err=posix_memalign((void**)&(blockGC->lambda), (size_t) DOUBLE_ALIGN, (size_t) blockparam->d_volume * sizeof(double complex *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  for(r=0; r<(blockparam->d_volume); r++)
     {
     err=posix_memalign((void**)&(blockGC->lambda[r]), (size_t) DOUBLE_ALIGN, (size_t) STDIM * sizeof(double complex));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }


  // initialize GC

  for(rb=0; rb<(blockparam->d_volume); rb++)
     {
     si_to_cart(blockcart, rb, blockparam);
     cart[0]=blockcart[0];
     for(i=1; i<STDIM; i++)
        {
        cart[i]=2*blockcart[i];
        }
     r=cart_to_si(cart, param);

     blockGC->lambda[rb][0]=GC->lambda[r][0];
     for(mu=1; mu<STDIM; mu++)
        {
        // this is the real point where the blocking is performed
        spatialblocking_singlesite(GC, geo, param, r, mu, &U);
        blockGC->lambda[rb][mu]=U;
        }
     }
  (void) geo;
  blockGC->update_index=GC->update_index;
  }


void block_measure_operators(Conf const * const GC,
                                 Geometry const *const geo,
                                 GParam const * const param,
                                 double complex ** polyakov_loop)
   {
   Conf blockGC;
   Conf blockGC2;
   Conf SmearedGC;
   Geometry blockgeo;
   Geometry blockgeo2;
   GParam blockparam;
   GParam blockparam2;
   int ind_op;

   //Polyakov correlator bare
   ind_op=0;
   polyakov_time_sliced(GC,geo,param,polyakov_loop,ind_op);
   ind_op+=3;

   polyakov_time_sliced(GC,geo,param,polyakov_loop,ind_op);
   //Create smeared copy
   copy_gauge_conf(&SmearedGC,GC,param);


   int n;
   for (n=0; n<param->smearing_steps;n++)
      {
      //fprintf(stdout, "(%s, %d)\n",  __FILE__, __LINE__);
      spatial_smearing(&SmearedGC,geo,param);
      //fprintf(stdout, "(%s, %d)\n",  __FILE__, __LINE__);
      polyakov_time_sliced(&SmearedGC,geo,param,polyakov_loop,ind_op);
      ind_op+=3;
      }

   //fprintf(stdout, "(%s, %d)\n",  __FILE__, __LINE__);
   if (param->numblock >0 ){
      //Create the first blocked lattice
      init_spatial_blocked_conf(&blockGC,&blockparam,&SmearedGC, geo,param);
      init_geometry(&blockgeo, &blockparam);
      
      //fprintf(stdout,"Initialized 0 \n");

      //Create a copy. A copy is necessary when blocking again
      blockparam2=blockparam;
      copy_gauge_conf(&blockGC2,&blockGC,&blockparam);
      init_geometry(&blockgeo2, &blockparam2);
   }
     

   free_conf_gauge(&SmearedGC,param);  
   
   int j;
   int k;
   for(j=0; j<param->numblock; j++)
      {
      for (k=0; k<param->smearing_steps; k++)
      {
      polyakov_time_sliced(&blockGC2,&blockgeo2,&blockparam2,polyakov_loop,ind_op);
      ind_op+=3;
      spatial_smearing(&blockGC2,&blockgeo2,&blockparam2);
      }
      polyakov_time_sliced(&blockGC2,&blockgeo2,&blockparam2,polyakov_loop,ind_op);
      ind_op+=3;

      if(j<(param->numblock-1))
         {
         //Free the previous blocked lattice
         free_conf_gauge(&blockGC,&blockparam);
         free_geometry(&blockgeo, &blockparam);

         //Create next blocked lattice from the smeared lattice
         init_spatial_blocked_conf(&blockGC,&blockparam,&blockGC2,&blockgeo2,&blockparam2);
         init_geometry(&blockgeo, &blockparam);

         //Free the last smeared lattice
         free_conf_gauge(&blockGC2,&blockparam2);
         free_geometry(&blockgeo2, &blockparam2);

         //Create the copy which will be smeared
         blockparam2=blockparam;
         copy_gauge_conf(&blockGC2,&blockGC,&blockparam);
         init_geometry(&blockgeo2, &blockparam2);

         }
      }
      //fprintf(stdout,"Delete \n");
      //Delete the blocked lattice
      free_conf_gauge(&blockGC2, &blockparam2);
      free_geometry(&blockgeo2, &blockparam2);

      //Delete the copy
       free_conf_gauge(&blockGC, &blockparam);
       free_geometry(&blockgeo, &blockparam);
}


void real_main(char *in_file)
    {
    Conf GC;
    Geometry geo;
    GParam param;

    long count;
    FILE *datafilep;

    time_t time1, time2;
    double acc_link, acc_site, acc_link_big;
    double acc_link_local, acc_site_local, acc_link_big_local;
    // read input file
    readinput(in_file, &param);

    // initialize random generator
    initrand(param.d_randseed);

    // open data_file
    init_data_file(&datafilep, &param);

    // initialize geometry
    init_geometry(&geo, &param);

    // initialize configuration
    init_conf(&GC, &param);

    // acceptance
    acc_link=0.0;
    acc_link_big=0.0;
    acc_site=0.0;

    // montecarlo
    time(&time1);

    double complex **polyakov_loop;
    int err;
    int ops;
    ops=((param.numblock+1)*(param.smearing_steps+1))*3;// +1 because of the unsmeared *3 because we have a basis of three operators for the polyakov
    err=posix_memalign((void**) &(polyakov_loop),(size_t) DOUBLE_ALIGN, (size_t) (ops) * sizeof(double complex *));
    if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the polyakov correlators! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
    int n;
    for(n=0; n<ops;n++)
       {
       err=posix_memalign((void**) &(polyakov_loop[n]), (size_t) DOUBLE_ALIGN, (size_t) param.d_size[0] * sizeof(double complex));
       if(err!=0)
          {
          fprintf(stderr, "Problems in allocating the polyakov correlators! (%s, %d)\n", __FILE__, __LINE__);
          exit(EXIT_FAILURE);
          }
       }


    /*double complex **glueball;
    ops=(param.numblock*param.smearing_steps+1)*1;// *1 because we have a basis of one for the glueballs
    err=posix_memalign((void**) &(glueball),(size_t) DOUBLE_ALIGN, (size_t) (ops) * sizeof(double complex *));
    if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the polyakov correlators! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }

    for(n=0; n<ops;n++)
       {
       err=posix_memalign((void**) &(glueball[n]), (size_t) DOUBLE_ALIGN, (size_t) param.d_size[0] * sizeof(double complex));
       if(err!=0)
          {
          fprintf(stderr, "Problems in allocating the polyakov correlators! (%s, %d)\n", __FILE__, __LINE__);
          exit(EXIT_FAILURE);
          }
       }
    (void) glueball;*/
    // count starts from 1 to avoid problems using %
    for(count=1; count < param.d_sample + 1; count++)
       {
       update(&GC, &geo, &param, &acc_site_local, &acc_link_local, &acc_link_big_local);
       if(count>param.d_thermal)
         {
         acc_site+=acc_site_local;
         acc_link+=acc_link_local;
         acc_link_big+=acc_link_big_local;
         }

       if(count<param.d_thermal)
         {
         if(acc_site_local>0.33)
           {
           param.d_epsilon_metro_site*=1.1;

           if(param.d_epsilon_metro_site>3.0)
             {
             param.d_epsilon_metro_site=3.0;
             }
           }
         else
           {
           param.d_epsilon_metro_site*=0.9;
           }

         if(acc_link_local>0.33)
           {
           param.d_epsilon_metro_link*=1.1;
      if(param.d_epsilon_metro_link>1.0)
        {
          param.d_epsilon_metro_link=1.0;
        }
           }
         else
           {
           param.d_epsilon_metro_link*=0.9;
           }
         }

       if(count % param.d_measevery ==0 && count > param.d_thermal)
         {
          double plaq;
          plaq=plaquette(&GC,&geo,& param);
          fprintf(datafilep, "%.12f ", plaq);

          block_measure_operators(&GC,&geo,&param,polyakov_loop);
          measure_polyakov_corr(&param,polyakov_loop,datafilep,ops);
          //measure_glueball_corr(&param,glueball,datafilep);
         }

       // save configuration for backup
       if(param.d_saveconf_back_every!=0)
         {
         if(count % param.d_saveconf_back_every == 0 )
           {
           // simple
           write_conf_on_file(&GC, &param);

           // backup copy
           write_conf_on_file_back(&GC, &param);
           }
         }
       }
    time(&time2);

    acc_site/=(double)(param.d_sample-param.d_thermal);
    acc_link/=(double)(param.d_sample-param.d_thermal);
    acc_link_big/=(double)(param.d_sample-param.d_thermal);



    // close data file
    fclose(datafilep);

    // save configuration
    if(param.d_saveconf_back_every!=0)
      {
      write_conf_on_file(&GC, &param);
      }

    print_parameters(&param, time1, time2, acc_site, acc_link, acc_link_big,0);


    // free configuration
    free_conf(&GC, &param);

    // free geometry
    free_geometry(&geo, &param);

    //free polyakov loops;
    for(n=0; n<ops; n++)
       {
       free(polyakov_loop[n]);
       }
    free(polyakov_loop);
    }


void print_template_input(void)
  {
  FILE *fp;

  fp=fopen("template_input.in", "w");

  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file template_input.in (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    fprintf(fp, "size 4 4 4\n");
    fprintf(fp,"\n");
    fprintf(fp, "J 5.705\n");
    fprintf(fp, "K 2.0\n");
    fprintf(fp, "masssq 1.0\n");
    fprintf(fp,"\n");
    fprintf(fp, "sample    10\n");
    fprintf(fp, "thermal   0\n");
    fprintf(fp, "overrelax 5\n");
    fprintf(fp, "measevery 1\n");

    fprintf(fp, "smearing_steps 2\n");
    fprintf(fp, "alpha_smearing 0.9\n");
    fprintf(fp, "numblock 1\n");
    fprintf(fp,"block_coeff 0.9\n");

    fprintf(fp,"\n");
    fprintf(fp, "epsilon_metro_site 0.8\n");
    fprintf(fp, "epsilon_metro_link 0.8\n");
    fprintf(fp,"\n");
    fprintf(fp, "start                   0  # 0=ordered  1=random  2=from saved configuration\n");
    fprintf(fp, "saveconf_back_every     5  # if 0 does not save, else save backup configurations every ... updates\n");
    fprintf(fp, "\n");
    fprintf(fp, "#output files\n");
    fprintf(fp, "conf_file  conf.dat\n");
    fprintf(fp, "data_file  dati.dat\n");
    fprintf(fp, "log_file   log.dat\n");
    fprintf(fp, "\n");
    fprintf(fp, "randseed 0    #(0=time)\n");
    fclose(fp);
    }
  }


int main (int argc, char **argv)
    {
    char in_file[50];

    if(argc != 2)
      {
      printf("\nPackage %s version %s\n", PACKAGE_NAME, PACKAGE_VERSION);
      printf("Claudio Bonati %s\n", PACKAGE_BUGREPORT);
      printf("Usage: %s input_file\n\n", argv[0]);

      printf("Compilation details:\n");
      printf("\tCHARGE: %d\n", CHARGE);
      printf("\tN_f (number of flavours): %d\n", NFLAVOUR);
      printf("\tST_dim (space-time dimensionality): %d\n", STDIM);
      printf("\n");
      printf("\tINT_ALIGN: %s\n", QUOTEME(INT_ALIGN));
      printf("\tDOUBLE_ALIGN: %s\n", QUOTEME(DOUBLE_ALIGN));

      #ifdef DEBUG
        printf("\n\tDEBUG mode\n");
      #endif

      #ifdef BIG_LINK
   printf("\n\t BIG LINK \n");
      #endif

      #ifdef CSTAR_BC
        printf("\n\tC^* BOUNDARY CONDITIONS\n");
      #else
        printf("\n\tPERIODIC BOUNDARY CONDITIONS\n");
      #endif

      #ifdef LINKS_FIXED_TO_ONE
         printf("\n\tLINKS FIXED TO 1\n");
      #endif

      #ifdef TEMPORAL_GAUGE
         printf("\n\tTEMPORAL GAUGE\n");
      #endif

     #ifdef GAUGE_FIX
        printf("\n\tGAUGE_FIX mode (stocastic gauge fixing)\n");
      #endif

      #ifdef DEBUG_GAUGE_FIX
        printf("\n\tDEBUG_GAUGE_FIX mode (stocastic gauge fixing)\n");
      #endif


      printf("\n");

      #ifdef __INTEL_COMPILER
        printf("\tcompiled with icc\n");
      #elif defined(__clang__)
        printf("\tcompiled with clang\n");
      #elif defined( __GNUC__ )
        printf("\tcompiled with gcc version: %d.%d.%d\n",
                __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
      #endif

      print_template_input();


      return EXIT_SUCCESS;
      }
    else
      {
      if(strlen(argv[1]) >= STD_STRING_LENGTH)
        {
        fprintf(stderr, "File name too long. Increse STD_STRING_LENGTH in include/macro.h\n");
        }
      else
        {
        strcpy(in_file, argv[1]);
        }
      }

    real_main(in_file);

    return EXIT_SUCCESS;
    }

#endif
