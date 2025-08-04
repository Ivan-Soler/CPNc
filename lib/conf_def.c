#ifndef CONF_DEF_C
#define CONF_DEF_C

#include"../include/macro.h"

#include</opt/openssl/openssl-1.1.1w/include/openssl/md5.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include"../include/endianness.h"
#include"../include/conf.h"
#include"../include/endianness.h"
#include"../include/flavour_matrix.h"
#include"../include/geometry.h"
#include"../include/gparam.h"
#include"../include/random.h"


void init_conf(Conf *GC,
               GParam const * const param)
  {
  long r, j;
  double theta;
  int err;

  // allocate the lattice
  err=posix_memalign((void**) &(GC->lambda), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(double complex *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  for(r=0; r<(param->d_volume); r++)
     {
     err=posix_memalign((void**)&(GC->lambda[r]), (size_t) DOUBLE_ALIGN, (size_t )STDIM * sizeof(double complex));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

	#ifdef GAUGE_FIX
  	err=posix_memalign((void**) &(GC->gauge), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(double complex));
  	if(err!=0)
  	{
  		fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
  		exit(EXIT_FAILURE);
  	}
  	for (r=0; r<(param->d_volume); r++)
  	{
  		err=posix_memalign((void**)&(GC->gauge[r]), (size_t) DOUBLE_ALIGN, (size_t )STDIM * sizeof(double complex));
  		     if(err!=0)
  		       {
  		       fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
  		       exit(EXIT_FAILURE);
  		       }
  		     }
	#endif

  	//allocating magnetic charge
   err=posix_memalign((void**) &(GC->charge), (size_t) INT_ALIGN, (size_t) param->d_volume * sizeof(int));
   if(err!=0)
     {
     fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
     }



  err=posix_memalign((void**) &(GC->phi), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(Vec));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  err=posix_memalign((void**) &(GC->Qh), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(FMatrix));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  // initialize lattice
  if(param->d_start==0) // ordered start
    {
    Vec v1, v2;

    GC->update_index=0;

    one_Vec(&v1);

    for(r=0; r<(param->d_volume); r++)
       {
       rand_rot_Vec(&v2, &v1, 0.05);
       equal_Vec(&(GC->phi[r]), &v2);

       for(j=0; j<STDIM; j++)
          {
          theta=0.05*(2.0*casuale()-1.0);
          theta=0;
          GC->lambda[r][j]=cos(theta)+I*sin(theta);
          }
      #ifdef DEBUG
       j=0;
       theta=PI*(2.0*casuale()-1.0);
       GC->lambda[r][j]=cos(theta)+I*sin(theta);
      #endif

       }

    /*#ifdef DEBUG
    //we introduce monopole anti-monopole pairx  for checks
    Geometry geo;
    init_geometry(&geo, param);
    int i,j;
    int k;
    i=0;
    j=1;
    k=2;
    long rr;
    long r;
    double theta,delta;

    delta=0.1;
    theta=PI-delta/4;

    r=0;

    fprintf(stderr, "Monopole inserted at: ");
    fprintf(stderr, "%ld \n ", r);
    rr=nnm(&geo,r,k);
    fprintf(stderr, "Anti-monopole inserted at: %ld",rr);

    GC->lambda[r][i]=cos(theta)-I*sin(theta);  // (1)
    GC->lambda[nnp(&geo, r, i)][j]=cos(delta)-I*sin(delta); // (2)
    GC->lambda[nnp(&geo, r, j)][i]=cos(theta)+I*sin(theta); // (3)
    GC->lambda[r][j]=cos(delta)+I*sin(delta); // (4)

    fprintf(stderr, "\n ");
    free_geometry(&geo, param);

    #endif*/
    }

  if(param->d_start==1)  // random start
    {
    Vec v1;

    GC->update_index=0;

    for(r=0; r<(param->d_volume); r++)
       {
       rand_vec_Vec(&v1);
       equal_Vec(&(GC->phi[r]), &v1);

       for(j=0; j<STDIM; j++)
          {
          theta=PI*(2.0*casuale()-1.0);
          GC->lambda[r][j]=cos(theta)+I*sin(theta);
          }
       }
    }

  if(param->d_start==2) // initialize from stored conf
    {
    read_conf(GC, param);
    }

  Geometry geo;
  init_geometry(&geo, param);
  for(r=0; r<param->d_volume;r++)
     {
     GC->charge[r]=monpoles_cube(GC,&geo,r,param);
     }
  free_geometry(&geo,param);



  #ifdef LINKS_FIXED_TO_ONE
    for(r=0; r<(param->d_volume); r++)
       {
       for(j=0; j<STDIM; j++)
          {
          GC->lambda[r][j]=1.0;
          }
       }
  #endif

  #ifdef TEMPORAL_GAUGE
    for(r=0; r<(param->d_volume); r++)
       {
       GC->lambda[r][0]=1.0;
       }
  #endif
  }

void update_charge(Conf *GC,
      Geometry const * const geo,
      int new_charges[4],
      long r,
      int i)
 {
   long rr;
   int j,k;

   j=(i+1)%STDIM;
   k=(j+1)%STDIM;
   GC->charge[r]=new_charges[0];

   rr=nnm(geo,r,k);
   GC->charge[rr]=new_charges[1];

   rr=nnm(geo,nnm(geo,r,k),j);
   GC->charge[rr]=new_charges[2];

   rr=nnm(geo,r,j);
   GC->charge[rr]=new_charges[3];

  }

void update_charge_test(Conf *GC,
      Geometry const * const geo,
      GParam const * const param)
   {
   int r;
   for(r=0; r<param->d_volume;r++)
      {
      GC->charge[r]=monpoles_cube(GC,geo,r,param);
      }
  }

void restart_gauge_conf(Conf *GC,
               GParam const * const param)
{
	int r;
	double theta;
	{
	for(r=0; r<(param->d_volume); r++)
	   {
	    //theta=PI*(2.0*casuale()-1.0);
		theta=0.05*(2.0*casuale()-1.0);
	    GC->gauge[r]=cos(theta)+I*sin(theta);
	    //GC->gauge[r]=1;
		}
	 }
}

void read_conf(Conf *GC, GParam const * const param)
  {
  FILE *fp;
  int i, dimension, tmp_i;
  int err, mu;
  long r;
  char md5sum_new[2*MD5_DIGEST_LENGTH+1];
  char md5sum_old[2*MD5_DIGEST_LENGTH+1];

  fp=fopen(param->d_conf_file, "r"); // open the configuration file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else // read the txt header of the configuration
    {
    err=fscanf(fp, "%d", &dimension);
    if(err!=1)
      {
      fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    if(dimension != STDIM)
      {
      fprintf(stderr, "The space time dimension of the configuration (%d) does not coincide with the one of the global parameter (%d)\n",
              dimension, STDIM);
      exit(EXIT_FAILURE);
      }

    for(i=0; i<STDIM; i++)
       {
       err=fscanf(fp, "%d", &tmp_i);
       if(err!=1)
         {
         fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }
       if(tmp_i != param->d_size[i])
         {
         fprintf(stderr, "The size (%d) of the configuration lattice (%d) does not coincide with the one of the global parameter (%d)\n",i,tmp_i,param->d_size[i]);
         exit(EXIT_FAILURE);
         }
       }

    err=fscanf(fp, "%ld %s\n", &(GC->update_index), md5sum_old);
    if(err!=2)
      {
      fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }

    fclose(fp);
    }

  fp=fopen(param->d_conf_file, "rb"); // open the configuration file in binary
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    // read again the header
    err=0;
    while(err!='\n')
         {
         err=fgetc(fp);
         }

    for(r=0; r<param->d_volume; r++)
       {
       err=read_from_binary_file_bigen_Vec(fp, &(GC->phi[r]));
       if(err!=0)
            {
            fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            }

       for(mu=0; mu<STDIM; mu++)
          {
          err=read_from_binary_file_bigen_doublecomplex(fp, &(GC->lambda[r][mu]));
          if(err!=0)
            {
            fprintf(stderr, "Error in reading the file %s (%s, %d)\n", param->d_conf_file, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            }
          }
       }
    fclose(fp);

    // compute the new md5sum and check for consistency
    compute_md5sum_conf(md5sum_new, GC, param);
    if(strncmp(md5sum_old, md5sum_new, 2*MD5_DIGEST_LENGTH+1)!=0)
      {
      fprintf(stderr, "The computed md5sum %s does not match the stored %s for the file %s (%s, %d)\n", md5sum_new, md5sum_old, param->d_conf_file, __FILE__, __LINE__);
      exit(EXIT_FAILURE);
      }
    }
  }


void free_conf(Conf *GC, GParam const * const param)
  {
  long i;

  for(i=0; i<(param->d_volume); i++)
     {
     free(GC->lambda[i]);
     }

  free(GC->lambda);

  free(GC->phi);

  free(GC->Qh);

  free(GC->charge);

  }

void free_conf_gauge(Conf *GC, GParam const * const param)
  {
  long i;

  for(i=0; i<(param->d_volume); i++)
     {
     free(GC->lambda[i]);
     }

  free(GC->lambda);

  }


void equal_conf(Conf const * const GC, Conf *GC2,
					GParam const * const param)
{
	long r;
	int i;

	GC2->update_index=GC->update_index;

	for(r=0; r<param->d_volume; r++)
	{
		equal_Vec(&GC2->phi[r], &GC->phi[r]);
		GC2->gauge[r]=GC->gauge[r];

		for(i=0; i<STDIM; i++)
		{
			GC2->lambda[r][i]=GC->lambda[r][i];
		}
	}
}

void equal_gauge_conf(Conf *GC2,
                     GParam const * const param,
                     Conf const * const GC)
{
   long r;
   int i;

   GC2->update_index=GC->update_index;

   for(r=0; r<param->d_volume; r++)
   {
      for(i=0; i<STDIM; i++)
      {
         GC2->lambda[r][i]=GC->lambda[r][i];
      }
   }
}

void copy_gauge_conf(Conf *GC2,
               Conf const * const GC,
               GParam const * const param)
  {
  long r;
  int err;

  // allocate the lattice
  err=posix_memalign((void**) &(GC2->lambda), (size_t) DOUBLE_ALIGN, (size_t) param->d_volume * sizeof(double complex *));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  for(r=0; r<(param->d_volume); r++)
     {
     err=posix_memalign((void**)&(GC2->lambda[r]), (size_t) DOUBLE_ALIGN, (size_t )STDIM * sizeof(double complex));
     if(err!=0)
       {
       fprintf(stderr, "Problems in allocating the lattice! (%s, %d)\n", __FILE__, __LINE__);
       exit(EXIT_FAILURE);
       }
     }

  equal_gauge_conf(GC2,param,GC);

  }

// save a configuration in ILDG-like format
void write_conf_on_file_with_name(Conf const * const GC,
                                  GParam const * const param,
                                  char const * const namefile)
  {
  long r;
  int i, mu, err;
  char md5sum[2*MD5_DIGEST_LENGTH+1];
  FILE *fp;

  compute_md5sum_conf(md5sum, GC, param);

  fp=fopen(namefile, "w"); // open the configuration file
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", namefile, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    fprintf(fp, "%d ", STDIM);
    for(i=0; i<STDIM; i++)
       {
       fprintf(fp, "%d ", param->d_size[i]);
       }
    fprintf(fp, "%ld %s\n", GC->update_index, md5sum);
    }
  fclose(fp);

  fp=fopen(namefile, "ab"); // open the configuration file in binary mode
  if(fp==NULL)
    {
    fprintf(stderr, "Error in opening the file %s (%s, %d)\n", namefile, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  else
    {
    for(r=0; r<param->d_volume; r++)
       {
       err=print_on_binary_file_bigen_Vec(fp, &(GC->phi[r]) );
       if(err!=0)
         {
         fprintf(stderr, "Error in writing the file %s (%s, %d)\n", namefile, __FILE__, __LINE__);
         exit(EXIT_FAILURE);
         }

       for(mu=0; mu<STDIM; mu++)
          {
          err=print_on_binary_file_bigen_doublecomplex(fp, GC->lambda[r][mu] );
          if(err!=0)
            {
            fprintf(stderr, "Error in writing the file %s (%s, %d)\n", namefile, __FILE__, __LINE__);
            exit(EXIT_FAILURE);
            }
          }
       }
    fclose(fp);
    }
  }


void write_conf_on_file(Conf const * const GC, GParam const * const param)
  {
  write_conf_on_file_with_name(GC, param, param->d_conf_file);
  }


void write_conf_on_file_back(Conf const * const GC, GParam const * const param)
  {
  char name[STD_STRING_LENGTH], aux[STD_STRING_LENGTH];
  static int counter=0;

  strcpy(name, param->d_conf_file);
  if(counter==0)
    {
    sprintf(aux, "_back0");
    }
  else
    {
    sprintf(aux, "_back1");
    }
  strcat(name, aux);

  write_conf_on_file_with_name(GC, param, name);

  counter=1-counter;
  }


// compute the md5sum of the configuration and save it in res, that is a char[2*MD5_DIGEST_LENGTH]
void compute_md5sum_conf(char *res, Conf const * const GC, GParam const * const param)
  {
  MD5_CTX mdContext;
  unsigned char c[MD5_DIGEST_LENGTH];
  long r;
  double a, b;
  int mu, k;

  MD5_Init(&mdContext);
  for(r=0; r<param->d_volume; r++)
     {
     for(k=0; k<NFLAVOUR; k++)
        {
        a=creal((GC->phi[r]).comp[k]);
        b=cimag((GC->phi[r]).comp[k]);

        if(endian()==0) // little endian
          {
          SwapBytesDouble(&a);
          SwapBytesDouble(&b);
          }
        MD5_Update(&mdContext, &a, sizeof(double));
        MD5_Update(&mdContext, &b, sizeof(double));
        }

     for(mu=0; mu<STDIM; mu++)
        {
        a=creal(GC->lambda[r][mu]);
        b=cimag(GC->lambda[r][mu]);

        if(endian()==0) // little endian
          {
          SwapBytesDouble(&a);
          SwapBytesDouble(&b);
          }
        MD5_Update(&mdContext, &a, sizeof(double));
        MD5_Update(&mdContext, &b, sizeof(double));
        }
     }
  MD5_Final(c, &mdContext);

  for(k = 0; k < MD5_DIGEST_LENGTH; k++)
     {
     sprintf(&(res[2*k]), "%02x", c[k]);
     }
  }



#endif
