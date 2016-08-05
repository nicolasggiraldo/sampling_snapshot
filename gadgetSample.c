#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sort_ulong.h>
//#include <gsl/gsl_interp.h>


#define GSL_MAX             0x00000000FFFFFFFF
#define MY_RANDOM64_MAX_INT 0xFFFFFFFFFFFFFFFF
#define MY_RANDOM42_MAX_INT 0x0000FFFFFFFFFFFF



/////////////////////////////////
//* Data structure definition *//
/////////////////////////////////
// General id's
typedef unsigned int  uint;
typedef unsigned long ulong;

typedef unsigned int idtype;
//typedef unsigned long idtype;


struct globalVariables
{
  char *SNAP_BASE;       // Path of the GADGET binary
  int  NSNAPS;
  int  GADGET_VERSION;   // GADGET version of the snapshot
  int  PERCENTAGE_TAKEN; // Perecentage of the snapshot to take
  unsigned long NPART_TOTAL;
  unsigned long NPART_SELECTED;
}GV;

struct gadget_head
{
  uint   npart[6];
  double mass[6];
  double time;
  double redshift;
  int    flag_sfr;
  int    flag_feedback;
  int    npartTotal[6];
  int    flag_cooling;
  int    num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int    flag_age;
  int    flag_metals;
  int    nallHW[6];
  char   fill[256-6*4-6*8-2*8-2*4-6*4-2*4-4*8-2*4-6*4]; // Fills to 256 Bytes
  //char   fill[256-6*4-6*8-2*8-2*4-6*4-2*4-4*8]; // Fills to 256 Bytes
};


struct part
{
  float pos[3];
  float vel[3];
  //float mass;
  idtype id;
}*Part;



////////////////////////
//* Global variables *//
////////////////////////

//int Npart_Total;
struct gadget_head Gheader;
unsigned long *IDs=NULL; 

int read_parameters(char param_file_name[])
{
  
  FILE *cfg =NULL;  // Stream to the parameter (config) file
  int   len = 200;   // Len of the read parameter
  char *buf =NULL; // buf variables to be used to read strings variables
  char *buf1=NULL;
  char *buf2=NULL; 
  char *dumb=NULL;
  
  
  if( (cfg=fopen(param_file_name,"r"))==NULL )
    {
      printf("%s not found.\n", param_file_name);
      // Value -1 means there is an error loading the param file
      return -1;
    }
  
  buf  = (char *) malloc( len*sizeof(char) );
  buf1 = (char *) malloc( len*sizeof(char) );
  buf2 = (char *) malloc( len*sizeof(char) );
  
  /* Reading SNAP_BASE parameter */  
  do{dumb=fgets(buf, len, cfg);}while(dumb[0]=='#');
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 )
    {
      printf("No 'SNAP_BASE' setting in configuration file.\n");
      return -2;
    }
  else
    {
      GV.SNAP_BASE = strdup(buf2);
      printf("Snapshot base name: %s\n", GV.SNAP_BASE);
    }
  
  /* Reading NSNAPS parameter */
  do{dumb=fgets(buf, len, cfg);}while(dumb[0]=='#');
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 )
    {
      printf("No 'NSNAPS' setting in configuration file.\n");
      return -2;
    }
  else
    {
      GV.NSNAPS=atoi(buf2);
      if(GV.NSNAPS>0)
	{
	  printf("Number of snapshots: %d\n", GV.NSNAPS);
	}
      else
	{
	  printf("Invalid 'NSNAPS' setting in configuration file.\n");
	  return -2;
	}
    }
  
  /* Reading GADGET_VERSION parameter */
  do{dumb=fgets(buf, len, cfg);}while(dumb[0]=='#');
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 )
    {
      printf("No 'GADGET_VERSION' setting in configuration file.\n");
      return -2;
    }
  else
    {
      GV.GADGET_VERSION=atoi(buf2);
      if(GV.GADGET_VERSION==1 || GV.GADGET_VERSION==2)
	{
	  printf("GADGET VERSION: %d\n", GV.GADGET_VERSION);
	}
      else
	{
	  printf("Invalid 'GADGET VERSION' setting in configuration file.\n");
	  return -2;
	}
    }

  /* Reading NPARTICLES parameter */
  do{dumb=fgets(buf, len, cfg);}while(dumb[0]=='#');
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 )
    {
      printf("No 'NPARTICLES' setting in configuration file.\n");
      return -2;
    }
  else
    {
      GV.NPART_TOTAL = atoi(buf2);
      if(GV.NPART_TOTAL > 0)
	{
	  GV.NPART_TOTAL= GV.NPART_TOTAL*GV.NPART_TOTAL*GV.NPART_TOTAL;
	  printf("NPARTICLES: %lu\n", GV.NPART_TOTAL);
	}
      else
	{
	  printf("Invalid 'NPARTICLES' setting in configuration file.\n");
	  return -2;
	}
    }

  /* Reading PERCENTAGE_TAKEN parameter */
  do{dumb=fgets(buf, len, cfg);}while(dumb[0]=='#');
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 )
    {
      printf("No 'PERCENTAGE_TAKEN' setting in configuration file.\n");
      return -2;
    }
  else
    {
      GV.PERCENTAGE_TAKEN = (int)round(atof(buf2));
      if(0<GV.PERCENTAGE_TAKEN && GV.PERCENTAGE_TAKEN<100)
	{
	  printf("PERCENTAGE TAKEN: %d\n", GV.PERCENTAGE_TAKEN);
	}
      else
	{
	  printf("Invalid 'PERCENTAGE_TAKEN' setting in configuration file.\n");
	  return -2;
	}
    }
  
  GV.NPART_SELECTED = (unsigned long) llround( (GV.PERCENTAGE_TAKEN/100.0)*GV.NPART_TOTAL );
  
  if(dumb==NULL){}
  
  fclose(cfg);
  free(buf);
  free(buf1);
  free(buf2);
  
  return 0;
}




unsigned long my_rng_uniform_int64(const gsl_rng *r, unsigned long n)
{
  unsigned long random;
  unsigned long range = 0xFFFFFFFFFFFFFFFF;
  unsigned long scale;
  
  if(n > MY_RANDOM64_MAX_INT || n == 0)
    {
      printf("Error: Invalid n, either 0 or exceeds maximum value of generator\n");
      exit(0);
    }

  scale = range / n;

  do
    {
      random = 
	( (gsl_rng_get(r) <<  0) & 0x00000000FFFFFFFF ) | 
	( (gsl_rng_get(r) << 32) & 0xFFFFFFFF00000000 ) ;
      
      random /= scale;

    }while(random >= n);

  return random;
}

unsigned long my_rng_uniform_int42(const gsl_rng *r, unsigned long n)
{
  
  unsigned long random;
  unsigned long range = 0x0000FFFFFFFFFFFF;
  unsigned long scale;
  
  if(n > MY_RANDOM42_MAX_INT || n == 0)
    {
      printf("Error: Invalid n, either 0 or exceeds maximum value of generator\n");
      exit(0);
    }
  
  scale = range / n;
  
  do
    {
      random = 
	( (gsl_rng_get(r) <<  0) & 0x00000000FFFFFFFF ) | 
	( (gsl_rng_get(r) << 32) & 0x0000FFFF00000000 ) ;
      
      random /= scale;
      
    }while(random >= n);
  
  return random;
}

int makeRandomIDs(ulong *randomIDs, ulong Npart, ulong len)
{  
  long long i,j;
  gsl_rng *r=NULL;
  long seed = time(NULL) * getpid();
  ulong (*rng_int)(const gsl_rng *, unsigned long) = NULL;

  if(Npart>GSL_MAX)
    {
      rng_int =  my_rng_uniform_int64;  
      printf("Using 64 bit random number generator\n");
    }
  else
    {
      rng_int = gsl_rng_uniform_int;
      printf("Using 32 bit random number generator\n");
    }
  
  
  // Initialize generation routines
  gsl_rng_env_setup();
  
  // Memory allocation
  r = gsl_rng_alloc(gsl_rng_mt19937);
  if( r==NULL )
    {
      printf("Error: Attempt to allocate instance of random number generator failed\n");
      return 1;
    }
  
  // Set the seed for our generator
  gsl_rng_set(r, seed);
  
  for(i=0; i<len; i++)
    {
      randomIDs[i] = rng_int(r, Npart);
    }

  do
    {
      gsl_sort_ulong(randomIDs, 1, len);
      j=0;
      for(i=1; i<len; i++)
	{
	  if( randomIDs[i-1]==randomIDs[i] )
	    {
	      randomIDs[i] = rng_int(r, Npart);
	      j++;
	    }
	}
    }while(j>0);
  
  gsl_rng_free(r);

  return 0;
}



int read_head0(char *infile, int gadgetVersion)
{  
  int dummi,i;
  char label[4];
  FILE *fp_inp;
  int err;
  
  if( (fp_inp=fopen(infile,"r"))==NULL )
    {
      printf("read_gadget: cannot open %s\n",infile);
      exit(0);
    }
  
  
  if(gadgetVersion==2)
    { 
      err=fread(&dummi, sizeof(dummi), 1, fp_inp);
      err=fread(&label, sizeof(char),  4, fp_inp);
      err=fread(&dummi, sizeof(dummi), 1, fp_inp);
      err=fread(&dummi, sizeof(dummi), 1, fp_inp);
    }

  err=fread(&dummi,   sizeof(dummi),   1, fp_inp);
  err=fread(&Gheader, sizeof(Gheader), 1, fp_inp);
  err=fread(&dummi,   sizeof(dummi),   1, fp_inp);
  
  fclose(fp_inp);
  
  
  for(i=0; i<6; i++)
    printf(" * %d Particles of class %d\n",Gheader.npart[i],i); 
  
  printf("\n"); 
  
  for(i=0; i<6; i++)
    { 
      if((Gheader.npart[i] != 0) && (Gheader.mass[i] != 0.0))
	printf(" * The mass of each particle is %d es %g\n",i,Gheader.mass[i]);
      
      if((Gheader.npart[i] != 0) && (Gheader.mass[i] == 0.0))
	printf(" * There are individual mases for this particle set %d\n",i);
      
    }     
  
  printf("\n");
  
  printf(" * Frame's Time... %g\n", Gheader.time); 
  printf(" * Redshift... %g\n",     Gheader.redshift);
  printf(" * Flagsfr... %d\n",      Gheader.flag_sfr);
  printf(" * Flagfed... %d\n",      Gheader.flag_feedback);
  
  //printf("\n");
  
  //GV.NPART_TOTAL = 0;
  //for(i=0; i<6; i++)
  //{
  //printf(" * Header nall[%d] is: %d\n",i,Gheader.npartTotal[i]);
  //GV.NPART_TOTAL += Gheader.npartTotal[i];
  //}

  //GV.NPART_SELECTED = (unsigned long) llround( (GV.PERCENTAGE_TAKEN/100.0)*GV.NPART_TOTAL );
  
  printf("\n");
  
  printf(" * Flagcool...       %d\n",  Gheader.flag_cooling);
  printf(" * numfiles...       %d\n",  Gheader.num_files);
  printf(" * Boxsize...        %g\n",  Gheader.BoxSize);
  printf(" * Omega0...         %g\n",  Gheader.Omega0);
  printf(" * OmageLa...        %g\n",  Gheader.OmegaLambda);
  printf(" * Hubbleparam...    %g\n",  Gheader.HubbleParam);
  printf(" * NPART_TOTAL...    %lu\n", GV.NPART_TOTAL);
  printf(" * NPART_SELECTED... %lu\n", GV.NPART_SELECTED);
  
  Part = (struct part *) malloc((size_t) GV.NPART_SELECTED*sizeof(struct part));
  if(Part == NULL)
    {
      printf("No memory available for load dar particles (%s)\n",infile);
      exit(0);
    }
  
  printf("Using %f Mb in memory\n",(GV.NPART_SELECTED*sizeof(struct part))/(1020*1024.0));

  if(err){}

  return 0;
}

int binary_search(ulong *arr, ulong val, ulong len)
{
  long long first, last, middle;
  
  first = 0;
  last = len - 1;
  middle = (first+last)/2;

  while(first <= last)
    {
      if(arr[middle] < val)
	{
	  first = middle + 1;
	}
      else if(arr[middle] == val)
	{
	  // If the routine finds a match then return 1
	  return 1;
	}
      else
	{
	  last = middle - 1;
	}
      
      middle = (first + last)/2;
    }
  // If the routine does not find a match then return 0
  return 0;
}



int main(int argc, char *argv[])
{
  int dummi,s,k,err;
  idtype id;
  unsigned long global_acum,Npart_snap,i;
  float pos[3],vel[3],aux[3];
  struct gadget_head header1;
  char label[4], buf[200], outfile[200];
  ulong blksize_pos,blksize_vel,blksize_ids;
  
  FILE *fp_pos  =NULL;
  FILE *fp_vel  =NULL;
  FILE *fp_ids  =NULL;
  FILE *fp_head =NULL;
  FILE *pf_out  =NULL;
  /////////////////////////////
  //FILE *out_big,*out_sample;
  //double xmin,xmax;
  //out_big = fopen ("outputgadget.txt","w");
  //out_sample = fopen ("outputsample.txt","w");
  //xmin=200.0;
  //xmax=205.0;
  /////////////////////////////////////////////////

  //gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  
  
  
  //////////////////////////////
  //* READING PARAMETER FILE *//
  //////////////////////////////
  if(argc<2){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("%s: You must specify the name of the parameter file\n",argv[0]);
    printf("For example: %s pathOfFile/parameterFile.txt\n",argv[0]);
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }
  
  // Reading parameter file and verifying there is no error.
  switch( read_parameters(argv[1]) )
    {
    case -1 :
      printf("\n***********************************");
      printf("***********************************\n");
      printf("Error: Bad path to the parameter file.\n" );
      printf("***********************************");
      printf("***********************************\n\n");
      exit(0);
    case -2 :
      printf("\n***********************************");
      printf("***********************************\n");
      printf("Error: Bad settings in the parameter file.\n" );
      printf("***********************************");
      printf("***********************************\n\n");
      exit(0);
    }
  


  ///////////////////////////////
  //* READING GADGET SNAPSHOT *//
  ///////////////////////////////
  global_acum = 0;
  // Reading the snapshot through its parts
  for(s=0; s<GV.NSNAPS; s++)
    {
      if(GV.NSNAPS>1)
	sprintf(buf,"%s.%d",GV.SNAP_BASE,s);
      else
	sprintf(buf,"%s",GV.SNAP_BASE);
      
      /* For the first part we take the header in order to get 
	 global information of the snapshot */
      if(s==0)
	{  
	  // Reading the header of the first snapshot
	  read_head0(buf,GV.GADGET_VERSION);
	  
	  // Storing in memory the array with the random ID's to select
	  IDs = (ulong *) malloc(sizeof(ulong)*GV.NPART_SELECTED);
	  if(IDs==NULL)
	    {
	      printf("Error: Attempt to allocate randomIDs array failed\n");
	      exit(0);
	    }
	  
	  // Setting which ID's are going to be selected
	  makeRandomIDs(IDs,GV.NPART_TOTAL,GV.NPART_SELECTED);
	}
      
      // Reading the positions, velocities and ID's of the snapshot
      if( (fp_head=fopen(buf,"r"))==NULL )
	{
	  printf("read_gadget cannot open %s",buf);
	  exit(0);
	}
      
      fp_pos = fopen(buf,"r");
      fp_vel = fopen(buf,"r");
      fp_ids = fopen(buf,"r");
      
      
      /* reading the header for this sub snap */
      if(GV.GADGET_VERSION==2)
	{
	  err=fread(&dummi, sizeof(dummi), 1, fp_head);
	  err=fread(&label, sizeof(char),  4, fp_head);
	  err=fread(&dummi, sizeof(dummi), 1, fp_head);
	  err=fread(&dummi, sizeof(dummi), 1, fp_head);
	}
      
      err=fread(&dummi,   sizeof(dummi),   1, fp_head);
      err=fread(&header1, sizeof(header1), 1, fp_head);
      err=fread(&dummi,   sizeof(dummi),   1, fp_head);
      
      fclose(fp_head);
      
      Npart_snap = 0;
      for(k=0; k<6; k++)
	{
	  printf(" * Header nsnap_%d[%d] is: %d\n",s,k,header1.npart[k]);
	  Npart_snap += header1.npart[k];
	}
      
      
      //  Begining with the groups
      if(GV.GADGET_VERSION==2)
	{
	  blksize_pos =  9*sizeof(int) + 2*4*sizeof(char) + sizeof(header1);
	  blksize_vel = 14*sizeof(int) + 3*4*sizeof(char) + sizeof(header1) +   3*Npart_snap*sizeof(float);
	  blksize_ids = 19*sizeof(int) + 4*4*sizeof(char) + sizeof(header1) + 2*3*Npart_snap*sizeof(float);
	}
      else
	{
	  blksize_pos = 3*sizeof(int) + sizeof(header1);
	  blksize_vel = 5*sizeof(int) + sizeof(header1) +   3*Npart_snap*sizeof(float);
	  blksize_ids = 7*sizeof(int) + sizeof(header1) + 2*3*Npart_snap*sizeof(float);
	}
    
    fseek(fp_pos, blksize_pos, SEEK_SET);
    fseek(fp_vel, blksize_vel, SEEK_SET);
    fseek(fp_ids, blksize_ids, SEEK_SET);
    
    for(i=0; i<Npart_snap ;i++)
      {
	err=fread(&pos[0], sizeof(float),  3, fp_pos);
	err=fread(&vel[0], sizeof(float),  3, fp_vel);
	err=fread(&id,     sizeof(idtype), 1, fp_ids);

	//////////////////////////////////////////
	//if(xmin <= pos[0] && pos[0] <= xmax)
	//{
	//fprintf(out_big,"%f %f %f\n", pos[0], pos[1], pos[2]);
	//}
	//////////////////////////////////////////
	
	// binary search
	if( (binary_search(IDs,(ulong)id,GV.NPART_SELECTED))==1 )
	  {
	    Part[global_acum].pos[0] = pos[0];
	    Part[global_acum].pos[1] = pos[1];
	    Part[global_acum].pos[2] = pos[2];
	    
	    Part[global_acum].vel[0] = vel[0];
	    Part[global_acum].vel[1] = vel[1];
	    Part[global_acum].vel[2] = vel[2];
      
	    Part[global_acum].id = id;
      
	    global_acum++;

	    //////////////////////////////////////
	    //if(xmin <= pos[0] && pos[0] <= xmax)
	    //{
	    //fprintf(out_sample,"%f %f %f\n", pos[0], pos[1], pos[2]);
	    //}
	    ///////////////////////////////////////
	  }
      }
    
    printf("global_acum    = %lu\n", global_acum);
    printf("NPART_SELECTED = %lu\n", GV.NPART_SELECTED);
    
    fclose(fp_pos);
    fclose(fp_vel);
    fclose(fp_ids);
    
    printf("=====================================================\n");
    printf(" *-* End of reading from gadget file %s\n",buf);
    printf("=====================================================\n");
    
    }
  //////////////////////////////////////////////
  //fclose(out_big);
  //fclose(out_sample);
  /////////////////////////////////////////
  
  
  /////////////////////////////////////////////
  /* Writing new snapshot in GADGET 1 format */
  /////////////////////////////////////////////
  
  //Gheader.npart[1] = Gheader.npartTotal[1];
  Gheader.npart[1]      = global_acum;
  Gheader.npartTotal[1] = global_acum;
  //random = 
  //( (gsl_rng_get(r) <<  0) & 0x00000000FFFFFFFF ) | 
  //( (gsl_rng_get(r) << 32) & 0xFFFFFFFF00000000 ) ;
  Gheader.num_files     = 1;
  
  sprintf(outfile,"%s%s",GV.SNAP_BASE,"_partial_gad1");
  
  printf("=====================================================\n");
  printf(" *-* Writing in %s\n",outfile);
  printf("=====================================================\n");
  
  if((pf_out=fopen(outfile,"w")) == NULL) {
    printf("cannot open %s\n",outfile);
    exit(0);
  }

  dummi = 256;

  
  // INIT GROUP HEADER
  fwrite(&dummi,   sizeof(dummi),   1, pf_out);
  fwrite(&Gheader, sizeof(Gheader), 1, pf_out);
  fwrite(&dummi,   sizeof(dummi),   1, pf_out);
  // END GROUP  HEADER

  //dummi = Npart_Total*3*sizeof(float);
  dummi = global_acum*3*sizeof(float);
  
  // INIT GROUP POSITIONS
  fwrite(&dummi, sizeof(dummi), 1, pf_out); 
  for(i=0; i<global_acum; i++)
    {
      for(k=0; k<3; k++) 
	aux[k] = Part[i].pos[k];
      
      fwrite(aux,sizeof(float),3,pf_out);
    }
  fwrite(&dummi,sizeof(dummi),1,pf_out);
  // END GROUP POSITIONS
  
  // INIT GROUP VELOCITIES
  //fwrite(&dummi, sizeof(dummi), 1, pf_out); 
  //for(i=0; i<global_acum; i++)
  //{
  //for(k=0; k<3; k++) 
  //aux[k] = Part[i].vel[k];
  //fwrite(aux,sizeof(float),3,pf_out);
  //}
  //fwrite(&dummi,sizeof(dummi),1,pf_out); 
  //END GROUP VELOCITIES
  
  //dummi = global_acum*sizeof(uint);

  // INIT GROUP ID'S
  //fwrite(&dummi,sizeof(dummi),1,pf_out); 
  //for(i=0; i<global_acum; i++)
  //{
  //uint iduint;
  //iduint = (uint)i;
  //fwrite(&Part[i].id,sizeof(idtype),1,pf_out);
  //fwrite(&iduint,sizeof(uint),1,pf_out);
  //}
  //fwrite(&dummi,sizeof(dummi),1,pf_out); 
  // END GROUP ID'S

  // INIT GROUP MASSES
  fwrite(&dummi,sizeof(dummi),1,pf_out); 
  fwrite(&dummi,sizeof(dummi),1,pf_out); 
  // END GROUP MASSES



  if(err){}
  
  fclose(pf_out);
  free(IDs); 
  //gsl_interp_accel_free(acc);

  return 0;
}
