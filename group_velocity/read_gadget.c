//
// Reads GADGET snapshot and distribute to all nodes
//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
//#include <mpi.h>
#include "particle.h"
#include "msg.h"

typedef struct {
  int      np[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  unsigned int np_total[6];
  int      flag_cooling;
  int      num_files;
  double   boxsize;
  double   omega0;
  double   omega_lambda;
  double   hubble_param;
  int flag_stellarage;
  int flag_metals;
  unsigned int np_total_highword[6];
  int  flag_entropy_instead_u;
  char fill[60];
} GadgetHeader;

static GadgetHeader* find_files(const char filename[]);

static void check_separator(FILE* fp, int expected_value)
{
  int n= -1;
  fread(&n, sizeof(int), 1, fp);
  
  if(n != expected_value) {
    msg_abort("Error: Unable to read snapshot correctly. Separator %d != %d\n",
	      n, expected_value);
  }  
}

int read_gadget_snapshot(const char filename[], Particles* const particles)
{
  msg_printf(msg_verbose, "Reading Gadget snapshot %s...\n", filename);

  int numfiles=0;
  {
    GadgetHeader* header= find_files(filename);
  
    if(header == 0) {
      msg_abort("Unable to find snapshot %s\n", filename);
    }

    msg_printf(msg_verbose, "%d files per snapshot.\n", header->num_files);

    particles->np_total= (((long long) header->np_total_highword[1]) << 32) +
                           (long long) header->np_total[1];

    particles->p= malloc(sizeof(Particle)*particles->np_total);
    msg_printf(msg_verbose, "Allocating %ld Mbytes for %ld particles.\n",
	       sizeof(Particle)*particles->np_total/(1024*1024),
	       particles->np_total);
    
    assert(particles->p);

    numfiles= header->num_files;
    free(header);
  }

  Particle* const p= particles->p;


  GadgetHeader header;
  float vfac_back; // Gadget convention back to km/s

  size_t ipart, nread= 0;
  for(int num = 0; num<numfiles; num++) {
    int np_snapshot= 0;

    char filename_i[256];

    if(numfiles > 1)
      sprintf(filename_i, "%s.%d", filename, num);
    else
      sprintf(filename_i, "%s", filename);
    
    FILE* fp= fopen(filename_i, "r");
    if(fp == 0) {
      msg_abort("Error: unable to open snapshot file %s (read.c)\n", 
		filename_i);
    }
    
    check_separator(fp, 256);
    fread(&header, sizeof(GadgetHeader), 1, fp);
    check_separator(fp, 256);
    
    np_snapshot = header.np[1];   // This code only reads type 1 dark matter
    vfac_back= sqrt(header.time);
    
    msg_printf(msg_verbose, "Reading %d particles read from %s.\n", 
	       np_snapshot, filename_i);
    
    float x[3], v[3];
    const float boxsize= (float) header.boxsize;
    
    // position
    ipart= nread;
    
    check_separator(fp, sizeof(float)*3*np_snapshot);
    for(int i=0; i<np_snapshot; i++) {
      int ret= fread(x, sizeof(float), 3, fp); 
      assert(ret == 3);
      
      for(int k=0; k<3; k++) {
	if(x[k] < 0.0f) x[k] += boxsize;
	if(x[k] >= boxsize) x[k] -= boxsize;
      }
      
      p[ipart].x[0]= x[0];
      p[ipart].x[1]= x[1];
      p[ipart].x[2]= x[2];
      ipart++;
    }
    check_separator(fp, sizeof(float)*3*np_snapshot);
    
    // velocity
    ipart= nread;
    check_separator(fp, sizeof(float)*3*np_snapshot);
    for(int i=0; i<np_snapshot; i++) {
      int ret= fread(v, sizeof(float), 3, fp); 
      assert(ret == 3);
      
      p[ipart].v[0]= vfac_back*v[0];
      p[ipart].v[1]= vfac_back*v[1];
      p[ipart].v[2]= vfac_back*v[2];
      
      ipart++;
    }
    check_separator(fp, sizeof(float)*3*np_snapshot);
    nread += np_snapshot;
    
    fclose(fp);
  }

  particles->a= header.time;
  particles->boxsize= (float) header.boxsize;
  particles->omega_m= header.omega0;
  particles->h= header.hubble_param;
  particles->np_local= nread;

  msg_printf(msg_verbose, "%lu particles read.\n", particles->np_local);

  return 1;
}


GadgetHeader* find_files(const char filename[])
{
  GadgetHeader* const header= malloc(sizeof(GadgetHeader));
  assert(header);

  msg_printf(msg_debug, "Trying %s\n", filename);
  FILE* fp= fopen(filename, "r");
  if(fp) {
    check_separator(fp, 256);
    fread(header, sizeof(header), 1, fp);
    check_separator(fp, 256);
    fclose(fp);

    return header;
  }


  char filename_i[256];
  sprintf(filename_i, "%s.0", filename);
  msg_printf(msg_debug, "Trying %s\n", filename_i);
  
  fp= fopen(filename_i, "r");
  if(fp) {
    check_separator(fp, 256);
    fread(header, sizeof(GadgetHeader), 1, fp);
    check_separator(fp, 256);
    fclose(fp);
    
    return header;
  }

  return 0;
}
