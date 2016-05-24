#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include "halo_file.h"

using namespace std;

void read_mock_text(const char filename[], vector<ParticleData>& v)
{
  ParticleData p;

  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    cerr << "Unable to open fof file: " << filename << endl;
    throw 1;
  }

  char buf[256];
  while(fgets(buf, 255, fp)) {
    if(buf[0] == '#')
      continue;
    int ret= sscanf(buf, "%e %e %e %e %e %e", 
	     p.x, p.x+1, p.x+2, p.v, p.v+1, p.v+2);
    if(ret == 6)
      v.push_back(p);
  }
}

void read_fof_text(const char filename[], vector<ParticleData>& v, 
		   const float m, const float logMmin, const float logMmax)
{
  ParticleData p;
  const float M_min= pow(10.0f, logMmin);
  const float M_max= pow(10.0f, logMmax);

  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    cerr << "Unable to open fof file: " << filename << endl;
    throw 1;
  }

  int nfof;
  char buf[256];
  while(fgets(buf, 255, fp)) {
    int ret= sscanf(buf, "%d %e %e %e %e %e %e", 
	     &nfof, p.x, p.x+1, p.x+2, p.v, p.v+1, p.v+2);
    if(ret == 7 && M_min <= m*nfof && m*nfof < M_max)
      v.push_back(p);
  }
}

void read_subhalo_text(const char filename[], vector<ParticleData>& v,
		       const float logMmin, const float logMmax)
{
  ParticleData p;
  
  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    cerr << "Unable to open fof file: " << filename << endl;
    throw 1;
  }

  const float Mmin= pow(10.0, logMmin);
  const float Mmax= pow(10.0, logMmax);
  
  char buf[256];
  while(fgets(buf, 255, fp)) {
    if(buf[0] == '#')
      continue;

    float m;
    int ret= sscanf(buf, "%e %e %e %e %e %e %e", 
		    &m, p.x, p.x+1, p.x+2, p.v, p.v+1, p.v+2);
    assert(ret == 7);

    if(Mmin <= m && m < Mmax)
      v.push_back(p);
  }
}

float read_fof_binary(const char filename[], vector<ParticleData>& v, const int nfof_min)
{
  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    cerr << "Unable to open fof file: " << filename << endl;
    throw "FOF File Error";
  }

  int id;
  float ll;
  int nhalo;
  float boxsize;

  fread(&id, sizeof(int), 1, fp);
  assert(id == 604);
  fread(&ll, sizeof(float), 1, fp);
  fread(&boxsize, sizeof(float), 1, fp);
  fread(&nhalo, sizeof(int), 1, fp); assert(nhalo);
  //cout << "# nhalo in file " << nhalo << endl;

  /*
  const double n = 0.25e6*pow(*boxsize/1000.0, 3.0);

  //const int n_start= (int)(nthquater*n);
  //const int n_read= (int) round(n);
  //const int n_end= (int)(nthquater*n) + (int) round(n);
  if(n_end > nhalo) {
    cerr << "Not enough particles for " << nthquater << "th quater: "
         << n_end << " > " << nhalo << " =nhalo";
    throw "FOF Read Error";
  }

  pv->reserve(n_read);
  */

  //const int block_size= sizeof(int)+6*sizeof(float);
  //fseek(fp, n_start*block_size, SEEK_CUR);

  int nfof;
  ParticleData p;

  for(int j=0; j<nhalo; ++j) {
    fread(&nfof, sizeof(int), 1, fp);
    fread(p.x, sizeof(float), 3, fp);
    fread(p.v, sizeof(float), 3, fp);
    if(nfof >= nfof_min)
      v.push_back(p);
  }
  
  //fseek(fp, (nhalo-n_start-n_read)*block_size, SEEK_CUR);
  int nhalo_check= 0;
  int ret= fread(&nhalo_check, sizeof(int), 1, fp); assert(ret == 1);
  assert(nhalo_check == nhalo);

  fclose(fp);

  // printf("# rank %d - %d; total %d\n", n_start, n_end, nhalo);

  return boxsize;
}

void read_subsample_binary(const char filename[], vector<ParticleData>& v, float* const boxsize)
{
  FILE* fp= fopen(filename, "r");
  if(fp == 0) {
    cerr << "Unable to open fof file: " << filename << endl;
    throw 1;
  }

  float header[6];
  int np=0, np_check=0;
  fread(header, sizeof(float), 6, fp);
  *boxsize= header[0];

  fread(&np, sizeof(int), 1, fp);
  v.reserve(np);

  ParticleData p;
  for(int i=0; i<np; ++i) {
    fread(p.x, sizeof(float), 3, fp);
    fread(p.v, sizeof(float), 3, fp);
    v.push_back(p);
  }
  
  fread(&np_check, sizeof(int), 1, fp);
  assert(np == np_check);

  fclose(fp);
}
