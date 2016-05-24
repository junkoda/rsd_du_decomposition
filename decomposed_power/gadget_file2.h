// gadget_file1 ver 1.0 April 7, 2005

// read() converts velocities to peculiar km/s iff
//  one set_velocity_conversion(true) and hubble_param > 0; Feb 2008
// write() requires internal velocity (no conversion)

// read() converts positions to physical kpc iff
//  one set_physical_coordinate(true) and hubble_param > 0; Feb 2008

#ifndef GADGET_FILE1_H
#define GADGET_FILE1_H

#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cassert>

//#define LONGID 1

typedef int int4;

#ifdef LONGID
typedef unsigned long long partid_t;
#else
typedef unsigned int partid_t;
#endif

struct gadget_header1{
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
  //char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];
           /* fills to 256 Bytes */
};

template<class Tgas, class Tcdm> class gadget_file{
 private:
  gadget_header1 header;
  std::vector<Tgas> *v_gas;
  std::vector<Tcdm> *v_cdm[4];
  
  int    max(int *a);
  int    fread_seperator(FILE * const fp, int i);
  int    read1(const char * const filename);
  bool   convert_position, convert_velocity;
 public:
  gadget_file();

  int    read_header(const char * const filename);
  int    read(const char * const filename);
  int    read_one(const char filename[]);
  int    write(const char * const filename);
  
  void   clear_header();

  gadget_header1* get_header(){return &header;};
  int    get_np(int type){return header.np[type];};
  int    get_np_total(int type){return header.np_total[type];};
  double get_mass(int type){return header.mass[type];};
  double get_time(){return header.time;};
  double get_redshift(){return header.redshift;};
  double get_boxsize(){return header.boxsize;};

  void   set_gas(std::vector<Tgas> *p_v_gas){v_gas= p_v_gas;};
  void   set_cdm(std::vector<Tcdm> *p_v_cdm, int type){
           v_cdm[type-1]= p_v_cdm;}; // type= 1,2,3,4
  void   set_number_of_files(const int nfiles){
                                              header.num_files= nfiles;};
  void   set_time(const double time){header.time= time;};
  void   set_cosmology(const double redshift, const double boxsize,
			   const double omega0, const double omegalambda,
			   const double hubble_param);
                           // This also sets time= 1/(1+z)
  void   set_mass(const double mgas, const double mhalo,
		      const double mdisk=0, const double mbulge=0,
		      const double mstars=0);
  void   set_n_particles(const int ngas, const int nhalo,
			 const int ndisk=0, const int nbulge=0,
			 const int nstars=0);
  void   set_n_particles_total(const int ngas, const int nhalo,
			 const int ndisk=0, const int nbulge=0,
			 const int nstars=0);
  void   print();
  void   set_velocity_conversion(const bool flag){ convert_velocity= flag; }
  void   set_physical_coordinate(const bool flag){ convert_position= flag; }
};

template<class Tgas, class Tcdm>
  gadget_file<Tgas, Tcdm>::gadget_file() :
  convert_position(false), convert_velocity(true)
{
  v_gas= 0;
  for(int k=0; k<4; k++)
    v_cdm[k]= 0;
  clear_header();
  
  if(sizeof(int4) != 4)
    std::cout << "Error: int4 is not 4byte but " << sizeof(int4) << "byte\n";
  if(sizeof(gadget_header1) != 256)
    std::cout << "Error: size of the header is not 256 bytes\n";

}

template<class Tgas, class Tcdm> void gadget_file<Tgas, Tcdm>::clear_header()
{
  const int jmax = (int)(256/sizeof(int));
  int *p = (int*) &header;
  for(int j=0; j<jmax; j++)
    p[j]= 0;
}

template<class Tgas, class Tcdm>
  void gadget_file<Tgas, Tcdm>::set_cosmology(
	  const double redshift, const double boxsize,
	  const double omega0, const double omega_lambda,
	  const double hubble_param){
  header.redshift= redshift;
  header.time    = 1.0/(1.0+redshift);
  header.boxsize = boxsize;
  header.omega0  = omega0;
  header.omega_lambda = omega_lambda;
  header.hubble_param = hubble_param;
}

template<class Tgas, class Tcdm>
  void gadget_file<Tgas, Tcdm>::set_mass(
	  const double mgas, const double mhalo,
	  const double mdisk, const double mbulge,
	  const double mstars){
  header.mass[0]= mgas;
  header.mass[1]= mhalo;
  header.mass[2]= mdisk;
  header.mass[3]= mbulge;
  header.mass[4]= mstars;
}

template<class Tgas, class Tcdm>
  void gadget_file<Tgas, Tcdm>::set_n_particles(
	  const int ngas, const int nhalo,
	  const int ndisk, const int nbulge,
	  const int nstars){
  header.np[0]= ngas;
  header.np[1]= nhalo;
  header.np[2]= ndisk;
  header.np[3]= nbulge;
  header.np[4]= nstars;
}

template<class Tgas, class Tcdm>
  void gadget_file<Tgas, Tcdm>::set_n_particles_total(
	  const int ngas, const int nhalo,
	  const int ndisk, const int nbulge,
	  const int nstars){
  header.np_total[0]= ngas;
  header.np_total[1]= nhalo;
  header.np_total[2]= ndisk;
  header.np_total[3]= nbulge;
  header.np_total[4]= nstars;
}



template<class Tgas, class Tcdm>
  int gadget_file<Tgas, Tcdm>::fread_seperator(FILE * const fp, int expected)
{
  int4 dummy;
  if(fread(&dummy, sizeof(dummy), 1, fp) != 1)
    return false;
  if(dummy != expected){
    std::cout << "Sperator Error: " << dummy << " not equal to "
	      << expected << "\n";
    return false;
  }
  return true;
}

template<class Tgas, class Tcdm>
  int gadget_file<Tgas, Tcdm>::max(int *a)
{
  int amax=a[0];
  for(int j=1; j<5; j++)
    amax = (amax > a[j] ? amax : a[j]);
  
  return amax;
}

template<class Tgas, class Tcdm>
  int gadget_file<Tgas, Tcdm>::read_header(const char * const filename)
{
  FILE *fp= 0;

  if(!(fp=fopen(filename, "rb"))) {
    char* filename1= (char*) malloc(sizeof(strlen(filename)+5)*sizeof(char));
    sprintf(filename1, "%s.0", filename);
    if(!(fp=fopen(filename1, "rb"))) {
      free(filename1);
      return false;
    }
    free(filename1);
    // std::cout << "can't open file " << filename << "\n";
  }

  fread_seperator(fp, sizeof(gadget_header1));
  fread(&header, sizeof(gadget_header1), 1, fp);
  fread_seperator(fp, sizeof(gadget_header1));
    
  fclose(fp);
  return true;
  
}

template<class Tgas, class Tcdm>
  int gadget_file<Tgas, Tcdm>::read_one(const char filename[])
{
  FILE* fp= 0;
  using namespace std;

  if(!(fp=fopen(filename, "rb"))) {
    //std::cerr << "can't open file " << filename << endl;
    return false;
  }

  int blklen= sizeof(gadget_header1);
  fread_seperator(fp, blklen);
  fread(&header, sizeof(gadget_header1), 1, fp);
  fread_seperator(fp, blklen);
  fclose(fp);
  
  size_t n_reserve[]= {0,0,0,0};
  // reserve memory and fill dummy elements to vectors
  if(v_gas && header.np > 0)
    v_gas->reserve(v_gas->size() + header.np[0]);

  for(int k=0; k<4; k++) {
    if(v_cdm[k] && header.np[k+1] > 0) {
      for(int j=0; j<=k; ++j) {	
	if(v_cdm[j] == v_cdm[k])
	  n_reserve[j] += header.np[k+1];
      }
    }
  }
  
  for(int k=0; k<4; k++) {
    if(v_cdm[k]) {
      v_cdm[k]->reserve(v_cdm[k]->size() + n_reserve[k]);
      //cout << "v" << k << ": " << n_reserve[k] << " reserved\n";
    }
  }
  
  return read1(filename);
}

template<class Tgas, class Tcdm>
  int gadget_file<Tgas, Tcdm>::read(const char * const filename)
{
  FILE* fp= 0;
  bool multiple_files= false;
  using namespace std;

  char* filename1= 0;
  if(!(fp=fopen(filename, "rb"))) {
    //std::cerr << "can't open file " << filename << endl;
    filename1= (char*) malloc(sizeof(char)*(strlen(filename)+5));
    sprintf(filename1, "%s.0", filename);
    if(!(fp=fopen(filename1, "rb"))) {
      std::cerr << "can't open file " << filename1 << " either\n";
      return false;
    }
    multiple_files= true;
  }

  int blklen= sizeof(gadget_header1);
  fread_seperator(fp, blklen);
  fread(&header, sizeof(gadget_header1), 1, fp);
  fread_seperator(fp, blklen);
  fclose(fp);
  
  size_t n_reserve[]= {0,0,0,0};
  // reserve memory and fill dummy elements to vectors
  if(v_gas && header.np_total[0] > 0)
    v_gas->reserve(v_gas->size() + header.np_total[0]);

  for(int k=0; k<4; k++) {
    if(v_cdm[k] && header.np_total[k+1] > 0) {
      for(int j=0; j<=k; ++j) {	
	if(v_cdm[j] == v_cdm[k])
	  n_reserve[j] += header.np_total[k+1];
      }
    }
  }
  
  for(int k=0; k<4; k++) {
    if(v_cdm[k]) {
      v_cdm[k]->reserve(v_cdm[k]->size() + n_reserve[k]);
      //cout << "v" << k << ": " << n_reserve[k] << " reserved\n";
    }
  }


  
  if(multiple_files) {
    int i=0;
    do {
      sprintf(filename1, "%s.%d", filename, i++);
      //cout << filename1 << endl;
    } while(read1(filename1));
    free(filename1);
  }
  else {
    return read1(filename);
  }
  return true;
}

template<class Tgas, class Tcdm>
  int gadget_file<Tgas, Tcdm>::read1(const char * const filename)
{
  //
  // Read one of the files filename.<i>
  //
  float  *pbuf;
  partid_t *pbuf_id;
  FILE   *fp= 0;

  Tgas gas0;
  Tcdm cdm0; // initialized?
  int n_gas_begin=0, n_gas_end=0;
  int n_cdm_begin[4], n_cdm_end[4];

  using namespace std;

  if(!(fp=fopen(filename, "rb")))
    return false;

  //std::cout << "reading " << filename << endl;

  int blklen= sizeof(gadget_header1);
  fread_seperator(fp, blklen);
  fread(&header, sizeof(gadget_header1), 1, fp);
  fread_seperator(fp, blklen);

  int ntot=0, ntot_withmasses= 0;
  for(int k=0; k<5; k++){
    ntot += header.np[k];
    if(header.mass[k]==0)
      ntot_withmasses+= header.np[k];
  }
  
  if(v_gas && header.np[0] > 0) {
                                                      assert(header.np[0] >= 0);
    n_gas_begin= (int) v_gas->size();
    v_gas->insert(v_gas->end(), header.np[0], gas0);
    n_gas_end= (int) v_gas->size();
  }

  for(int k=0; k<4; k++) {
    if(v_cdm[k] && header.np_total[k+1]>0) {
                                                    assert(header.np[k+1] >= 0);
      n_cdm_begin[k]= v_cdm[k]->size();
      v_cdm[k]->insert(v_cdm[k]->end(), header.np[k+1], cdm0);
      n_cdm_end[k]= v_cdm[k]->size();
    }
  }
  

  // Allocate memory for buf
  const int nbuf = 3*max(header.np);
  assert(3*sizeof(float) >= sizeof(partid_t));
  float * const buf = (float*) malloc(sizeof(float)*nbuf);
  if(buf == 0) {
    std::cerr << "Unable to allocate buffer memory: " << nbuf*sizeof(float)
         << " bytes\n";
  }

  blklen= 3*ntot*sizeof(float);
  fread_seperator(fp, blklen);
  
  // Position
  const float x_factor= (convert_position && header.hubble_param > 0.0) ?
                        header.time : 1.0f;
  assert(x_factor > 0.0f && x_factor <= 1.01f);

  if(header.np[0]>0){
    if(v_gas && (*v_gas)[0].get_position()){
      fread(pbuf=buf, sizeof(float), 3*header.np[0], fp);
      for(int n=n_gas_begin; n!=n_gas_end; n++){
	(*v_gas)[n].get_position()[0]= x_factor*pbuf[0];
	(*v_gas)[n].get_position()[1]= x_factor*pbuf[1];
	(*v_gas)[n].get_position()[2]= x_factor*pbuf[2];
	pbuf += 3;
      }
    }
    else
      fseek(fp, 3*sizeof(float)*header.np[0], SEEK_CUR);
  }
  for(int k=0; k<4; k++){
    if(header.np[k+1]>0){
      if(v_cdm[k] && (*v_cdm[k])[0].get_position()){
	fread(pbuf=buf, sizeof(float), 3*header.np[k+1], fp);
	for(int n=n_cdm_begin[k]; n!=n_cdm_end[k]; n++){
	  (*v_cdm[k])[n].get_position()[0]= x_factor*pbuf[0];
	  (*v_cdm[k])[n].get_position()[1]= x_factor*pbuf[1];
	  (*v_cdm[k])[n].get_position()[2]= x_factor*pbuf[2];
	  pbuf += 3;
	}
      }
      else
	fseek(fp, 3*sizeof(float)*header.np[k+1], SEEK_CUR);
 
    }
  }
  fread_seperator(fp, blklen);
  fread_seperator(fp, blklen);

  // Velocity
  const float v_factor= (convert_velocity && header.hubble_param > 0.0) ?
                        sqrt(header.time) : 1.0f;
  //std::cerr << "v_factor " << v_factor << endl;
  assert(v_factor > 0.0f && v_factor <= 1.01f);
  

  if(header.np[0]>0){
    if(v_gas && (*v_gas)[0].get_velocity()){
      fread(pbuf=buf, sizeof(float), 3*header.np[0], fp);
      for(int n=n_gas_begin; n!=n_gas_end; n++){
	(*v_gas)[n].get_velocity()[0]= pbuf[0]*v_factor;
	(*v_gas)[n].get_velocity()[1]= pbuf[1]*v_factor;
	(*v_gas)[n].get_velocity()[2]= pbuf[2]*v_factor;
	pbuf += 3;
      }
    }
    else
      fseek(fp, 3*sizeof(float)*header.np[0], SEEK_CUR);
  } 
 
  for(int k=0; k<4; k++){
    if(header.np[k+1]>0){
      if(v_cdm[k] &&  (*v_cdm[k])[0].get_velocity()){
	fread(pbuf=buf, sizeof(float), 3*header.np[k+1], fp);
	for(int n=n_cdm_begin[k]; n!=n_cdm_end[k]; n++){
	  (*v_cdm[k])[n].get_velocity()[0]= pbuf[0]*v_factor;
	  (*v_cdm[k])[n].get_velocity()[1]= pbuf[1]*v_factor;
	  (*v_cdm[k])[n].get_velocity()[2]= pbuf[2]*v_factor;
	  pbuf += 3;
	}
      }
      else
	fseek(fp, 3*sizeof(float)*header.np[k+1], SEEK_CUR);
    }
  }

  fread_seperator(fp, blklen);
  blklen= ntot*sizeof(partid_t);
  fread_seperator(fp, blklen);

  // ID
  if(header.np[0]>0){
    if(v_gas && (*v_gas)[0].get_id()){
      fread(pbuf_id=(partid_t*)buf, sizeof(partid_t), header.np[0], fp);
      for(int n=n_gas_begin; n!=n_gas_end; n++){
	*(*v_gas)[n].get_id()= *pbuf_id;
	pbuf_id++;
      }
    }
    else
      fseek(fp, sizeof(float)*header.np[0], SEEK_CUR);
  }
 
  for(int k=0; k<4; k++){
    if(header.np[k+1]>0){
      if(v_cdm[k] &&  (*v_cdm[k])[0].get_id()){
	fread(pbuf_id=(partid_t*)buf, sizeof(partid_t), header.np[k+1], fp);
	for(int n=n_cdm_begin[k]; n!=n_cdm_end[k]; n++){
	  *(*v_cdm[k])[n].get_id()= *pbuf_id;
	  pbuf_id++;
	}
      }
      else
	fseek(fp, sizeof(float)*header.np[k+1], SEEK_CUR);
    }
  }
  fread_seperator(fp, blklen);
  blklen= ntot_withmasses * sizeof(float);
  if(ntot_withmasses>0)
    fread_seperator(fp, blklen);

  // Mass
  if(header.np[0]>0){
    if(v_gas && (*v_gas)[0].get_mass()){
      if(header.mass[0] == 0.){
	fread(pbuf=buf, sizeof(float), header.np[0], fp);
	for(int n=n_gas_begin; n!=n_gas_end; n++){
	  *(*v_gas)[n].get_mass()= *pbuf;
	  pbuf++;
	}
      }
      else{
	for(int n=n_gas_begin; n!=n_gas_end; n++)
	  *(*v_gas)[n].get_mass()= header.mass[0];
      }
    }
    else{
      if(header.mass[0] == 0)
	fseek(fp, sizeof(float)*header.np[0], SEEK_CUR);
    }
  }
 
  for(int k=0; k<4; k++){
    if(header.np[k+1]>0){
      if(v_cdm[k] &&  (*v_cdm[k])[0].get_mass()){
	if(header.mass[k+1] == 0.){
	  fread(pbuf=buf, sizeof(float), header.np[k+1], fp);
	  for(int n=n_cdm_begin[k]; n!=n_cdm_end[k]; n++){
	    *(*v_cdm[k])[n].get_mass()= *pbuf;
	    pbuf++;
	  }
	}
	else{
	  for(int n=n_cdm_begin[k]; n!=n_cdm_end[k]; n++)
	    *(*v_cdm[k])[n].get_mass()= header.mass[k+1];
	}          
      }
      else{
	if(header.mass[k+1] == 0.)
	  fseek(fp, sizeof(float)*header.np[k+1], SEEK_CUR);
      }
    }
  }

  if(ntot_withmasses>0)
    fread_seperator(fp, blklen);
    
  blklen= header.np[0]*sizeof(float);
  //std::cout << "eof " << feof(fp) << endl;
  if(header.np[0]>0 && fread_seperator(fp,blklen)){
    // Internal Energy
    if(v_gas && (*v_gas)[0].get_internal_energy()){
      fread(pbuf=buf, sizeof(float), header.np[0], fp);
      for(int n=n_gas_begin; n!=n_gas_end; n++){
	*(*v_gas)[n].get_internal_energy()= *pbuf;
	pbuf++;
      }
    }
    else
      fseek(fp, sizeof(float)*header.np[0], SEEK_CUR);
    
    fread_seperator(fp, blklen);

    fread_seperator(fp, blklen);
    // Density
    if(v_gas && (*v_gas)[0].get_density()){
      fread(pbuf=buf, sizeof(float), header.np[0], fp);
      for(int n=n_gas_begin; n!=n_gas_end; n++){
	*(*v_gas)[n].get_density()= *pbuf;
	pbuf++;
      }
    }
    else
      fseek(fp, sizeof(float)*header.np[0], SEEK_CUR);  
    fread_seperator(fp, blklen);

    // Electron Density (if cooling)
    if(header.flag_cooling){
      fread_seperator(fp, 15);
      if(v_gas && (*v_gas)[0].get_electron_density()){
	fread(pbuf=buf, sizeof(float), header.np[0], fp);
	for(int n=n_gas_begin; n!=n_gas_end; n++){
	  *(*v_gas)[n].get_electron_density()= *pbuf;
	  pbuf++;
	}
      }
      else
	fseek(fp, sizeof(float)*header.np[0], SEEK_CUR);  
      fread_seperator(fp, blklen);
    }

    fread_seperator(fp, blklen);
    // Smoothing Length Hsml
    if(v_gas && (*v_gas)[0].get_smoothing_length()){
      fread(pbuf=buf, sizeof(float), header.np[0], fp);
      for(int n=n_gas_begin; n!=n_gas_end; n++){
	*(*v_gas)[n].get_smoothing_length()= *pbuf;
	pbuf++;
      }
    }
    else
      fseek(fp, sizeof(float)*header.np[0], SEEK_CUR);  
    fread_seperator(fp, blklen);

  }

  fclose(fp);
  free(buf);

  // std::cout << "read() done.\n";

  return true;
}


template<class Tgas, class Tcdm>
  int gadget_file<Tgas, Tcdm>::write(const char * const filename)
{
  FILE *fp;
  int blklen;
  using namespace std;

  // Check that the vectors have enough data
  if(v_gas){
  /*
    if(!((*v_gas)[0].get_position() && (*v_gas)[0].get_velocity()
       && (*v_gas)[0].get_internal_energy() && (*v_gas)[0].get_density()
       &&  (*v_gas)[0].get_electron_density())){
  */ 
    if(!((*v_gas)[0].get_position() && (*v_gas)[0].get_velocity())){
      std::cout << "Need position, velocity for gas particle in order to write\n";
      return false;
    }
    else if(header.mass[0]!=0 && (*v_gas)[0].get_mass()==0){
      std::cout << "Need mass for gas particles when header.mass[0] == 0\n";
      return false;
    }
    header.np[0]= v_gas->size();
    // std::cout << "v_gas has \n" << header.np[0] << " particles\n";
  }
  
  for(int k=0; k<4; k++){
    if(v_cdm[k]){
      if(!((*v_cdm[k])[0].get_position()
	   && (*v_cdm[k])[0].get_velocity())){
	std::cout << "Need position, velocity for type " << k+1
	     << " particles in order to write\n";
	return false;
      }
      else if(header.mass[k+1]!=0 && (*v_cdm[k])[0].get_mass()==0){
	std::cout << "Need mass for type " << k+1 
	     << "particles when header.mass is zero\n";
	return false;
      }
      header.np[k+1]= v_cdm[k]->size();
      // std::cout << "v_cdm Type " << k+1 << " has " << header.np[k+1] <<
      //         " particles" << endl;
    } 
  }
  
  if(header.num_files == 0)
    std::cout << "Warning! num_files is zero\n";

  if(!(fp=fopen(filename,"wb"))){
    // printf("can't open file `%s`\n", filename);
    return 1;
  }

  // printf("writing `%s' ...", filename);

  blklen= sizeof(gadget_header1);
  fwrite(&blklen, sizeof(blklen), 1, fp);
  fwrite(&header, sizeof(gadget_header1), 1, fp);
  fwrite(&blklen, sizeof(blklen), 1, fp);

  int ntot_withmasses=0, ntot=0;
  for(int k=0; k<5; k++){
    ntot += header.np[k];
    if(header.mass[k]==0)
      ntot_withmasses+= header.np[k];
  }

  // position
  blklen= ntot*sizeof(float)*3;
  fwrite(&blklen, sizeof(blklen), 1, fp);
  for(int n=0; n<header.np[0]; n++){
      fwrite((*v_gas)[n].get_position(), sizeof(float), 3, fp);
  }
  for(int k=0; k<4; k++){
    for(int n=0; n<header.np[k+1]; n++)
      fwrite((*v_cdm[k])[n].get_position(), sizeof(float), 3, fp);
  }
 
  fwrite(&blklen, sizeof(blklen), 1, fp);
  
  // velocity
  fwrite(&blklen, sizeof(blklen), 1, fp);
  for(int n=0; n<header.np[0]; n++)
      fwrite((*v_gas)[n].get_velocity(), sizeof(float), 3, fp);

  for(int k=0; k<4; k++){
    for(int n=0; n<header.np[k+1];n++)
      fwrite((*v_cdm[k])[n].get_velocity(), sizeof(float), 3, fp);
  }
  fwrite(&blklen, sizeof(blklen), 1, fp);
    
  // id
  const int default_id= 0;
  blklen= ntot*sizeof(int);
  fwrite(&blklen, sizeof(blklen), 1, fp);

  if(v_gas && (*v_gas)[0].get_id()){
    for(int n=0; n<header.np[0]; n++){
      fwrite((*v_gas)[n].get_id(), sizeof(int), 1, fp);
    }
  }
  else
    for(int n=0; n<header.np[0]; n++){
      fwrite(&default_id, sizeof(int), 1, fp);
    }
  
  for(int k=0; k<4; k++){
    if(v_cdm[k] && (*v_cdm[k])[0].get_id()){
      for(int n=0; n<header.np[k+1]; n++)
	fwrite((*v_cdm[k])[n].get_id(), sizeof(int), 1, fp);
    }
    else
      for(int n=0; n<header.np[k+1]; n++)
	fwrite(&default_id, sizeof(int), 1, fp);
  }
  fwrite(&blklen, sizeof(blklen), 1, fp);

  // mass
  blklen= ntot_withmasses*sizeof(float);
  if(ntot_withmasses>0)
    fwrite(&blklen, sizeof(blklen), 1, fp);

  if(header.mass[0] == 0.){
    for(int n=0; n<header.np[0]; n++)
	fwrite((*v_gas)[n].get_mass(), sizeof(float), 1, fp);
  }
  for(int k=0; k<4; k++){
    if(header.mass[k+1]==0.)
      for(int n=0; n<header.np[k+1]; n++)
	fwrite((*v_cdm[k])[n].get_mass(), sizeof(float), 1, fp);
  }
  
  if(ntot_withmasses>0)
    fwrite(&blklen, sizeof(blklen), 1, fp);
  
  if(header.np[0]>0){
    // Internal Energy
    blklen= header.np[0]*sizeof(float);

    fwrite(&blklen, sizeof(blklen), 1, fp);
    if((*v_gas)[0].get_internal_energy()){
      for(int n=0; n<header.np[0]; n++)
        fwrite((*v_gas)[n].get_internal_energy(), sizeof(float), 1, fp);
    }
    else { // Internal Energy block must exists. Fill with zeros.
      float internal_energy= 0.0f;
      for(int n=0; n<header.np[0]; n++)
        fwrite(&internal_energy, sizeof(float), 1, fp);

    }
    fwrite(&blklen, sizeof(blklen), 1, fp);

    if((*v_gas)[0].get_density()){
      fwrite(&blklen, sizeof(blklen), 1, fp);
      for(int n=0; n<header.np[0];n++)
        fwrite((*v_gas)[n].get_density(), sizeof(float), 1, fp);
      fwrite(&blklen, sizeof(blklen), 1, fp);
    }

    if(header.flag_cooling && (*v_gas)[0].get_electron_density()){
      fwrite(&blklen, sizeof(blklen), 1, fp);
      for(int n=0; n<header.np[0];n++)
	fwrite((*v_gas)[n].get_electron_density(), sizeof(float), 1, fp);
    
      fwrite(&blklen, sizeof(blklen), 1, fp);
    } 
  }
  
  int file_state = !ferror(fp);
  fclose(fp);
  // printf("done\n");
  
  return file_state;
}



template<class Tgas, class Tcdm>
  void gadget_file<Tgas, Tcdm>::print()
{
  using namespace std;
  std::cout << endl;
  std::cout << "Number of Particles       "
       << header.np[0] << " " << header.np[1] << " " << header.np[2] << " "
       << header.np[3] << " " << header.np[4] << " " << header.np[5] << endl;
  std::cout << "Total Number of Particles "
       << header.np_total[0] << " " 
       << header.np_total[1] << " " << header.np_total[2] << " " 
       << header.np_total[3] << " " <<  header.np_total[4] << " "
       << header.np_total[5] << endl;
  std::cout << "Number of Files           " << header.num_files << "\n";
  std::cout << "Mass                      " << header.mass[0] << " "
       << header.mass[1] << " " << header.mass[2] << " "
       << header.mass[3] << " " << header.mass[4] << " "
       << header.mass[5] << endl;

  std::cout << "Time                      " << header.time << endl;
  std::cout << "Red Shift                 " << header.redshift << endl;
  std::cout << "Box Size                  " << header.boxsize << endl;
  std::cout << "Omega Matter              " << header.omega0 << endl;
  std::cout << "Omega Lambda              " << header.omega_lambda << endl;
  std::cout << "Hubble Param              " << header.hubble_param << endl;

  if(header.np[0]>0 && v_gas){
    printf("Type 0: gas\n");
    for(int n=0; n<5; n++){
      if((*v_gas)[0].get_position())
        printf("x=(% 6.1f,% 6.1f,% 6.1f) ",
	        (*v_gas)[n].get_position()[0],
	        (*v_gas)[n].get_position()[1],
	        (*v_gas)[n].get_position()[2]);
      if((*v_gas)[0].get_velocity())
        printf("v=(% 6.1f,% 6.1f,% 6.1f) ",
	        (*v_gas)[n].get_velocity()[0],
	        (*v_gas)[n].get_velocity()[1],
	        (*v_gas)[n].get_velocity()[2]);
      if((*v_gas)[n].get_mass())
        printf("m=% 8.2e ", *(*v_gas)[n].get_mass());
	
      if((*v_gas)[n].get_id())
        printf("id=%7d ", *(*v_gas)[n].get_id());
      printf("\n");
    }
  }
  printf("\n");

  for(int k=0; k<4; k++){
    if(header.np[k+1]>0 && v_cdm[k]){
      printf("Type %d\n", k+1);
      for(int n=0; n<5; n++){
	if((*v_cdm[k])[0].get_position())
	  printf("x=(% 6.1f,% 6.1f,% 6.1f) ",
		 (*v_cdm[k])[n].get_position()[0],
		 (*v_cdm[k])[n].get_position()[1],
		 (*v_cdm[k])[n].get_position()[2]);
	if((*v_cdm[k])[0].get_velocity())
	  printf("v=(% 6.1f,% 6.1f,% 6.1f) ",
		 (*v_cdm[k])[n].get_velocity()[0],
		 (*v_cdm[k])[n].get_velocity()[1],
		 (*v_cdm[k])[n].get_velocity()[2]);
	if((*v_cdm[k])[n].get_mass())
	  printf("m=% 8.2e ", *(*v_cdm[k])[n].get_mass());
	
	if((*v_cdm[k])[n].get_id())
	  printf("id=%7d ", *(*v_cdm[k])[n].get_id());
	printf("\n");
      }
      printf("\n");
    }
  }
  std::cout << endl;
}

// Default particle_data_class


class particle_data_sph_all{
 public:
  float x[3];
  float v[3];
  partid_t  id;
  float m;
  float u;
  float rho;
  float rk;
  float* get_position(){return x;}
  float* get_velocity(){return v;}
  partid_t*  get_id(){return &id;}
  float* get_mass(){return &m;}
  float* get_internal_energy(){return &u;}
  float* get_density(){return &rho;}
  float* get_electron_density(){return 0;}
  float* get_smoothing_length(){return &rk;}
};

class particle_data_cdm_all{
 public:
  float x[3];
  float v[3];
  partid_t   id;
  float m;
  float* get_position(){return x;};
  float* get_velocity(){return v;};
  partid_t*   get_id(){return &id;};
  float* get_mass(){return &m;};
};

#endif
