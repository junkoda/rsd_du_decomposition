//
// Friends of Friends halo finder
//
// This code is based on a serial FoF halo finder by
// University of Washington N-BODY SHOP
// http://www-hpcc.astro.washington.edu/tools/fof.html
//


/*
						FOF v1.1

			A Group Finder for N-body Simulations

					October 26, 1994
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "particle.h"
#include "msg.h"
#include "fof.h"

static const int nBucket= 16;

#define ROOT	   1
#define LOWER(i)   (i<<1)
#define UPPER(i)   ((i<<1)+1)
#define PARENT(i)  (i>>1)
#define SIBLING(i) ((i&1)?i-1:i+1)
#define SETNEXT(i) { while (i&1) i=i>>1; ++i; }

#define LEFT  2
#define RIGHT 1
#define DUAL  3
#define GLOBAL 4

typedef struct bndBound {
  float fMin[3];
  float fMax[3];
} BND;

typedef struct kdNode {
  float fSplit;
  BND bnd;
  int iDim;
  int pLower;
  int pUpper;
} KDN;

typedef struct kdContext {
  int nBucket;
  int nActive;
  float fPeriod[3];
  int nLevels;
  int nNodes;
  int nSplit;
  Particle* p;
  int* iGroup;
  KDN *kdNodes;
  int nGroup;
} *KD;

#define INTERSECT(c,cp,fBall2,lx,ly,lz,x,y,z,sx,sy,sz)\
{\
  float dx,dy,dz,dx1,dy1,dz1,fDist2,fMax2;	\
  dx = c[cp].bnd.fMin[0]-x;			\
  dx1 = x-c[cp].bnd.fMax[0];			\
  dy = c[cp].bnd.fMin[1]-y;			\
  dy1 = y-c[cp].bnd.fMax[1];			\
  dz = c[cp].bnd.fMin[2]-z;			\
  dz1 = z-c[cp].bnd.fMax[2];			\
  if (dx > 0.0) {				\
    if (dx1+lx < dx) {				\
      dx1 += lx;				\
      dx -= lx;					\
      sx = x+lx;				\
      fDist2 = dx1*dx1;				\
      fMax2 = dx*dx;				\
    }						\
    else {					\
      sx = x;					\
      fDist2 = dx*dx;				\
      fMax2 = dx1*dx1;				\
    }							\
    if (fDist2 > fBall2) goto GetNextCell;		\
  }							\
  else if (dx1 > 0.0) {					\
    if (dx+lx < dx1) {					\
      dx += lx;						\
      dx1 -= lx;					\
      sx = x-lx;					\
      fDist2 = dx*dx;					\
      fMax2 = dx1*dx1;					\
    }							\
    else {						\
      sx = x;						\
      fDist2 = dx1*dx1;					\
      fMax2 = dx*dx;					\
    }							\
    if (fDist2 > fBall2) goto GetNextCell;		\
  }							\
  else {						\
    sx = x;						\
    fDist2 = 0.0;					\
    if (dx < dx1) fMax2 = dx*dx;			\
    else fMax2 = dx1*dx1;				\
  }							\
  if (dy > 0.0) {					\
    if (dy1+ly < dy) {					\
      dy1 += ly;					\
      dy -= ly;						\
      sy = y+ly;					\
      fDist2 += dy1*dy1;				\
      fMax2 += dy*dy;					\
    }							\
    else {						\
      sy = y;						\
      fDist2 += dy*dy;					\
      fMax2 += dy1*dy1;					\
    }							\
    if (fDist2 > fBall2) goto GetNextCell;		\
  }							\
  else if (dy1 > 0.0) {					\
    if (dy+ly < dy1) {					\
      dy += ly;						\
      dy1 -= ly;					\
      sy = y-ly;					\
      fDist2 += dy*dy;					\
      fMax2 += dy1*dy1;					\
    }							\
    else {						\
      sy = y;						\
      fDist2 += dy1*dy1;				\
      fMax2 += dy*dy;					\
    }							\
    if (fDist2 > fBall2) goto GetNextCell;		\
  }							\
  else {						\
    sy = y;						\
    if (dy < dy1) fMax2 += dy*dy;			\
    else fMax2 += dy1*dy1;				\
  }							\
  if (dz > 0.0) {					\
    if (dz1+lz < dz) {					\
      dz1 += lz;					\
      dz -= lz;						\
      sz = z+lz;					\
      fDist2 += dz1*dz1;				\
      fMax2 += dz*dz;					\
    }							\
    else {						\
      sz = z;						\
      fDist2 += dz*dz;					\
      fMax2 += dz1*dz1;					\
    }							\
    if (fDist2 > fBall2) goto GetNextCell;		\
  }							\
  else if (dz1 > 0.0) {					\
    if (dz+lz < dz1) {					\
      dz += lz;						\
      dz1 -= lz;					\
      sz = z-lz;					\
      fDist2 += dz*dz;					\
      fMax2 += dz1*dz1;					\
    }							\
    else {						\
      sz = z;						\
      fDist2 += dz1*dz1;				\
      fMax2 += dz*dz;					\
    }							\
    if (fDist2 > fBall2) goto GetNextCell;		\
  }							\
  else {						\
    sz = z;						\
    if (dz < dz1) fMax2 += dz*dz;			\
    else fMax2 += dz1*dz1;				\
  }							\
  if (fMax2 < fBall2) goto ContainedCell;		\
}

void kdTime(KD,int *,int *);
void kdBuildTree(KD);
int kdFoF(KD,float);
void compute_halo_xv();

//
// Added part for parallel version
//


static float BoxSize, HalfBoxSize;

static KDN* KdNodes;
static int  NKdNodeAlloc;
static int *igrp, *Fifo, *Map;
static HaloInfo* Halo;
static int NHalo, NHaloAlloc;
static void* Buf; static size_t BufSize;



static inline void halo_particle(HaloInfo* const h, Particle const * const p)
{
  if(NHalo >= NHaloAlloc)
    msg_abort("Not enough space for halo\n");

  if(h->nfof == 0) {
    for(int i=0; i<3; i++)
      h->x0[i]= p->x[i];
  }
  else {
    for(int i=0; i<3; i++) {
      float dx= p->x[i] - h->x0[i];
      if(dx > HalfBoxSize)
	dx -= BoxSize;
      else if(dx < -HalfBoxSize)
	dx += BoxSize;

      h->dx_sum[i] += dx;
    }
  }

  for(int i=0; i<3; i++)
    h->v_sum[i] += p->v[i];

  h->nfof++;
}
  
static inline void halo_clear(HaloInfo* const h)
{
  h->nfof= 0;

  for(int i=0; i<3; i++) {
    h->x0[i]= 0.0f;
    h->dx_sum[i]= 0.0f;
    h->v_sum[i]= 0.0f;
  }
}

void kdSelect(KD kd,int d,int k,int l,int r)
{
  Particle* p = kd->p;
  while(r > l) {
    double v = p[k].x[d];
    Particle t = p[r];
    p[r] = p[k];
    p[k] = t;
    int i = l - 1;
    int j = r;
    while(1) {
      while(i < j) if (p[++i].x[d] >= v) break;
      while(i < j) if (p[--j].x[d] <= v) break;
      t = p[i];
      p[i] = p[j];
      p[j] = t;
      if (j <= i) break;
    }
    p[j] = p[i];
    p[i] = p[r];
    p[r] = t;
    if (i >= k) r = i - 1;
    if (i <= k) l = i + 1;
  }
}


void kdCombine(KDN const * const p1, KDN const * const p2, KDN * const pOut)
{
  // Combine the bounds.
  for (int j=0; j<3; j++) {
    if (p2->bnd.fMin[j] < p1->bnd.fMin[j])
      pOut->bnd.fMin[j] = p2->bnd.fMin[j];
    else
      pOut->bnd.fMin[j] = p1->bnd.fMin[j];
    if (p2->bnd.fMax[j] > p1->bnd.fMax[j])
      pOut->bnd.fMax[j] = p2->bnd.fMax[j];
    else
      pOut->bnd.fMax[j] = p1->bnd.fMax[j];
  }
}


void kdUpPass(KD kd, int iCell)
{
  int l, u, pj,j;

  KDN* c = kd->kdNodes;
  if (c[iCell].iDim != -1) {
    l = LOWER(iCell);
    u = UPPER(iCell);
    kdUpPass(kd, l);
    kdUpPass(kd, u);
    kdCombine(&c[l], &c[u], &c[iCell]);
  }
  else {
    l = c[iCell].pLower;
    u = c[iCell].pUpper;
    for (j=0;j<3;++j) {
      c[iCell].bnd.fMin[j] = kd->p[u].x[j];
      c[iCell].bnd.fMax[j] = kd->p[u].x[j];
    }
    for (pj=l;pj<u;++pj) {
      for (j=0;j<3;++j) {
	if (kd->p[pj].x[j] < c[iCell].bnd.fMin[j])
	  c[iCell].bnd.fMin[j] = kd->p[pj].x[j];
	if (kd->p[pj].x[j] > c[iCell].bnd.fMax[j])
	  c[iCell].bnd.fMax[j] = kd->p[pj].x[j];
      }
    }
  }
}

void kdBuildTree(KD kd)
{
  BND bnd;
  
  int n = kd->nActive;
  kd->nLevels = 1;
  int l = 1;
  while(n > kd->nBucket) {
    n = n>>1;
    l = l<<1;
    ++kd->nLevels;
  }
  kd->nSplit = l;
  kd->nNodes = l<<1;
  kd->kdNodes= KdNodes;
  assert(NKdNodeAlloc >= kd->nNodes);
  assert(kd->kdNodes != NULL);

  // Calculate Bounds.
  for (int j=0; j<3; ++j) {
    bnd.fMin[j] = kd->p[0].x[j];
    bnd.fMax[j] = kd->p[0].x[j];
  }
  for (int i=1; i<kd->nActive; ++i) {
    for (int j=0; j<3; ++j) {
      if (bnd.fMin[j] > kd->p[i].x[j]) 
	bnd.fMin[j] = kd->p[i].x[j];
      else if (bnd.fMax[j] < kd->p[i].x[j])
	bnd.fMax[j] = kd->p[i].x[j];
    }
  }

  // Set up ROOT node

  KDN* c = kd->kdNodes;
  c[ROOT].pLower = 0;
  c[ROOT].pUpper = kd->nActive-1;
  c[ROOT].bnd = bnd;
  int i = ROOT;
  while(1) {
    assert(c[i].pUpper - c[i].pLower + 1 > 0);
    if (i < kd->nSplit && (c[i].pUpper - c[i].pLower) > 0) {
      int d = 0;
      for (int j=1; j<3; ++j) {
	if (c[i].bnd.fMax[j]-c[i].bnd.fMin[j] > 
	    c[i].bnd.fMax[d]-c[i].bnd.fMin[d]) d = j;
      }
      c[i].iDim = d;
      
      int m = (c[i].pLower + c[i].pUpper)/2;
      kdSelect(kd,d,m,c[i].pLower,c[i].pUpper);
      
      c[i].fSplit = kd->p[m].x[d];
      c[LOWER(i)].bnd = c[i].bnd;
      c[LOWER(i)].bnd.fMax[d] = c[i].fSplit;
      c[LOWER(i)].pLower = c[i].pLower;
      c[LOWER(i)].pUpper = m;
      c[UPPER(i)].bnd = c[i].bnd;
      c[UPPER(i)].bnd.fMin[d] = c[i].fSplit;
      c[UPPER(i)].pLower = m+1;
      c[UPPER(i)].pUpper = c[i].pUpper;
      int diff = (m-c[i].pLower+1)-(c[i].pUpper-m);
      assert(diff == 0 || diff == 1);
      i = LOWER(i);
    }
    else {
      c[i].iDim = -1;
      SETNEXT(i);
      if (i == ROOT) break;
    }
  }
  kdUpPass(kd,ROOT);
}


int kdFoF(KD kd,float fEps)
{
  float dx,dy,dz,x,y,z,sx,sy,sz,fDist2;
  
  Particle* const p = kd->p;
  KDN* const c= kd->kdNodes;
  const float lx = kd->fPeriod[0];
  const float ly = kd->fPeriod[1];
  const float lz = kd->fPeriod[2];
  const float fEps2 = fEps*fEps;

  const int nFifo = kd->nActive;


  int iHead= 0;
  int iTail= 0;
  int iGroup= 0;
  
  assert(Map);

  NHalo= 0;

  assert(igrp);
  for (int pn=0; pn<kd->nActive; ++pn)
    igrp[pn]= 0;

  for (int pn=0; pn<kd->nActive; ++pn) {
    if (igrp[pn]) continue;

    ++iGroup;

    halo_clear(&Halo[NHalo]);

    //
    // Mark it and add to the do-fifo.
    //
    igrp[pn]= iGroup; halo_particle(&Halo[NHalo], &p[pn]);
    Fifo[iTail++] = pn;
    if (iTail == nFifo) iTail = 0;

    while (iHead != iTail) {
      int pi = Fifo[iHead++];
      if (iHead == nFifo) iHead = 0;

      //
      // Now do an fEps-Ball Gather!
      //
      x = p[pi].x[0];
      y = p[pi].x[1];
      z = p[pi].x[2];
      int cp = ROOT;
      while (1) {
	INTERSECT(c,cp,fEps2,lx,ly,lz,x,y,z,sx,sy,sz);
	//
	// We have an intersection to test.
	//
	if (c[cp].iDim >= 0) {
	  cp = LOWER(cp);
	  continue;
	}
	else {
	  for (int pj=c[cp].pLower; pj<=c[cp].pUpper; ++pj) {
	    if (igrp[pj]) continue;
	    dx = sx - p[pj].x[0];
	    dy = sy - p[pj].x[1];
	    dz = sz - p[pj].x[2];
	    fDist2 = dx*dx + dy*dy + dz*dz;
	    if (fDist2 < fEps2) {
	      //
	      // Mark it and add to the do-fifo.
	      //
	      igrp[pj]= iGroup; halo_particle(&Halo[NHalo], &p[pj]);
	      Fifo[iTail++] = pj;
	      if (iTail == nFifo) iTail = 0;
	    }
	  }
	  SETNEXT(cp);
	  if (cp == ROOT) break;
	  continue;
	}
      ContainedCell:
	for (int pj=c[cp].pLower; pj<=c[cp].pUpper; ++pj) {
	  if (igrp[pj]) continue;
	  //
	  // Mark it and add to the do-fifo.
	  //
	  igrp[pj] = iGroup; halo_particle(&Halo[NHalo], &p[pj]);
	  Fifo[iTail++] = pj;
	  if (iTail == nFifo) iTail = 0;
	}
      GetNextCell:
	SETNEXT(cp);
	if (cp == ROOT) break;
      }
    }

    Map[iGroup]= NHalo;
    NHalo++;
  }

  // remap igrp
  for (int pn=0; pn<kd->nActive; ++pn)
    igrp[pn]= Map[igrp[pn]];

  kd->nGroup = NHalo; // *** +plus1 ? well not used anymore...anyway
  

  return NHalo;
}

//
// Main interface called from main.c
//

static int n_kd_nodes(int n)
{
  int l = 1;
  while(n > nBucket) {
    n= n >> 1;
    l= l << 1;
  }
  return l<<1;
}

size_t fof_calc_memory(const int np_alloc)
{
  int nNodes= n_kd_nodes(np_alloc);
  size_t size= sizeof(KDN)*nNodes;

  size += sizeof(int)*np_alloc*3; // igrp, Map, Fifo

  int n_halo_alloc= np_alloc;
  size += sizeof(HaloInfo)*n_halo_alloc;

  return size;
}

HaloInfo* fof_init(const int np_alloc)
{
  size_t bytes= 0;

  size_t mem_size= fof_calc_memory(np_alloc);
  void* mem= malloc(mem_size); assert(mem);

  msg_printf(msg_info, "%d Mbytes allocated for FOF halo finding for %d haloes.\n",
	     mem_size/(1024*1024), np_alloc);

  
  Buf= mem;
  BufSize= mem_size;

  NHaloAlloc= np_alloc; // depends on collapsed fraction (set to 0.2)
  Halo= (HaloInfo*) mem; mem= Halo + NHaloAlloc;
  bytes += sizeof(HaloInfo)*NHaloAlloc;

  int nNodes= n_kd_nodes(np_alloc);
  KdNodes= (KDN*) mem; mem= KdNodes + nNodes;
  NKdNodeAlloc= nNodes;  
  bytes += sizeof(KDN)*nNodes;
  
  igrp= (int*) mem;        
  Fifo= igrp + np_alloc;
  Map= Fifo + np_alloc;

  mem= Map + np_alloc; bytes += sizeof(int)*3*np_alloc;


  msg_printf(msg_info, "Going to use %d Mbytes allocated for FOF halo finding.\n",
	     bytes/(1024*1024));

  assert(bytes <= mem_size);
  return Halo;
}

size_t fof_find_halos(Particles* particles, const float linking_param)
{
  BoxSize= particles->boxsize; assert(BoxSize > 0.0f);
  HalfBoxSize= 0.5f*BoxSize;

  const float ll= linking_param*BoxSize/
                  pow((double) particles->np_total, 1.0/3.0);

  msg_printf(msg_verbose, "FOF halo finding, ll= %e...\n", ll);

  struct kdContext kdcontext;
  KD kd= &kdcontext;
  kd->nBucket= nBucket;
  kd->p= particles->p;
  kd->iGroup= igrp;
  kd->kdNodes= 0;
  kd->nActive= particles->np_local;
  for (int j=0; j<3; j++) 
    kd->fPeriod[j] = BoxSize;

  kdBuildTree(kd);
  msg_printf(msg_verbose, "KD tree built.\n");
  int nGroup= kdFoF(kd, ll);
  msg_printf(msg_debug, "Number of haloes: %d\n", nGroup);

  compute_halo_xv();

  return nGroup;
}

void fof_assign_halo_velocity(Particles* const particles)
{
  Particle* const p= particles->p;
  
  const size_t n= particles->np_local;
  for(size_t i=0; i<n; ++i) {
    assert(0 <= igrp[i] && igrp[i] < NHalo);
    
    p[i].v[0] = Halo[igrp[i]].v_sum[0];
    p[i].v[1] = Halo[igrp[i]].v_sum[1];
    p[i].v[2] = Halo[igrp[i]].v_sum[2];
  }
}

void fof_write_halos(char filename[])
{
  int nhalo_written= 0;

  FILE* fp= fopen(filename, "w");
  if(fp == 0)
    msg_abort("Error: Unable to write halo to file %s\n", filename);
    
  fclose(fp);

  msg_printf(msg_info, "%d halos written to %s\n", nhalo_written, filename);
}
  
void compute_halo_xv()
{
  HaloInfo* const h= Halo;
  
  for(int i=0; i<NHalo; ++i) {
    assert(h[i].nfof > 0);
    
    for(int k=0; k<3; k++) {
      float x= h[i].x0[k] + h[i].dx_sum[k]/h[i].nfof;
      if(x < 0.0f) x += BoxSize;
      if(x >= BoxSize) x -= BoxSize; 
      
      h[i].x0[k]= x;
      
      h[i].v_sum[k]= h[i].v_sum[k]/h[i].nfof;
    }    
  }
}
