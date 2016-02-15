#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "parallel.h"
using namespace std;

// **************************************************************
//    ADJACENCY ARRAY REPRESENTATION
// **************************************************************

struct EdgeInfo {
    intE weight;
    intE capacity;
    intE lower;
    intE flow;
};

typedef EdgeInfo edgeInfo;

struct vertex {
#ifndef WEIGHTED
  uintE* inNeighbors, *outNeighbors;
#else
  intE* inNeighbors, *outNeighbors;
#endif
  uintT outDegree;
  uintT inDegree;
  void del() {free(inNeighbors); free(outNeighbors);}
#ifndef WEIGHTED
vertex(uintE* iN, uintE* oN, uintT id, uintT od) 
#else
vertex(intE* iN, intE* oN, uintT id, uintT od) 
#endif
: inNeighbors(iN), outNeighbors(oN), inDegree(id), outDegree(od) {}
#ifndef WEIGHTED
  uintE* getInNeighbors () { return inNeighbors; }
  uintE* getOutNeighbors () { return outNeighbors; }
  uintE getInNeighbor(uintT j) { return inNeighbors[j]; }
  uintE getOutNeighbor(uintT j) { return outNeighbors[j]; }
  void setInNeighbors(uintE* _i) { inNeighbors = _i; }
  void setOutNeighbors(uintE* _i) { outNeighbors = _i; }
#else 
  intE* getInNeighbors () { return inNeighbors; }
  intE* getOutNeighbors () { return outNeighbors; }
  intE getInNeighbor(uintT j) { return inNeighbors[2*j]; }
  intE getOutNeighbor(uintT j) { return outNeighbors[2*j]; }
  intE getInInfo(uintT j) { return inNeighbors[2*j+1]; }
  intE getOutInfo(uintT j) { return outNeighbors[2*j+1]; }
  void setInNeighbors(intE* _i) { inNeighbors = _i; }
  void setOutNeighbors(intE* _i) { outNeighbors = _i; }
#endif
  uintT getInDegree() { return inDegree; }
  uintT getOutDegree() { return outDegree; }
  void setInDegree(uintT _d) { inDegree = _d; }
  void setOutDegree(uintT _d) { outDegree = _d; }
  void flipEdges() { swap(inNeighbors,outNeighbors); swap(inDegree,outDegree); }
};

struct graph {
  vertex *V;
  edgeInfo *E;
  long n;
  long m;
  double* p;
  intE* supply;
#ifndef WEIGHTED
  uintE* allocatedInplace, * inEdges;
#else
  intE* allocatedInplace, * inEdges;
#endif
  uintE* flags;
  bool transposed;
  graph(vertex* VV, long nn, long mm) 
  : V(VV), n(nn), m(mm), allocatedInplace(NULL), flags(NULL), transposed(false) {}
#ifndef WEIGHTED
  graph(vertex* VV, long nn, long mm, uintE* ai, uintE* _inEdges = NULL) 
#else
  graph(vertex* VV, long nn, long mm, intE* ai, intE* _inEdges = NULL, edgeInfo* ei = NULL) 
#endif
  : V(VV), n(nn), m(mm), allocatedInplace(ai), inEdges(_inEdges), flags(NULL), transposed(false), E(ei) {}
  void del() {
    if (flags != NULL) free(flags);
    if (allocatedInplace == NULL) 
      for (long i=0; i < n; i++) V[i].del();
    else free(allocatedInplace);
    free(V);
    free(E);
    free(p);
    free(supply);
    if(inEdges != NULL) free(inEdges);
  }
  void transpose() {
    parallel_for(long i=0;i<n;i++) {
      V[i].flipEdges();
    }
    transposed = !transposed;
  }
  edgeInfo& getInInfo(uintT v, uintT inEdgeIndex) {
      return E[V[v].getInInfo(inEdgeIndex)];
  }
  edgeInfo& getOutInfo(uintT v, uintT outEdgeIndex) {
      return E[V[v].getOutInfo(outEdgeIndex)];
  }
};
#endif