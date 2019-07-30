// -*- outline-regexp: "/\\*-+";  -*-
/*-  Identification and Changes  */

/*
  cube.cxx -- Written by Randall Gray 
  Initial coding: 
  Date: 2007.09.06
  Location: odin.valhalla.asgard:/usr/home/gray/Projects/Ningaloo/cube.cxx

  Compile with: g++ -ggdb -DDEBUGGING -Wall -o cube cube.cxx
  History:


  2019-07-30-07:21:00 -- minor explanatory note 

  This code is, at least as far as I know, still novel. There is very
  little work on multi-contaminant sources of mortality.  While this
  approach is couched in terms of contaminants, it would apply equally
  well to other aquatic/marine conditions such as hypoxia or
  overheating.  It necessarily assumes the linear independence of the
  sources of mortality in its calculation, which is likely to be a bit
  of a stretch when it comes to things like hypoxia.  This might be
  adjusted by using a non-linear function to map "cause" onto "effect"
  -- rather than use dissolved oxygen, we might use a function that maps 
  the oxygen levels to something appropriate.

  The way this approach works is that we first construct a unit
  hypercube whose axes take values in [0,1] and identify the the
  complement of its volume in the unit hypercube with the number of
  individuals in a population. To begin with, the cube has zero
  volume, so the population is unaffected. When an individual dies of
  natural causes (predation, age, ...) an appropriate ordinate in the
  point which defines the cube is incremented so that the product of
  the reference population and the size of the complement of the cube
  in the unit cube is the new population level.


  $Log: cube.cxx,v $
  Revision 1.9  2009/03/18 03:54:18  gray
  Code writers committed for psychiatric assessment w.r.t. stasis/migration

  Revision 1.8  2009/02/20 01:27:40  gray
  Whopping big patch bomb

  Revision 1.7  2008/08/29 04:11:21  gray
  More buggering with indenting!

  Revision 1.6  2008/08/26 08:21:39  gray
  Massive mods.  Lotsa space fixing and such. Merged Beths stuff in

  Revision 1.5  2007/10/24 00:21:03  gray
  Minor emacs tweaks -- no code changes yet

  Revision 1.4  2007/10/12 02:58:12  gray
  muchous state save/restore stuff

  Revision 1.3  2007/10/09 01:29:04  gray
  Reworked init_contaminants to do it incrementally.  Will need to test with a taxon change or some such thing

  Revision 1.2  2007/10/08 03:12:34  gray
  Why dont you come with me while we poison the fishes in the park, la da de da da de da

  Revision 1.1  2007/09/21 00:44:12  gray
  poisoning support


*/

/*-  Copyright  */

/*
  (C) 2007 CSIRO Australia
  All rights reserved
*/

/*-  Discussion  */

/*-  Configuration stuff  */

#ifndef __cube_cxx
#define __cube_cxx
#endif

/*-  Included files  */

#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#include "cube.hxx"
#include "memchk.h"
#include "memchk.h"

/*-  Local variables, constants, and defines  */

/*-  Code  */

// Index 0 is the natural mortality thing.  Indices 1...  correspond to the *other* bits


/*-- Serialisation */
/*--- GetState(int *sz) -- */

void *Cube::GetState(int *sz) {
	double *d = (double*)Calloc(2+n, sizeof(double));
	if (!d) abort();
	d[0] = (double)n;
	d[1] = value;
	for (int i = 0; i < n; i++) {
		d[i+2] = v[i];
	}

	assert(sz);
	*sz = (n+2)*sizeof(double);
	return d;
}

/*--- SetState(void *data, int sz) -- */

void Cube::SetState(void *data, int sz) {
	double *d = (double*)data;
	assert(data);

	n = (int)(d[0]);
	assert((unsigned)sz == sizeof(double)*(n+2));

	value = d[1];
	if (v) Free(v);
	v = (double *)Calloc(n, sizeof(double));
	if (!v) abort();

	memcpy(v, d+2, (sizeof(double)*n));
}

// This class gets handed the new load at the end of each time step, and it updates the contact cube 

/*-- Constructors & Destructor */
/*--- Cube(int N) -- make a pristine (empty) cube */

Cube::Cube(int N) { // 
	n = N;
	value = 1.0;
	v = (double *)Calloc(n, sizeof(double));
	if (!v) abort();
};

/*--- Cube(int N, double val) -- make it and set the value associated with each axis */
/* This would not usually be called, unless  we were modelling a known population
   with a known, preexisting mortality due to contaminants  */

Cube::Cube(int N, double val) {
	n = N;
	value = val;
	v = (double *)Calloc(n, sizeof(double));
	if (!v) abort();
};

/*--- ~Cube() -- */

Cube::~Cube() {
	if (v) Free(v);
};





/*-- proportion_of_box(double *v, int dim) -- */

double Cube::proportion_of_box(double *v, int dim) {
	double prod = 1.0;
	int i = 0;

	assert(dim >= 0 && dim <= n);

	for (i = 0; i < dim; i++) {
		if (v[i] > 1 || v[i] < 0) abort();
		prod *= (1.0 - v[i]);
	}
	return prod;
}

/*-- add_dimension() -- */

int Cube::add_dimension() {
	n = n+1;

	v = (double *)Realloc(v, sizeof(*v) * (n+1));
	if (!v) return 0;

	v[n-1] = 0.0;
	return 1;
}

/*-- service calls: getting data and changing data */

/*--- Value() -- */

double Cube::Value() {
	return ceil(value * proportion_of_box(v, n)); 
};

/*--- LValue() -- */

double Cube::LValue() {
	return value * proportion_of_box(v, n); 
};

/*--- AdjustN(double K, int I) -- */

double Cube::AdjustN(double K, int I) { // Removes a number against axis I
	double Q = 1.0;
	int i = 0;
	double cv = LValue();

	if (!v || !value) return 0.0; 

	Q = 1.0;  
	for (i = 0; i < n; i++) {
		if (i == I) continue;
		if (v[i] > 1 || v[i] < 0) abort();

		Q *= (1.0 - v[i]);
	} // Q = Surviorship w.r.t. other axes

	if (Q > 1) {
		abort();
		Q = 1;
	}

	if (K > value*Q) K = value*Q;

	if (Q > 0 && cv > 0) {
		double d = v[I] + K / (value * Q);
		if (d < 1.0) v[I] = d;
		else v[I] = 1.0;
	}


	K = cv - LValue();
	K = floor(cv - Value());

	if (K < 0) return 0;
	else return K;
}

#if 0

/*--- AdjustLevel(double level, int I) -- */

double Cube::AdjustLevel(double level, int I) { // sets a level if "level" is greater than the previous one
	double Q = 1.0;
	int i = 0;
	double cv = Value();
	
	if (level < v[i]) return 0.0; // its already been higher
	else v[i] = level;

	K = Max(0.0, cv - Value()); // number actually decreased

	return K;
}
#endif

/*--- AdjustLevels(double *level, int base, int n) */
// Adjusts the indicated levels [base ... base+n] by the proportion of the remaining range

double Cube::AdjustLevels(double *level, int base, int n) {
	int i = 0;
	double K;

	for (i = 0; i < n; i++) {
		v[i+base] = v[i+base] + level[i] * (1.0 - v[i+base]);
	}

	K = Value();
	if (K < 0) return 0;
	else return K;
}



/*-  The End  */

